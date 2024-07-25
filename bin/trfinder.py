#!/usr/bin/env python

import argparse
import os
import re
import sys

import gzip

def parse_args(args=None):
    Description = "Identify terminal repeats in FASTA sequences."
    Epilog = "Example usage: python trfinder.py -i sequences.fasta.gz -o filtered_viruses.fna"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input FASTA (gzipped) to be searched for terminal repeats.",
    )
    parser.add_argument(
        "--min_dtr",
        help="Minimum DTR length to be considered a direct terminal repeat.",
        default=20,
        type=int
    )
    parser.add_argument(
        "--min_itr",
        help="Minimum ITR length to be considered an inverse terminal repeat.",
        default=20,
        type=int
    )
    parser.add_argument(
        "--max_itr",
        help="Maximum ITR length to be considered an inverse terminal repeat.",
        default=1000,
        type=int
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="Prefix to be used in renaming contigs and naming output files.",
    )
    return parser.parse_args(args)


def parse_fasta(handle):
    try:
        id = next(handle).split()[0][1:]
    except:
        return
    seq = ""
    for line in handle:
        if line[0]==">":
            yield id, seq
            id = line.split()[0][1:]
            seq = ""
        else:
            seq += line.rstrip().upper()
    yield id, seq


def fetch_dtr(fullseq, min_length):
    startseq = fullseq[:min_length]
    # Find the first occurrence of startseq in the second half of fullseq
    startpos = fullseq.find(startseq, len(fullseq) // 2)
    while startpos != -1:
        endseq = fullseq[startpos:]
        # Check if the match extends to the contig end
        if fullseq.startswith(endseq):
            return endseq
        # Find the next occurrence of startseq in the second half
        startpos = fullseq.find(startseq, startpos + 1)
    return None

def reverse_complement(seq):
    trans = str.maketrans("ACTG", "TGAC")
    return seq[::-1].translate(trans)

def fetch_itr(seq, min_len, max_len):
    rev = reverse_complement(seq)
    # see if minimal substring occurs at end
    if seq[:min_len] == rev[:min_len]:
        # extend to maximum substring, up to <max_len>
        i = min_len + 1
        while seq[:i] == rev[:i] and i <= max_len:
            i += 1
        return seq[: i - 1]
    # no match
    else:
        return None

def fetch_trs(input, prefix, min_dtr, min_itr, max_itr):
    trs = []
    try:
        if not os.path.exists(input):
            return trs
        with gzip.open(input, "rt") as infile:
            handle = parse_fasta(infile)
            for contig_index, contig in enumerate(handle):
                contig_id, contig_seq = contig
                dtr = fetch_dtr(contig_seq, min_dtr)
                itr = fetch_itr(contig_seq, min_itr, max_itr)
                if not dtr and not itr:
                    continue
                tr = {
                    "sample_id":prefix,
                    "contig_id":prefix+"_"+str(contig_index+1),
                    "contig_name":contig_id,
                    "contig_seq":contig_seq,
                    "contig_len":len(contig_seq),
                }
                if dtr:
                    tr.update({"tr_type":"DTR", "tr_seq":dtr})
                elif itr:
                    tr.update({"tr_type":"ITR", "tr_seq":itr})
                tr_nt_counts = [tr["tr_seq"].count(_) for _ in list("ACGT")]
                tr["tr_len"] = len(tr["tr_seq"])
                tr["tr_nt_acgt_count"] = sum(tr_nt_counts)
                tr["tr_nt_max_freq"] = round(100.0*max(tr_nt_counts)/sum(tr_nt_counts),2) if sum(tr_nt_counts) > 0 else None
                tr["tr_nt_n_count"] = (tr["tr_len"]-sum(tr_nt_counts))
                trs.append(tr)
    except:
        pass
    return trs

def write_batch(sample, results):
    with open(f"./{sample}_trfinder.tsv", "w") as tsv_outfile:
        with open(f"./{sample}_trfinder.fasta", "w") as fasta_outfile:
            fields = ['contig_id',
                    'contig_name',
                    # 'contig_seq',
                    'contig_len',
                    'tr_type',
                    'tr_seq',
                    'tr_len',
                    'tr_nt_acgt_count',
                    'tr_nt_n_count',
                    'tr_nt_max_freq'
                    ]
            tsv_outfile.write("\t".join(fields)+"\n")
            for _ in results:
                tsv_outfile.write("\t".join([str(_[f]) for f in fields])+"\n")
                fasta_outfile.write(">" + str(_["contig_name"]) +"\n" + str(_["contig_seq"]) + "\n")

def trfinder(
    input,
    prefix,
    min_dtr,
    min_itr,
    max_itr,
):
    results = fetch_trs(
        input,
        prefix,
        min_dtr,
        min_itr,
        max_itr
    )
    write_batch(prefix, results)


def main(args=None):
    args = parse_args(args)
    trfinder(
        args.input,
        args.prefix,
        args.min_dtr,
        args.min_itr,
        args.max_itr,
    )

if __name__ == "__main__":
    sys.exit(main())
