#!/usr/bin/env python

import argparse
import sys
import gzip

import Bio.SeqIO

def parse_args(args=None):
    Description = "Calculate nucleotide stats (GC %, CDS density, N count) on input FASTA."
    Epilog = "Example usage: python nucleotid_stats.py -i sequences.fasta.gz -o nucleotide_stat.tsv"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input FASTA (gzipped) to be searched for terminal repeats.",
    )
    parser.add_argument(
        "-p",
        "--proteins",
        help="Path to input FAA file (gzipped) containing protein sequences associated with input FASTA.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output TSV file containing nucleotide statistics.",
    )
    return parser.parse_args(args)


file = sys.argv[1]

def calc_gc(x):
    x = x.upper()
    gc = round(100.0*(x.count("G") + x.count("C")) / len(x), 2)
    return gc

def nucleotide_stats(
    input,
    proteins,
    output
):
    stats = {}
    with gzip.open(input, "rt") as input_gunzipped:
        for r in Bio.SeqIO.parse(input_gunzipped, "fasta"):
            stats[r.id] = {
                "contig_id":r.id,
                "contig_length":len(r.seq),
                "cds_length":0,
                "gc_content":calc_gc(str(r.seq)),
                "n_count":str(r.seq).upper().count("N")
            }

    with gzip.open(proteins, "rt") as proteins_gunzipped:
        for line in proteins_gunzipped:
            if line[0] == ">":
                r = line.split()
                contig_id = r[0][1:].rsplit("_", 1)[0]
                if contig_id in stats:
                    stats[contig_id]["cds_length"] += int(r[4]) - int(r[2]) + 1
                else:
                    stats[contig_id + "_1"]["cds_length"] += int(r[4]) - int(r[2]) + 1

    with open(output, "w") as out:
        fields = ["contig_id", "contig_length", "cds_length", "cds_density", "gc_content", "n_count"]
        out.write("\t".join(fields)+"\n")
        for r in stats.values():
            r["cds_density"] = round(100*r["cds_length"]/r["contig_length"],2)
            row = [str(r[f]) for f in fields]
            out.write("\t".join(row)+"\n")


def main(args=None):
    args = parse_args(args)
    nucleotide_stats(
        args.input,
        args.proteins,
        args.output
    )

if __name__ == "__main__":
    sys.exit(main())
