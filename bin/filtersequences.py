#!/usr/bin/env python

#--------------------------------------
# Imports
#--------------------------------------
import argparse
import os
import gzip
import sys

import pandas as pd

from Bio import SeqIO
from pathlib import Path


def parse_args(args=None):
    Description = "Filter viral sequences based on aggregated quality information."
    Epilog = """
    Example usage:
    python filtersequences.py \\
        --input input.fasta \\
        --genomad_summary genomad_summary.tsv \\
        --quality_summary quality_summary.tsv \\
        --tantan tantan.tsv \\
        --nucleotide_stats nuc_stats.tsv \\
        --output_fasta test_filtered.fasta \\
        --output_tsv test_filtering_data.tsv \\
        ${args}
    """

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-i",
        "--input",
        help="Path to input FASTA (gzipped) containing terminal repeat sequences.",
    )
    parser.add_argument(
        "-g",
        "--genomad_summary",
        help="Path to input TSV file containing geNomad summary file.",
    )
    parser.add_argument(
        "-q",
        "--quality_summary",
        help="Path to input TSV file containing CheckV's quality_summary.",
    )
    parser.add_argument(
        "-t",
        "--tantan",
        help="Path to input TSV file containing tantan's output.",
    )
    parser.add_argument(
        "-nt",
        "--nucleotide_stats",
        help="Path to input TSV file containing nucleotide stats (GC %, CDS density, and 'N' percentage).",
    )
    parser.add_argument(
        "-k",
        "--contigs_to_keep",
        help="Path to input TSV file containing contigs to keep, regardless of whether they pass filters.",
    )
    parser.add_argument(
        "--min_length",
        type=int,
        default=10000,
        help="Minimum length to keep a sequence regardless of completeness.",
    )
    parser.add_argument(
        "--min_completeness",
        type=float,
        default=50,
        help="Minimum CheckV completeness to keep a sequence.",
    )
    parser.add_argument(
        "--keep_proviruses",
        action="store_true",
        help="Keep proviruses.",
    )
    parser.add_argument(
        "--ignore_warnings",
        action="store_true",
        help="Ignore warnings.",
    )
    parser.add_argument(
        "--min_cds_density",
        type=float,
        default=70,
        help="Minimum coding density keep a sequence.",
    )
    parser.add_argument(
        "--max_tantan_frac",
        type=float,
        default=5,
        help="Maximum tantan percentage to keep a sequence.",
    )
    parser.add_argument(
        "--max_kmer_freq",
        type=float,
        default=1.15,
        help="Maximum kmer frequency to keep a sequence.",
    )
    parser.add_argument(
        "--max_n_frac",
        type=float,
        default=5,
        help="Maximum 'N' percentage to keep a sequence.",
    )
    parser.add_argument(
        "-s",
        "--sample",
        help="Sample prefix to add to contig_id.",
    )
    parser.add_argument(
        "-of",
        "--output_fasta",
        help="Output filtered FASTA file.",
    )
    parser.add_argument(
        "-ot",
        "--output_tsv",
        help="Output TSV file containing data used for filtering.",
    )
    return parser.parse_args(args)


#--------------------------------------
# Load and merge virus data
#--------------------------------------
def combine_virus_data(
    genomad_summary: Path,
    quality_summary: Path,
    tantan: Path,
    nucleotide_stats: Path,
) -> Path:
    """
    Combine virus data from multiple sources into a single TSV file.

    Args:
        genomad_summary (Path)      : Path to input TSV file containing geNomad summary file.
        quality_summary (Path)      : Path to input TSV file containing CheckV's quality_summary.
        tantan (Path)               : Path to input TSV file containing tantan's output.
        nucleotide_stats (Path)     : Path to input TSV file containing nucleotide stats (GC %, CDS density, and 'N' percentage).

    Returns:
        Path: Path to output TSV file containing all pertinent virus data combined across sources.
    """
    # load genomad virus summary data
    if os.stat(genomad_summary).st_size > 0:
        genomad_summary_df = pd.read_csv(
            genomad_summary,
            sep="\t",
            header=0,
            index_col="seq_name",
            usecols=["seq_name", "topology", "virus_score", "fdr", "taxonomy"]
        )
    else:
        genomad_summary_df = pd.DataFrame(columns=["topology", "virus_score", "fdr", "taxonomy"])

    # load checkv quality summary data
    if os.stat(quality_summary).st_size > 0:
        quality_summary_df = pd.read_csv(
            quality_summary,
            sep="\t",
            header=0,
            index_col="contig_id",
            usecols=["contig_id", "contig_length", "provirus", "completeness", "completeness_method", "kmer_freq", "warnings"]
        )
    else:
        quality_summary_df = pd.DataFrame(columns=["contig_length", "provirus", "completeness", "completeness_method", "kmer_freq", "warnings"])

    # load tantan data
    if os.stat(tantan).st_size > 0:
        tantan_dict = dict([[contig_id, 0] for contig_id in genomad_summary_df.index])
        for line in open(tantan):
            contig_id, start, end = line.split()
            length = int(end) - int(start) + 1
            if contig_id in tantan_dict:
                tantan_dict[contig_id] += length
            else:
                tantan_dict[contig_id.rpartition('_')[0]] = length
        tantan_df = pd.DataFrame.from_dict(tantan_dict, orient="index")
        tantan_df.columns = ["tantan_len"]
    else:
        tantan_df = pd.DataFrame(columns=["tantan_len"])

    # load nucleotide stats
    nucleotide_stats_df = pd.read_csv(
        nucleotide_stats,
        sep="\t",
        header=0,
        index_col="contig_id",
        usecols=["contig_id", "cds_length", "cds_density", "gc_content", "n_count"]
    )
    # combine quality data from all sources by contig_id (index)
    comb_virus_data_df = pd.concat([
        genomad_summary_df,
        quality_summary_df,
        tantan_df,
        nucleotide_stats_df
    ], axis=1)

    comb_virus_data_df.insert(11, "tantan_perc", (comb_virus_data_df['tantan_len'] / comb_virus_data_df['contig_length']))
    comb_virus_data_df.insert(len(comb_virus_data_df.columns), "n_perc", (comb_virus_data_df['n_count'] / comb_virus_data_df['contig_length']))

    return comb_virus_data_df


#--------------------------------------
# Filter sequences based on classification, completeness, and sequence quality
#--------------------------------------
def filter_sequences(
    input_fasta: Path,
    contigs_to_keep: Path,
    comb_virus_data_df: pd.DataFrame,
    min_length: int,
    min_completeness: float,
    keep_proviruses: bool,
    ignore_warnings: bool,
    min_cds_density: float,
    max_tantan_freq: float,
    max_kmer_freq: float,
    max_n_freq: float,
    sample: str,
    output_fasta: Path,
    output_tsv: Path
) -> None:
    """
    Filter viral sequences based on classification and write to output FASTA/TSV files.

    Args:
        input_fasta (Path)              : Path to input FASTA file containing terminal repeat sequences.
        contigs_to_keep (Path)          : Path to input TSV file (no header) containing contigs to keep (1 contig id per row), regardless of whether they pass filters.
        comb_virus_data_df (DataFrame)  : DataFrame containing all pertinent virus data combined across sources.
        min_genomad_score (float)      : Minimum geNomad viral score to classify a sequence as viral.
        max_genomad_fdr (float)         : Maximum geNomad FDR to classify a sequence as viral.
        min_length (int)                : Minimum length to keep a sequence.
        keep_proviruses (bool)          : Whether to keep proviruses.
        ignore_warnings (bool)          : Whether to ignore warnings.
        min_completeness (float)        : Minimum CheckV completeness to keep a sequence.
        min_hmm_completeness (float)    : Minimum HMM completeness to keep a sequence.
        min_cds_density (float)         : Minimum coding density to keep a sequence.
        max_tantan_freq (float)         : Maximum tantan frequency to keep a sequence.
        max_kmer_freq (float)           : Maximum kmer frequency to keep a sequence.
        max_n_freq (float)              : Maximum 'N' frequency to keep a sequence.
        sample (str)                    : Sample prefix to add to contig_id.
        output_fasta (Path)             : Path to output FASTA file containing filtered viral sequences.
        output_tsv (Path)               : Path to output TSV file containing filtered viral data.

    Returns:
        Outputs a TSV file with combined virus data, and a FASTA file containing the filtered viral sequences.
    """
    # load contigs to keep
    if os.path.exists(str(contigs_to_keep)):
        contigs_to_keep_df = pd.read_csv(
            contigs_to_keep,
            sep="\t",
            header=None,
            names=["contig_name"],
        )
        contigs_to_keep_set = set(contigs_to_keep_df["contig_name"])
    else:
        contigs_to_keep_set = set()

    if not keep_proviruses:
        comb_virus_data_df = comb_virus_data_df.query("`provirus` == 'No' and `topology` != 'Provirus'")

    if not ignore_warnings:
        comb_virus_data_df = comb_virus_data_df.query("`warnings` != 'contig >1.5x longer than expected genome length'")

    # filter viral sequences by classification, completeness, and quality
    is_hq_filter = "(`completeness` >= @min_completeness or `contig_length` >= @min_length)"
    seq_filter = "(`cds_density` >= @min_cds_density and `tantan_perc` <= @max_tantan_freq and `kmer_freq` <= @max_kmer_freq and `n_perc` <= @max_n_freq)"
    contigs_passing_filters_set = set(
        comb_virus_data_df.query(
            f"{is_hq_filter} and {seq_filter}"
        ).index
    )

    filtered_contigs_set = contigs_to_keep_set.union(contigs_passing_filters_set)
    # write filtered viral sequences to output FASTA
    with gzip.open(input_fasta, "rt") as input_gunzipped, open(output_fasta, "w") as output:
        for r in SeqIO.parse(input_gunzipped, "fasta"):
            if r.id in filtered_contigs_set:
                r.id = "sample_" + sample + "|" + r.id
                SeqIO.write(r, output, "fasta")
    # write combined metadata to output TSV
    comb_virus_data_df.insert(0, "contig_id", comb_virus_data_df.index)
    comb_virus_data_df.to_csv(output_tsv, sep="\t", index=False)
    return


def main(args=None):
    args = parse_args(args)
    combined_virus_data_df = combine_virus_data(
        args.genomad_summary,
        args.quality_summary,
        args.tantan,
        args.nucleotide_stats,
    )
    filter_sequences(
        args.input,
        args.contigs_to_keep,
        combined_virus_data_df,
        args.min_length,
        args.min_completeness,
        args.keep_proviruses,
        args.ignore_warnings,
        args.min_cds_density,
        args.max_tantan_frac,
        args.max_kmer_freq,
        args.max_n_frac,
        args.sample,
        args.output_fasta,
        args.output_tsv
    )


if __name__ == "__main__":
    sys.exit(main())
