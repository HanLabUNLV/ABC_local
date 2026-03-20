import argparse
import os
import subprocess

import pandas as pd

from neighborhoods import apply_tss_bins_qnorm, count_tss_bins, read_gene_bed_file, process_gene_bed


MARKS = ["H3K4me1", "H3K4me3", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3"]


def parseargs():
    class formatter(
        argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
    ):
        pass

    parser = argparse.ArgumentParser(
        description="Generate strand-aware TSS bins for histone marks",
        formatter_class=formatter,
    )
    parser.add_argument(
        "--genes", required=True,
        help="Processed genes BED file (BED6 + Ensembl_ID + gene_type columns)",
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Directory to write output files to",
    )
    parser.add_argument(
        "--chrom_sizes", required=True,
        help="Genome chromosome sizes file (2-column TSV: chr, size)",
    )
    parser.add_argument(
        "--chrom_sizes_bed", required=True,
        help="BED file of chromosome sizes (used by bedtools)",
    )
    parser.add_argument(
        "--gene_name_annotations", default="symbol",
        help="Comma-delimited name annotations in the gene BED name field",
    )
    parser.add_argument(
        "--primary_gene_identifier", default="symbol",
        help="Primary gene identifier (must be present in gene_name_annotations)",
    )
    parser.add_argument(
        "--bin_size", type=int, default=128,
        help="Bin size in bp (default: 128, matching EPInformer sequence encoder)",
    )
    parser.add_argument(
        "--slop", type=int, default=10240,
        help="Half-window size in bp around TSS (default: 10240 -> 160 bins of 128bp)",
    )
    parser.add_argument(
        "--use_secondary_counting_method", action="store_true",
        help="Use slower secondary counting method (passed to count_bam)",
    )
    parser.add_argument(
        "--qnorm_reference", default=None,
        help="Path to a TSS-bins quantile-normalization reference file "
             "(built with build_tss_bins_qnorm_ref.py). When provided, adds "
             "'{mark}.normalized_RPM' columns to GeneTSSbins.txt.",
    )
    parser.add_argument(
        "--skip_counting", action="store_true",
        help="Skip BAM counting and apply --qnorm_reference directly to an "
             "existing GeneTSSbins.txt in --outdir. Requires --qnorm_reference.",
    )

    for mark in MARKS:
        parser.add_argument(
            f"--{mark}", default="", nargs="?",
            help=f"Comma-delimited {mark} BAM file paths",
        )

    return parser.parse_args()


def get_features(args):
    features = {}
    for mark in MARKS:
        val = getattr(args, mark, "") or ""
        bams = [b.strip() for b in val.split(",") if b.strip()]
        if bams:
            features[mark] = bams
    return features


def index_bams(features):
    for bam_list in features.values():
        for bam in bam_list:
            if bam.endswith(".bam") and not os.path.exists(bam + ".bai"):
                print(f"Indexing {bam} ...", flush=True)
                subprocess.check_call(["samtools", "index", bam])


def main():
    args = parseargs()
    os.makedirs(args.outdir, exist_ok=True)

    if args.skip_counting:
        if not args.qnorm_reference:
            raise ValueError("--skip_counting requires --qnorm_reference.")
        existing = os.path.join(args.outdir, "GeneTSSbins.txt")
        if not os.path.exists(existing):
            raise FileNotFoundError(f"--skip_counting set but {existing} not found.")
        print(f"Skipping BAM counting; reading {existing}", flush=True)
        result_df = pd.read_csv(existing, sep="\t")
        result_df = apply_tss_bins_qnorm(result_df, args.qnorm_reference)
        result_df.to_csv(existing, sep="\t", index=False, float_format="%.6f")
        print(f"Updated: {existing}", flush=True)
        return

    features = get_features(args)
    if not features:
        raise ValueError("At least one histone mark BAM must be provided.")

    index_bams(features)

    # Load and process genes -- same pipeline as run.neighborhoods.py
    bed = read_gene_bed_file(args.genes)
    genes = process_gene_bed(
        bed,
        name_cols=args.gene_name_annotations,
        main_name=args.primary_gene_identifier,
        chrom_sizes=args.chrom_sizes,
        fail_on_nonunique=False,
    )

    chrom_sizes_map = pd.read_csv(
        args.chrom_sizes, sep="\t", header=None, index_col=0
    ).to_dict()[1]

    count_tss_bins(
        genes=genes,
        genome_sizes=args.chrom_sizes,
        genome_sizes_bed=args.chrom_sizes_bed,
        chrom_sizes_map=chrom_sizes_map,
        features=features,
        outdir=args.outdir,
        bin_size=args.bin_size,
        slop=args.slop,
        use_fast_count=(not args.use_secondary_counting_method),
        qnorm_reference=args.qnorm_reference,
    )

    print("TSS bin counting complete.", flush=True)


if __name__ == "__main__":
    main()
