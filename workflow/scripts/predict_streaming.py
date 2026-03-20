import argparse
import os
import os.path
import time

# isort:skip_file
from predictor import make_predictions  # hicstraw must be imported before pandas
import pandas as pd
from tools import determine_expressed_genes, write_params


def get_model_argument_parser():
    class formatter(
        argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
    ):
        pass

    parser = argparse.ArgumentParser(
        description="Predict enhancer relative effects.", formatter_class=formatter
    )
    readable = argparse.FileType("r")

    # Basic parameters
    parser.add_argument(
        "--enhancers",
        required=True,
        help="Candidate enhancer regions. Formatted as the EnhancerList.txt file produced by run.neighborhoods.py",
    )
    parser.add_argument(
        "--genes",
        required=True,
        help="Genes to make predictions for. Formatted as the GeneList.txt file produced by run.neighborhoods.py",
    )
    parser.add_argument(
        "--score_column",
        default="ABC.Score",
        help="Column name of score to use for thresholding",
    )
    parser.add_argument(
        "--accessibility_feature",
        default=None,
        nargs="?",
        help="If both ATAC and DHS are provided, this flag must be set to either 'DHS' or 'ATAC' signifying which datatype to use in computing activity",
    )
    parser.add_argument("--outdir", required=True, help="output directory")
    parser.add_argument("--cellType", help="Name of cell type")
    parser.add_argument("--chrom_sizes", required=True, help="Chromosome sizes file")

    # hic
    parser.add_argument(
        "--hic_file", default=None, help="HiC file: (file, web link, or directory)"
    )
    parser.add_argument("--hic_resolution", type=int, help="HiC resolution")
    parser.add_argument(
        "--hic_pseudocount_distance",
        type=int,
        required=True,
        help="A pseudocount is added equal to the powerlaw fit at this distance",
    )
    parser.add_argument(
        "--hic_type",
        default="hic",
        choices=["hic", "juicebox", "bedpe", "avg"],
        help="format of hic files",
    )
    parser.add_argument(
        "--hic_is_doubly_stochastic",
        action="store_true",
        help="If hic matrix is already doubly stochastic, can skip this step",
    )

    # Power law values
    parser.add_argument(
        "--scale_hic_using_powerlaw",
        action="store_true",
        help="Quantile normalize Hi-C values using powerlaw relationship. This parameter will rescale Hi-C contacts from the input Hi-C data (specified by --hic_gamma and --hic_scale) to match the power-law relationship of a reference cell type (specified by --hic_gamma_reference)",
    )
    parser.add_argument(
        "--hic_gamma",
        type=float,
        help="Powerlaw exponent (gamma) to scale to. Must be positive",
    )
    parser.add_argument(
        "--hic_scale",
        type=float,
        help="scale of hic data. Must be positive",
    )
    parser.add_argument(
        "--hic_gamma_reference",
        type=float,
        default=0.87,
        help="Powerlaw exponent (gamma) to scale to. Must be positive",
    )

    # Genes to run through model
    parser.add_argument(
        "--expression_cutoff",
        type=float,
        default=1,
        help="Make predictions for genes with expression higher than this value. Use of this parameter is not recommended.",
    )
    parser.add_argument(
        "--promoter_activity_quantile_cutoff",
        type=float,
        default=0.30,
        help="Quantile cutoff on promoter activity. Used to consider a gene 'expressed' in the absence of expression data",
    )

    # Output formatting
    parser.add_argument(
        "--make_all_putative",
        action="store_true",
        help="Make big file with concatenation of all genes file",
    )
    parser.add_argument(
        "--use_hdf5",
        action="store_true",
        help="Write AllPutative file in hdf5 format instead of tab-delimited",
    )

    # Parameters used in development of ABC model that should not be changed for most use cases
    parser.add_argument(
        "--window",
        type=int,
        default=5000000,
        help="Consider all candidate elements within this distance of the gene's TSS when computing ABC scores. This was a parameter optimized during development of ABC and should not typically be changed.",
    )
    parser.add_argument(
        "--tss_hic_contribution",
        type=float,
        default=100,
        help="Diagonal bin of Hi-C matrix is set to this percentage of the maximum of its neighboring bins. Default value (100%) means that the diagonal bin of the Hi-C matrix will be set exactly to the maximum of its two neighboring bins. This is a parameter used in development of the ABC model that should not typically be changed, unless experimenting with using different types of input 3D contact data that require different handling of the diagonal bin.",
    )

    # Other
    parser.add_argument(
        "--tss_slop",
        type=int,
        default=500,
        help="Distance from tss to search for self-promoters",
    )
    parser.add_argument(
        "--chromosomes",
        default="all",
        help="chromosomes to make predictions for. Defaults to intersection of all chromosomes in --genes and --enhancers",
    )
    parser.add_argument(
        "--include_chrY",
        "-y",
        action="store_true",
        help="Make predictions on Y chromosome",
    )

    return parser


def get_predict_argument_parser():
    parser = get_model_argument_parser()
    return parser


def main():
    parser = get_predict_argument_parser()
    args = parser.parse_args()

    validate_args(args)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    write_params(args, os.path.join(args.outdir, "parameters.predict.txt"))

    print("reading genes")
    genes = pd.read_csv(args.genes, sep="\t")
    genes = determine_expressed_genes(
        genes, args.expression_cutoff, args.promoter_activity_quantile_cutoff
    )
    print("reading enhancers")
    enhancers_full = pd.read_csv(args.enhancers, sep="\t")
    enhancers_column_names = ["chr", "start", "end", "name", "class", "activity_base"]
    if args.accessibility_feature not in {"ATAC", "DHS"}:
        raise ValueError("The feature has to be either ATAC or DHS!")
    normalized_activity_col = f"normalized_{args.accessibility_feature.lower()}"
    normalized_h3k27ac = "normalized_h3k27ac"

    genes_columns_to_subset = [
        "chr",
        "symbol",
        "tss",
        "Expression",
        "PromoterActivityQuantile",
        "isExpressed",
        "Ensembl_ID",
        f"{args.accessibility_feature}.RPKM.quantile.TSS1Kb",
    ]
    new_genes_column_names = [
        "chr",
        "TargetGene",
        "TargetGeneTSS",
        "TargetGeneExpression",
        "TargetGenePromoterActivityQuantile",
        "TargetGeneIsExpressed",
        "TargetGeneEnsembl_ID",
        f"{normalized_activity_col}_prom",
    ]
    if "H3K27ac.RPKM.quantile.TSS1Kb" in genes.columns:
        genes_columns_to_subset.append("H3K27ac.RPKM.quantile.TSS1Kb")
        new_genes_column_names.append(f"{normalized_h3k27ac}_prom")

    genes = genes.loc[:, genes_columns_to_subset]
    genes.columns = new_genes_column_names

    enhancers = enhancers_full.loc[:, enhancers_column_names]
    enhancers["activity_base_enh"] = enhancers_full["activity_base"]
    enhancers["activity_base_squared_enh"] = enhancers["activity_base_enh"] ** 2
    enhancers[f"{normalized_activity_col}_enh"] = enhancers_full[
        f"{normalized_activity_col}"
    ]
    if "normalized_h3K27ac" in enhancers_full.columns:
        enhancers[f"{normalized_h3k27ac}_enh"] = enhancers_full["normalized_h3K27ac"]

    # Output file paths (.tmp = uncompressed during run, gzipped at the end)
    all_pred_file_expressed = os.path.join(
        args.outdir, "EnhancerPredictionsAllPutative.tsv.gz"
    )
    all_pred_file_nonexpressed = os.path.join(
        args.outdir, "EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz"
    )
    variant_overlap_file = os.path.join(
        args.outdir,
        "EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp.tsv.gz",
    )
    tmp_expressed = all_pred_file_expressed.replace(".tsv.gz", ".tsv.tmp")
    tmp_nonexpressed = all_pred_file_nonexpressed.replace(".tsv.gz", ".tsv.tmp")
    tmp_variant = variant_overlap_file.replace(".tsv.gz", ".tsv.tmp")

    # Make predictions
    if args.chromosomes == "all":
        chromosomes = set(genes["chr"]).intersection(set(enhancers["chr"]))
        if not args.include_chrY:
            chromosomes.discard("chrY")
        chromosomes = sorted(chromosomes)
    else:
        chromosomes = args.chromosomes.split(",")

    chrom_sizes_map = pd.read_csv(
        args.chrom_sizes, sep="\t", header=None, index_col=0
    ).to_dict()[1]

    csv_kwargs = dict(
        sep="\t",
        index=False,
        float_format="%.6f",
        na_rep="NaN",
    )

    first_chrom = True
    for chromosome in chromosomes:
        print("Making predictions for chromosome: {}".format(chromosome))
        t = time.time()
        this_enh = enhancers.loc[enhancers["chr"] == chromosome, :].copy()
        this_genes = genes.loc[genes["chr"] == chromosome, :].copy()

        this_chr = make_predictions(
            chromosome,
            this_enh,
            this_genes,
            args,
            args.hic_gamma,
            args.hic_scale,
            chrom_sizes_map,
        )

        if len(this_chr) == 0:
            print("Completed chromosome: {}. Elapsed time: {} \n".format(
                chromosome, time.time() - t
            ))
            first_chrom = False
            continue

        this_chr["CellType"] = args.cellType
        if args.hic_file:
            this_chr["hic_contact_squared"] = this_chr["hic_contact"] ** 2

        # Write expressed and non-expressed rows directly, appending after first chrom.
        this_chr.loc[this_chr.TargetGeneIsExpressed, :].to_csv(
            tmp_expressed, mode="a" if not first_chrom else "w",
            header=first_chrom, **csv_kwargs,
        )
        this_chr.loc[~this_chr.TargetGeneIsExpressed, :].to_csv(
            tmp_nonexpressed, mode="a" if not first_chrom else "w",
            header=first_chrom, **csv_kwargs,
        )

        # Write variant overlap rows (same filters as test_variant_overlap).
        score_t = this_chr[args.score_column] > 0.015
        not_promoter = this_chr["class"] != "promoter"
        is_promoter = this_chr["class"] == "promoter"
        score_one = this_chr[args.score_column] > 0.1
        variant_overlap = this_chr[(score_t & not_promoter) | (is_promoter & score_one)]
        variant_overlap = variant_overlap.dropna(subset=[args.score_column])
        variant_overlap = variant_overlap.loc[variant_overlap["distance"] <= 2000000]
        variant_overlap.to_csv(
            tmp_variant, mode="a" if not first_chrom else "w",
            header=first_chrom, **csv_kwargs,
        )

        del this_chr, variant_overlap
        first_chrom = False

        print("Completed chromosome: {}. Elapsed time: {} \n".format(
            chromosome, time.time() - t
        ))

    # Gzip the tmp files now that all chromosomes are written.
    print("Compressing output files...")
    import subprocess
    for tmp, final in [
        (tmp_expressed, all_pred_file_expressed),
        (tmp_nonexpressed, all_pred_file_nonexpressed),
        (tmp_variant, variant_overlap_file),
    ]:
        if os.path.exists(tmp):
            subprocess.run(["gzip", "-f", tmp], check=True)
            os.rename(tmp + ".gz", final)

    print("Done.")


def validate_args(args):
    if args.hic_file and (args.hic_type == "juicebox" or args.hic_type == "hic"):
        assert (
            args.hic_resolution is not None
        ), "HiC resolution must be provided if hic_type is hic or juicebox"

    if not args.hic_file:
        print(
            "WARNING: Hi-C not provided. Model will only compute ABC score using powerlaw!"
        )


if __name__ == "__main__":
    main()
