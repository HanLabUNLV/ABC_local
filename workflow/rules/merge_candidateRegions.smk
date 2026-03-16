rule merge_candidate_regions:
    input:
        candidateRegions = expand(
            os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed"),
            biosample=BIOSAMPLES_CONFIG["biosample"]
        )
    params:
        chrom_sizes = config['ref']['chrom_sizes']
    output:
        merged = os.path.join(RESULTS_DIR, "global", "Peaks", "merged.candidateRegions.bed")
    shell:
        """
        cat {input.candidateRegions} \
            | bedtools sort -faidx {params.chrom_sizes} -i stdin \
            | bedtools merge -i stdin \
            > {output.merged}
        """
