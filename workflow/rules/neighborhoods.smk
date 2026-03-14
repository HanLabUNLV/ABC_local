rule create_neighborhoods:
	input:
		candidateRegions = os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted.candidateRegions.bed"),
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed')
	params:
		DHS = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "DHS"] or '',
		ATAC = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "ATAC"] or '',
		default = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'default_accessibility_feature'],
		H3K27ac = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K27ac"] or '',
		H3K4me1 = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K4me1"] or '' if "H3K4me1" in BIOSAMPLES_CONFIG.columns else '',
		H3K4me3 = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K4me3"] or '' if "H3K4me3" in BIOSAMPLES_CONFIG.columns else '',
		H3K27me3 = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K27me3"] or '' if "H3K27me3" in BIOSAMPLES_CONFIG.columns else '',
		H3K36me3 = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K36me3"] or '' if "H3K36me3" in BIOSAMPLES_CONFIG.columns else '',
		H3K9me3 = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K9me3"] or '' if "H3K9me3" in BIOSAMPLES_CONFIG.columns else '',
		genes = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'genes'],
		ubiquitous_genes = config['ref']['ubiquitous_genes'],
		chrom_sizes = config['ref']['chrom_sizes'],
		qnorm = f"--qnorm {config['ref']['qnorm']}" if config['params_neighborhoods']['use_qnorm'] else "",
		scripts_dir = SCRIPTS_DIR
	conda:
		"../envs/abcenv.yml"
	output:
		enhList = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "EnhancerList.txt"),
		geneList = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "GeneList.txt"),
		neighborhoodDirectory = directory(os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods")),
		processed_genes_file = os.path.join(RESULTS_DIR, "{biosample}", "processed_genes_file.bed"),
		supplementary_features = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "supplementary_features.tsv"),
	resources:
		tmpdir='/tmp',
		mem_mb=32*1000
	shell:
		"""
		# get sorted & unique gene list
		# intersect first to remove alternate chromosomes
		bedtools intersect -u -a {params.genes} -b {input.chrom_sizes_bed} | \
		bedtools sort -faidx {params.chrom_sizes} -i stdin | \
		uniq > {output.processed_genes_file}

		# Generate supplementary_features.tsv (only include marks with non-empty paths)
		mkdir -p $(dirname {output.supplementary_features})
		printf 'feature_name\tfile\n' > {output.supplementary_features}
		[[ -n "{params.H3K4me1}" ]] && printf 'H3K4me1\t{params.H3K4me1}\n' >> {output.supplementary_features} || true
		[[ -n "{params.H3K4me3}" ]] && printf 'H3K4me3\t{params.H3K4me3}\n' >> {output.supplementary_features} || true
		[[ -n "{params.H3K27me3}" ]] && printf 'H3K27me3\t{params.H3K27me3}\n' >> {output.supplementary_features} || true
		[[ -n "{params.H3K36me3}" ]] && printf 'H3K36me3\t{params.H3K36me3}\n' >> {output.supplementary_features} || true
		[[ -n "{params.H3K9me3}" ]] && printf 'H3K9me3\t{params.H3K9me3}\n' >> {output.supplementary_features} || true

		python {params.scripts_dir}/run.neighborhoods.py \
			--candidate_enhancer_regions {input.candidateRegions} \
			--DHS {params.DHS} \
			--ATAC {params.ATAC} \
			--default_accessibility_feature {params.default} \
			--chrom_sizes {params.chrom_sizes} \
			--chrom_sizes_bed {input.chrom_sizes_bed} \
			--outdir {output.neighborhoodDirectory} \
			--genes {output.processed_genes_file} \
			--ubiquitously_expressed_genes {params.ubiquitous_genes} \
			--H3K27ac {params.H3K27ac} \
			--supplementary_features {output.supplementary_features} \
			{params.qnorm}
		"""
