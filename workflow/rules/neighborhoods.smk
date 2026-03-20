rule create_neighborhoods:
	input:
		candidateRegions = os.path.join(RESULTS_DIR, "global", "Peaks", "merged.candidateRegions.bed"),
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
		supplementary_features = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, 'supplementary_features'] if 'supplementary_features' in BIOSAMPLES_CONFIG.columns else '',
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
	resources:
		tmpdir='/tmp',
		mem_mb=32*1000
	shell:
		"""
		# get sorted & unique gene list
		# intersect first to remove alternate chromosomes
		for bam in $(echo "{params.DHS},{params.ATAC},{params.H3K27ac},{params.H3K4me1},{params.H3K4me3},{params.H3K27me3},{params.H3K36me3},{params.H3K9me3}" | tr ',' ' '); do
			[[ -z "${{bam}}" ]] || [[ -f "${{bam}}.bai" ]] || samtools index "${{bam}}"
		done

		bedtools intersect -u -a {params.genes} -b {input.chrom_sizes_bed} | \
		bedtools sort -faidx {params.chrom_sizes} -i stdin | \
		uniq > {output.processed_genes_file}

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
			$([ -n "{params.supplementary_features}" ] && echo "--supplementary_features {params.supplementary_features}") \
			{params.qnorm}
		"""


rule create_tss_bins:
	input:
		processed_genes = os.path.join(RESULTS_DIR, "{biosample}", "processed_genes_file.bed"),
		chrom_sizes_bed = os.path.join(RESULTS_DIR, "tmp", os.path.basename(config['ref']['chrom_sizes']) + '.bed'),
	params:
		H3K4me1  = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K4me1"]  or '' if "H3K4me1"  in BIOSAMPLES_CONFIG.columns else '',
		H3K4me3  = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K4me3"]  or '' if "H3K4me3"  in BIOSAMPLES_CONFIG.columns else '',
		H3K9me3  = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K9me3"]  or '' if "H3K9me3"  in BIOSAMPLES_CONFIG.columns else '',
		H3K27ac  = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K27ac"]  or '' if "H3K27ac"  in BIOSAMPLES_CONFIG.columns else '',
		H3K27me3 = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K27me3"] or '' if "H3K27me3" in BIOSAMPLES_CONFIG.columns else '',
		H3K36me3 = lambda wildcards: BIOSAMPLES_CONFIG.loc[wildcards.biosample, "H3K36me3"] or '' if "H3K36me3" in BIOSAMPLES_CONFIG.columns else '',
		chrom_sizes    = config['ref']['chrom_sizes'],
		outdir         = lambda wildcards: os.path.join(RESULTS_DIR, wildcards.biosample, "Neighborhoods"),
		scripts_dir    = SCRIPTS_DIR,
		qnorm_reference = config.get('tss_bins_qnorm_reference', ''),
	conda:
		"../envs/abcenv.yml"
	output:
		tss_bins = os.path.join(RESULTS_DIR, "{biosample}", "Neighborhoods", "GeneTSSbins.txt"),
	resources:
		tmpdir = '/tmp',
		mem_mb = 32*1000,
	shell:
		"""
		for bam in $(echo "{params.H3K4me1},{params.H3K4me3},{params.H3K9me3},{params.H3K27ac},{params.H3K27me3},{params.H3K36me3}" | tr ',' ' '); do
			[[ -z "${{bam}}" ]] || [[ -f "${{bam}}.bai" ]] || samtools index "${{bam}}"
		done

		python {params.scripts_dir}/run.tss_bins.py \
			--genes {input.processed_genes} \
			--outdir {params.outdir} \
			--chrom_sizes {params.chrom_sizes} \
			--chrom_sizes_bed {input.chrom_sizes_bed} \
			$([ -n "{params.H3K4me1}"  ] && echo "--H3K4me1 {params.H3K4me1}")  \
			$([ -n "{params.H3K4me3}"  ] && echo "--H3K4me3 {params.H3K4me3}")  \
			$([ -n "{params.H3K9me3}"  ] && echo "--H3K9me3 {params.H3K9me3}")  \
			$([ -n "{params.H3K27ac}"  ] && echo "--H3K27ac {params.H3K27ac}")  \
			$([ -n "{params.H3K27me3}" ] && echo "--H3K27me3 {params.H3K27me3}") \
			$([ -n "{params.H3K36me3}" ] && echo "--H3K36me3 {params.H3K36me3}") \
			$([ -n "{params.qnorm_reference}" ] && echo "--qnorm_reference {params.qnorm_reference}")
		"""
