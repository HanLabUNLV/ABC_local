MERGED_PEAKS_DIR = os.path.join(RESULTS_DIR, "global", "Peaks")

rule cat_and_sort_peaks:
	input:
		narrowPeaks = expand(
			os.path.join(RESULTS_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted"),
			biosample=BIOSAMPLES_CONFIG["biosample"]
		)
	params:
		chrom_sizes = config['ref']['chrom_sizes']
	output:
		catSorted = os.path.join(MERGED_PEAKS_DIR, "all_peaks.narrowPeak.sorted")
	shell:
		"cat {input.narrowPeaks} | bedtools sort -faidx {params.chrom_sizes} -i stdin > {output.catSorted}"

rule bedtools_merge_peaks:
	input:
		catSorted = os.path.join(MERGED_PEAKS_DIR, "all_peaks.narrowPeak.sorted")
	output:
		collapsed = os.path.join(MERGED_PEAKS_DIR, "merged.collapsed.narrowPeak")
	shell:
		"""
		bedtools merge -i {input.catSorted} \
			-c 2,3,5,7,8,9,10 \
			-o collapse,collapse,mean,mean,collapse,collapse,collapse \
			> {output.collapsed}
		"""

rule process_merged_peaks:
	input:
		collapsed = os.path.join(MERGED_PEAKS_DIR, "merged.collapsed.narrowPeak")
	params:
		mergepeaks_script = os.path.join(ABC_DIR_PATH, "workflow", "mergepeaks.py")
	output:
		merged = os.path.join(MERGED_PEAKS_DIR, "merged.narrowPeak")
	shell:
		"python {params.mergepeaks_script} --input {input.collapsed} --output {output.merged}"
