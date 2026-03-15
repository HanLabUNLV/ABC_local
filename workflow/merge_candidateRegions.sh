#!/bin/bash
# Usage: bash merge_candidateRegions.sh <results_dir> <chrom_sizes>
# Output: <results_dir>/global/Peaks/merged.candidateRegions.bed

RESULTS_DIR=${1:-results}
CHROM_SIZES=${2:?Must provide chrom_sizes path as second argument}

OUT_DIR="${RESULTS_DIR}/global/Peaks"
mkdir -p "${OUT_DIR}"

FILES=$(ls "${RESULTS_DIR}"/*/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed 2>/dev/null)

if [[ -z "${FILES}" ]]; then
    echo "Error: No candidateRegions files found in ${RESULTS_DIR}" >&2
    exit 1
fi

cat ${FILES} \
    | bedtools sort -faidx "${CHROM_SIZES}" -i stdin \
    | bedtools merge -i stdin \
    > "${OUT_DIR}/merged.candidateRegions.bed"

echo "Merged candidateRegions written to ${OUT_DIR}/merged.candidateRegions.bed"
