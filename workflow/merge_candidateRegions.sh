cat K562/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed GM12878/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed > common.candidateRegions.bed
bedtools sort -i common.candidateRegions.bed > common.candidateRegions.sorted.bed 
bedtools merge -i common.candidateRegions.sorted.bed > common.candidateRegions.sorted.merged.bed
