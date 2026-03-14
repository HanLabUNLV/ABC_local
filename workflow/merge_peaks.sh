cat K562/Peaks/macs2_peaks.narrowPeak GM12878/Peaks/macs2_peaks.narrowPeak > K562_GM12878/Peaks/common.macs2_peaks.narrowPeak
bedtools sort -i K562_GM12878/Peaks/common.macs2_peaks.narrowPeak > K562_GM12878/Peaks/common.macs2_peaks.sorted.narrowPeak 
#bedtools merge -i Peaks/common.macs2_peaks.sorted.narrowPeak -c 2,3,5,7,8,9,10 -o collapse,collapse,mean,mean,collapse,collapse,collapse > Peaks/merged.collapsed.narrowPeak
