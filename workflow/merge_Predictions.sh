#!/usr/bin/bash
set -uex
export LC_ALL=C

cat results/K562/Predictions/EnhancerPredictionsAllPutative.tsv  > results/K562/Predictions/AllPutative.txt
tail -n +2 results/K562/Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv  >> results/K562/Predictions/AllPutative.txt
cat results/GM12878/Predictions/EnhancerPredictionsAllPutative.tsv  > results/GM12878/Predictions/AllPutative.txt
tail -n +2 results/GM12878/Predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv >> results/GM12878/Predictions/AllPutative.txt

export TMPDIR='/tmp/'
export TEMP='/tmp/'
## join two cell types
tail -n +2 results/K562/Predictions/AllPutative.txt |  sort -k1,1 -k12,12n -k2,2n  > ${TEMP}/K562.sorted
tail -n +2 results/GM12878/Predictions/AllPutative.txt |  sort -k1,1 -k12,12n -k2,2n  > ${TEMP}/GM12878.sorted
paste ${TEMP}/K562.sorted ${TEMP}/GM12878.sorted > ${TEMP}/AllPutative.sorted

head -n 1 results/K562/Predictions/AllPutative.txt | sed '{s/\t/\tK562./g}' > ${TEMP}/K562.header
head -n 1 results/K562/Predictions/AllPutative.txt | sed '{s/\t/\tGM12878./g}' > ${TEMP}/GM12878.header
paste ${TEMP}/K562.header ${TEMP}/GM12878.header  > ${TEMP}/AllPutative.header

cat ${TEMP}/AllPutative.header ${TEMP}/AllPutative.sorted > results/K562_GM12878/Predictions/AllPutative.K562_GM12878.txt


