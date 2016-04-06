#!/bin/bash

PROTEINS=(E2F1 E2F4)
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
GENOMES=(hg19 hg38)
THRESH=0.207

PREDICTIONS_DIR=/Users/dcl9/Data/tf-dna-predictions/hardac-results
PYTHON_DIR=/Users/dcl9/Code/python/SVR_models
for GENOME in ${GENOMES[@]}; do
  for PROTEIN in ${PROTEINS[@]}; do
    for CHROM in ${CHROMS[@]}; do
      OUTDIR=${PREDICTIONS_DIR}/${GENOME}/central20bp/${PROTEIN} 
      mkdir -p $OUTDIR
      COMBINED=$(mktemp)
      echo     python ${PYTHON_DIR}/combine_predictions-sql.py ${PREDICTIONS_DIR}/${GENOME}-${PROTEIN}-${CHROM}-*-predictions.bed
      python ${PYTHON_DIR}/combine_predictions_sql.py ${PREDICTIONS_DIR}/${GENOME}-${PROTEIN}-${CHROM}-*-predictions.bed > $COMBINED
      FILTERED=$(mktemp)
      python ${PYTHON_DIR}/filter.py $COMBINED $THRESH > $FILTERED
      python ${PYTHON_DIR}/resize_regions.py --width 20 $FILTERED ${OUTDIR}/${GENOME}-${CHROM}-${PROTEIN}-central20bp-SVR-scores.bed
      rm $FILTERED
      rm $COMBINED
    done
  done
done
