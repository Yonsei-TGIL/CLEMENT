#!/bin/bash

#PATH#
INPUT_PATH=$1
OUTPUT_PATH=$2

PON=1000g_pon.hg38.M${num}.rem.vcf.gz

REF=/reference/human/NCBI/GRCh38_GATK/BWAIndex/genome.fa
INTERVAL=/TargetRegion/SureSelect_v5.hg38.new.bed
library_name=$3

##CODE
#1. Mutect2
gatk Mutect2 \
-R ${REF} \
-I ${INPUT_PATH}'/'${library_name}'_RF_il_WES.RGadded.marked.realigned.fixed.recal.LA.woClp.bam' \
-L ${INTERVAL} \
--panel-of-normals ${PON} \
-O ${OUTPUT_PATH}/${library_name}'.MT.whole.vcf.gz' \
--tmp-dir ${TMP_PATH}

#2. FilterMutectCalls
gatk FilterMutectCalls \
-R ${REF} \
-V ${OUTPUT_PATH}/${library_name}'.MT.whole.vcf.gz' \
-O ${FMC_PATH}/${library_name}'.MT.whole.filtered.vcf.gz'

cat ${FMC_PATH}/${library_name}'.MT.whole.filtered.vcf.gz' | zgrep -a -E "PASS|#" > ${FMC_PATH}/${library_name}'.MT.whole.filtered.PASS.vcf'

gzip ${FMC_PATH}/${library_name}'.MT.whole.filtered.PASS.vcf'
