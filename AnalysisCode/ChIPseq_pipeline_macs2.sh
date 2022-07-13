#!/bin/bash
#$ -q 1-day
#$ -o /logs/$JOB_NAME.stdout
#$ -e /logs/$JOB_NAME.stderr
#$ -M wazim.ismail@gmail.com
#$ -m abe
#$ -pe threaded 16
#$ -l h_vmem=8G

. $HOME/.bash_profile

macs2 callpeak -t ${OUTDIR}/${TARGET}_NoDup.bam -c ${OUTDIR}/${CONTROL}_NoDup.bam -f BAM -g ${GSIZE} -n ${CONTROL} --outdir ${OUTDIR} --keep-dup all --nomodel --bdg --SPMR 2> ${OUTDIR}/${CONTROL}.macs2.log
awk '{print $1,$2,$3,$4}' OFS='\t' ${OUTDIR}/${CONTROL}_peaks.narrowPeak > ${OUTDIR}/${CONTROL}_NP.bed

macs2 bdgcmp -t ${OUTDIR}/${CONTROL}_treat_pileup.bdg -c ${OUTDIR}/${CONTROL}_control_lambda.bdg -m FE --o-prefix ${OUTDIR}/${CONTROL}
bedClip ${OUTDIR}/${CONTROL}_FE.bdg ${CHROMSIZES} ${OUTDIR}/${CONTROL}_FE_clipped.bdg
sort -k1,1 -k2,2n ${OUTDIR}/${CONTROL}_FE_clipped.bdg > ${OUTDIR}/${CONTROL}_FE_clipped_sorted.bdg
bedGraphToBigWig ${OUTDIR}/${CONTROL}_FE_clipped_sorted.bdg ${CHROMSIZES} ${OUTDIR}/${CONTROL}_FE.bw
