#!/bin/bash
#$ -q 1-day
#$ -o /logs/$JOB_NAME.stdout
#$ -e /logs/$JOB_NAME.stderr
#$ -M wazim.ismail@gmail.com
#$ -m abe
#$ -pe threaded 16
#$ -l h_vmem=8G

. $HOME/.bash_profile

trim_galore ${DATADIR}/${SAMPLE}.${EXT} --fastqc --output_dir ${OUTDIR}

bowtie2 -p 8 -x ${BOWTIE2INDEX} -U ${OUTDIR}/${SAMPLE}_trimmed.fq.gz -S ${OUTDIR}/${SAMPLE}.sam -q 2> ${OUTDIR}/${SAMPLE}_bowtie2.log

grep "chr" ${OUTDIR}/${SAMPLE}.sam > ${OUTDIR}/${SAMPLE}-b.sam
grep -v "chrM" ${OUTDIR}/${SAMPLE}-b.sam > ${OUTDIR}/${SAMPLE}-M.sam

samtools view -bS -q 30 ${OUTDIR}/${SAMPLE}-M.sam | samtools sort -O BAM -o ${OUTDIR}/${SAMPLE}_q30.bam

java -jar ${PICARD} MarkDuplicates I=${OUTDIR}/${SAMPLE}_q30.bam O=${OUTDIR}/${SAMPLE}_NoDup.bam M=${OUTDIR}/${SAMPLE}_picard.log REMOVE_DUPLICATES=true

macs2 callpeak --nomodel -t ${OUTDIR}/${SAMPLE}_NoDup.bam -g ${GSIZE} -n ${SAMPLE} --outdir ${OUTDIR} -q 5e-5 --nolambda --keep-dup all --slocal 10000 2> ${OUTDIR}/${SAMPLE}_macs2.log
awk '{print $1,$2,$3,$4}' OFS='\t' ${OUTDIR}/${SAMPLE}_peaks.narrowPeak > ${OUTDIR}/${SAMPLE}_NP.bed

samtools index ${OUTDIR}/${SAMPLE}_NoDup.bam
bamCoverage --outFileFormat=bigwig --binSize 10 --scaleFactor 0.5 --skipNonCoveredRegions --normalizeUsing RPKM -b ${OUTDIR}/${SAMPLE}_NoDup.bam -o ${OUTDIR}/${SAMPLE}.bw

# QC
java -jar ${PICARD} MarkDuplicates I=${OUTDIR}/${SAMPLE}_q30.bam O=${OUTDIR}/${SAMPLE}_DupMarked.bam M=${OUTDIR}/${SAMPLE}_picard_1.log REMOVE_DUPLICATES=false
samtools index ${OUTDIR}/${SAMPLE}_DupMarked.bam
bedtools bamtobed -i ${OUTDIR}/${SAMPLE}_DupMarked.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep --color=auto -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${OUTDIR}/${SAMPLE}_DupMarked_LC

# Cleanup
rm ${OUTDIR}/${SAMPLE}.sam ${OUTDIR}/${SAMPLE}-b.sam ${OUTDIR}/${SAMPLE}-M.sam ${OUTDIR}/${SAMPLE}_q30.bam ${OUTDIR}/${SAMPLE}_NoDup.bam ${OUTDIR}/${SAMPLE}_picard_1.log
