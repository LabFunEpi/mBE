#!/bin/bash
#$ -q 1-day
#$ -o /logs/$JOB_NAME.stdout
#$ -e /logs/$JOB_NAME.stderr
#$ -M wazim.ismail@gmail.com
#$ -m abe
#$ -pe threaded 16
#$ -l h_vmem=8G

. $HOME/.bash_profile

trim_galore --output_dir ${OUTDIR} ${DATADIR}/${SAMPLE}.${EXT} --fastqc

bowtie -p 8 -k 1 -m 20 --best -S -n 2 -l 65 --chunkmbs 200 ${BOWTIE1INDEX} -q ${OUTDIR}/${SAMPLE}_trimmed.fq.gz ${OUTDIR}/${SAMPLE}.sam 2> ${OUTDIR}/${SAMPLE}.bowtie.log

grep "chr" ${OUTDIR}/${SAMPLE}.sam > ${OUTDIR}/${SAMPLE}-b.sam
grep -v "chrM" ${OUTDIR}/${SAMPLE}-b.sam > ${OUTDIR}/${SAMPLE}-M.sam

samtools view -bS -q 30 ${OUTDIR}/${SAMPLE}-M.sam | samtools sort -O BAM -o ${OUTDIR}/${SAMPLE}_q30.bam

java -jar ${PICARD} MarkDuplicates I=${OUTDIR}/${SAMPLE}_q30.bam O=${OUTDIR}/${SAMPLE}_NoDup.bam M=${OUTDIR}/${SAMPLE}_picard.log REMOVE_DUPLICATES=true

# QC
java -jar ${PICARD} MarkDuplicates I=${OUTDIR}/${SAMPLE}_q30.bam O=${OUTDIR}/${SAMPLE}_DupMarked.bam M=${OUTDIR}/${SAMPLE}_picard_1.log REMOVE_DUPLICATES=false
samtools index ${OUTDIR}/${SAMPLE}_DupMarked.bam
bedtools bamtobed -i ${OUTDIR}/${SAMPLE}_DupMarked.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep --color=auto -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${OUTDIR}/${SAMPLE}_DupMarked_LC

# Cleanup
rm ${OUTDIR}/${SAMPLE}.sam ${OUTDIR}/${SAMPLE}-b.sam ${OUTDIR}/${SAMPLE}-M.sam ${OUTDIR}/${SAMPLE}_q30.bam ${OUTDIR}/${SAMPLE}_DupMarked.bam ${OUTDIR}/${SAMPLE}_picard_1.log
