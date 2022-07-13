#!/bin/bash
#$ -q 1-day
#$ -o /logs/$JOB_NAME.stdout
#$ -e /logs/$JOB_NAME.stderr
#$ -M wazim.ismail@gmail.com
#$ -m abe
#$ -pe threaded 16
#$ -l h_vmem=8G

. $HOME/.bash_profile

myhome=/
softwares=${myhome}/softwares
references=${myhome}/references

cr_path=${softwares}/cellranger-atac-2.0.0/bin

ref=${references}/transcriptomes/refdata-cellranger-arc-GRCh38-2020-A-2.0.0

cd ${outdir}
${cr_path}/cellranger-atac count --id=${id} \
                                 --reference=${ref} \
                                 --fastqs=${readsdir} \
                                 --sample=${sample} \
                                 --localcores=16 \
                                 --localmem=115

