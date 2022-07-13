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

cr_arc_path=${softwares}/cellranger-arc-2.0.0

ref=${references}/transcriptomes/refdata-cellranger-arc-mm10-2020-A-2.0.0

cd $outdir
${cr_arc_path}/cellranger-arc aggr --id=${id} \
                                   --reference=${ref} \
                                   --csv=${libraries} \
                                   --localcores=16 \
                                   --localmem=115
