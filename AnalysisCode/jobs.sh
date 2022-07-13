##### ATAC-seq processing pipelines ########

pipeline=/mBE/AnalysisCode/ATAC_pipeline_paired.sh
PICARD=/softwares/picard/2.9.0/picard.jar
BOWTIE2INDEX=/Indexes/Bowtie2Index/hg19
GSIZE=hs
DATADIR=/data_dir
OUTDIR=/out_dir
EXT=fastq.gz

qsub -N J_2G_ATAC -v BOWTIE2INDEX=${BOWTIE2INDEX},GSIZE=${GSIZE},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=2G_ATAC ${pipeline}
qsub -N J_2M_ATAC -v BOWTIE2INDEX=${BOWTIE2INDEX},GSIZE=${GSIZE},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=2M_ATAC ${pipeline}
qsub -N HMEC_ATAC -v BOWTIE2INDEX=${BOWTIE2INDEX},GSIZE=${GSIZE},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=HMEC_ATAC ${pipeline}

ataqv --peak-file /out_dir/2G_ATAC_NP.bed --tss-file /public_data/tss-annotation/TSS/hg19.refGene.TSS.bed human /out_dir/2G_ATAC_DupMarked.bam > /out_dir/2G_ATAC_QC &
ataqv --peak-file /out_dir/2M_ATAC_NP.bed --tss-file /public_data/tss-annotation/TSS/hg19.refGene.TSS.bed human /out_dir/2M_ATAC_DupMarked.bam > /out_dir/2M_ATAC_QC &
ataqv --peak-file /out_dir/HMEC_ATAC_NP.bed --tss-file /public_data/tss-annotation/TSS/hg19.refGene.TSS.bed human /out_dir/HMEC_ATAC_DupMarked.bam > /out_dir/HMEC_ATAC_QC &

#######

BOWTIE2INDEX=/Indexes/Bowtie2Index/mm9
GSIZE=mm
qsub -N mDF_omniATAC_2 -v BOWTIE2INDEX=${BOWTIE2INDEX},GSIZE=${GSIZE},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=mDF-omniATAC-5.FCHTYGHBBXX_L8 ${pipeline}

ataqv --peak-file /out_dir/mDF-omniATAC-5.FCHTYGHBBXX_L8_NP.bed --tss-file /public_data/tss-annotation/TSS/mm9.refGene.TSS.bed human /out_dir/mDF-omniATAC-5.FCHTYGHBBXX_L8_DupMarked.bam > /out_dir/mDF-omniATAC-5.FCHTYGHBBXX_L8_QC &

#######

pipeline=/mBE/AnalysisCode/ATAC_pipeline_single.sh
PICARD=/softwares/picard/2.9.0/picard.jar
BOWTIE2INDEX=/Indexes/Bowtie2Index/mm9
GSIZE=mm
DATADIR=/data_dir
OUTDIR=/out_dir
EXT=fastq.gz

qsub -N c36 -v BOWTIE2INDEX=${BOWTIE2INDEX},GSIZE=${GSIZE},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=c36_MaSC_ATACseq ${pipeline}
qsub -N c39 -v BOWTIE2INDEX=${BOWTIE2INDEX},GSIZE=${GSIZE},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=c39_MaSC_ATACseq ${pipeline}
qsub -N d92 -v BOWTIE2INDEX=${BOWTIE2INDEX},GSIZE=${GSIZE},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=d92_MaSC_ATACseq ${pipeline}
qsub -N d93 -v BOWTIE2INDEX=${BOWTIE2INDEX},GSIZE=${GSIZE},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=d93_MaSC_ATACseq ${pipeline}

ataqv --peak-file /out_dir/c36_MaSC_ATACseq_NP.bed --tss-file /public_data/tss-annotation/TSS/mm9.refGene.TSS.bed mouse /out_dir/c36_MaSC_ATACseq_DupMarked.bam > /out_dir/c36_MaSC_ATACseq_QC &
ataqv --peak-file /out_dir/c39_MaSC_ATACseq_NP.bed --tss-file /public_data/tss-annotation/TSS/mm9.refGene.TSS.bed mouse /out_dir/c39_MaSC_ATACseq_DupMarked.bam > /out_dir/c39_MaSC_ATACseq_QC &
ataqv --peak-file /out_dir/d92_MaSC_ATACseq_NP.bed --tss-file /public_data/tss-annotation/TSS/mm9.refGene.TSS.bed mouse /out_dir/d92_MaSC_ATACseq_DupMarked.bam > /out_dir/d92_MaSC_ATACseq_QC &
ataqv --peak-file /out_dir/d93_MaSC_ATACseq_NP.bed --tss-file /public_data/tss-annotation/TSS/mm9.refGene.TSS.bed mouse /out_dir/d93_MaSC_ATACseq_DupMarked.bam > /out_dir/d93_MaSC_ATACseq_QC &

##### ChIP-seq processing pipelines ##########

pipeline=/mBE/AnalysisCode/ChIPseq_pipeline_single.sh
PICARD=/softwares/picard/2.9.0/picard.jar
BOWTIE1INDEX=/Indexes/Bowtie1_Index/hg19
DATADIR=/data_dir
OUTDIR=/out_dir
EXT=fastq.gz

qsub -N HMEC_Input -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=HMEC_Input_ATCACG_L002_R1_001 ${pipeline}
qsub -N HMEC_m1 -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=HMEC_m1_CGATGT_L003_R1_001 ${pipeline}
qsub -N HMEC_m2 -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=HMEC_m2_ATCACG_L003_R1_001 ${pipeline}

qsub -N MCF7_Input -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=Input ${pipeline}
qsub -N MCF7_m1 -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=mH2A1 ${pipeline}
qsub -N MCF7_m2 -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=mH2A2 ${pipeline}

qsub -N NHM_Input -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=NHM_Input_R1 ${pipeline}
qsub -N NHM_m1 -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=NHM_mH2A1_R1 ${pipeline}
qsub -N NHM_m2 -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=NHM_mH2A2_R1 ${pipeline}

qsub -N G_Input -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=231LG_Input ${pipeline}
qsub -N G_m1 -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=231LG_mH2A1 ${pipeline}
qsub -N G_m2 -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=231LG_mH2A2 ${pipeline}
qsub -N G_K27m3 -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=231LG_H3K27me3 ${pipeline}
qsub -N M_m2 -v BOWTIE1INDEX=${BOWTIE1INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=231LM_mH2A2 ${pipeline}

pipeline=/mBE/AnalysisCode/ChIPseq_pipeline_paired.sh
PICARD=/softwares/picard/2.9.0/picard.jar
BOWTIE2INDEX=/Indexes/Bowtie2Index/hg19
DATADIR=/data_dir
OUTDIR=/out_dir
EXT=fastq.gz

qsub -N G_Input_2 -v BOWTIE2INDEX=${BOWTIE2INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=2G_Input-10 ${pipeline}
qsub -N G_K27ac_2 -v BOWTIE2INDEX=${BOWTIE2INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=2G_K27Ac-2 ${pipeline}
qsub -N G_K4me1_2 -v BOWTIE2INDEX=${BOWTIE2INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=2G_K4me1-12 ${pipeline}
qsub -N G_p300_2 -v BOWTIE2INDEX=${BOWTIE2INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=2G_p300-16 ${pipeline}
qsub -N G_BRD4_2 -v BOWTIE2INDEX=${BOWTIE2INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=2G_Brd4-6 ${pipeline}
qsub -N M_K27ac_2 -v BOWTIE2INDEX=${BOWTIE2INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=2M_K27Ac-1 ${pipeline}
qsub -N M_K4me1_2 -v BOWTIE2INDEX=${BOWTIE2INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=2M_K4me1-11 ${pipeline}
qsub -N M_p300_2 -v BOWTIE2INDEX=${BOWTIE2INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=2M_p300-15 ${pipeline}
qsub -N M_BRD4_2 -v BOWTIE2INDEX=${BOWTIE2INDEX},PICARD=${PICARD},DATADIR=${DATADIR},OUTDIR=${OUTDIR},EXT=${EXT},SAMPLE=2M_Brd4-5 ${pipeline}

# Wait for all previous jobs before running these Macs2 calls

pipeline=/mBE/AnalysisCode/ChIPseq_pipeline_macs2.sh
CHROMSIZES=/public_data/hg19.chrom.sizes
GSIZE=hs
qsub -N HMEC_m1 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=HMEC_m1_CGATGT_L003_R1_001,CONTROL=HMEC_Input_ATCACG_L002_R1_001 ${pipeline}
qsub -N HMEC_m2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=HMEC_m2_ATCACG_L003_R1_001,CONTROL=HMEC_Input_ATCACG_L002_R1_001 ${pipeline}

qsub -N MCF7_m1 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=mH2A1,CONTROL=Input ${pipeline}
qsub -N MCF7_m2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=mH2A2,CONTROL=Input ${pipeline}

qsub -N NHM_m1 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=NHM_mH2A1_R1,CONTROL=NHM_Input_R1 ${pipeline}
qsub -N NHM_m2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=NHM_mH2A2_R1,CONTROL=NHM_Input_R1 ${pipeline}

qsub -N G_m1 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=231LG_mH2A1,CONTROL=231LG_Input ${pipeline}
qsub -N G_m2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=231LG_mH2A2,CONTROL=231LG_Input ${pipeline}
qsub -N G_K27m3 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=231LG_H3K27me3,CONTROL=231LG_Input ${pipeline}
qsub -N M_m2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=231LM_mH2A2,CONTROL=231LG_Input ${pipeline}

qsub -N G_K27ac_2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=2G_K27Ac-2,CONTROL=2G_Input-10 ${pipeline}
qsub -N G_K4me1_2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=2G_K4me1-12,CONTROL=2G_Input-10 ${pipeline}
qsub -N G_p300_2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=2G_p300-16,CONTROL=2G_Input-10 ${pipeline}
qsub -N G_BRD4_2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=2G_Brd4-6,CONTROL=2G_Input-10 ${pipeline}
qsub -N M_K27ac_2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=2M_K27Ac-1,CONTROL=2G_Input-10 ${pipeline}
qsub -N M_K4me1_2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=2M_K4me1-11,CONTROL=2G_Input-10 ${pipeline}
qsub -N M_p300_2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=2M_p300-15,CONTROL=2G_Input-10 ${pipeline}
qsub -N M_BRD4_2 -v OUTDIR=${OUTDIR},CHROMSIZES=${CHROMSIZES},GSIZE=${GSIZE},TARGET=2M_Brd4-5,CONTROL=2G_Input-10 ${pipeline}
