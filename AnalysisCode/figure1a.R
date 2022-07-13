# https://egg2.wustl.edu/roadmap/web_portal/imputed.html#imp_sig
# For imputing epigenomic data sets, we used a new method, ChromImpute (Ernst and Kellis, Nature Biotech 2015), that to predict a target mark in a target reference epigenome combines information about other marks mapped in the target reference epigenome, and the target mark at the same position in similar reference epigenomes through an ensemble of regression trees to make predictions about unobserved datasets. Additionally an imputed version of each observed data set was generated without using the corresponding observed data. Each directory corresponds to a Mark and contains a set of BigWig files one for each reference epigenome. The signal tracks for the Histone Modifications and DNase are based on the p-value signal tracks.

# https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html#scale-regions
# This tool calculates scores per genome regions. In the scale-regions mode, all regions in the BED file are stretched or shrunken to the length (in bases) indicated by the user. We use 50 bp. --binSize 50, Length, in bases, of the non-overlapping bins for averaging the score over the regions length; --regionBodyLength 50, Distance in bases to which all regions will be fit. We get one "average" score per BED region. 

# Median signal scores for each chromatin state across all genomic regions in the state, were plotted as heatmap using the R program heatmap.2 with parameters (scale='column'). Mann-Whitney U test with Bonferroni correction was performed to calculate the statistical significance of the difference in scores in each state compared to those in all other states, for each histone mark or variant. Marks that were not significantly different in a state (i.e., p > 0.05) are marked with a circle in the heatmap. 

######### Downloading Roadmap Data for HMEC, NHM and HepG2 ############## 

# wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/E059_25_imputed12marks_stateno.bed.gz -O /dev/stdout | gunzip | awk '{print $1, $2, $3, $4, "0", "+"}' | awk '{$1=$1}1' OFS="\t" | head -n -1 > E059_25_imputed12marks_stateno.bed
# wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/E118_25_imputed12marks_stateno.bed.gz -O /dev/stdout | gunzip | awk '{print $1, $2, $3, $4, "0", "+"}' | awk '{$1=$1}1' OFS="\t" | head -n -1 > E118_25_imputed12marks_stateno.bed
# wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/E119_25_imputed12marks_stateno.bed.gz -O /dev/stdout | gunzip | awk '{print $1, $2, $3, $4, "0", "+"}' | awk '{$1=$1}1' OFS="\t" | head -n -1 > E119_25_imputed12marks_stateno.bed

# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H2BK12ac/E059-H2BK12ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H2BK5ac/E059-H2BK5ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H2BK120ac/E059-H2BK120ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E059-H3K27ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27me3/E059-H3K27me3.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K36me3/E059-H3K36me3.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K79me2/E059-H3K79me2.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H2A.Z/E059-H2A.Z.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K4me1/E059-H3K4me1.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K4me3/E059-H3K4me3.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K9me3/E059-H3K9me3.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K9ac/E059-H3K9ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H4K8ac/E059-H4K8ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H4K20me1/E059-H4K20me1.imputed.pval.signal.bigwig

# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H2BK12ac/E118-H2BK12ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H2BK5ac/E118-H2BK5ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H2BK120ac/E118-H2BK120ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E118-H3K27ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27me3/E118-H3K27me3.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K36me3/E118-H3K36me3.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K79me2/E118-H3K79me2.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H2A.Z/E118-H2A.Z.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K4me1/E118-H3K4me1.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K4me3/E118-H3K4me3.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K9me3/E118-H3K9me3.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K9ac/E118-H3K9ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H4K8ac/E118-H4K8ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H4K20me1/E118-H4K20me1.imputed.pval.signal.bigwig

# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H2BK12ac/E119-H2BK12ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H2BK5ac/E119-H2BK5ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H2BK120ac/E119-H2BK120ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27ac/E119-H3K27ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K27me3/E119-H3K27me3.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K36me3/E119-H3K36me3.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K79me2/E119-H3K79me2.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H2A.Z/E119-H2A.Z.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K4me1/E119-H3K4me1.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K4me3/E119-H3K4me3.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K9me3/E119-H3K9me3.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H3K9ac/E119-H3K9ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H4K8ac/E119-H4K8ac.imputed.pval.signal.bigwig
# wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/H4K20me1/E119-H4K20me1.imputed.pval.signal.bigwig

# HMEC_mH2A1_FE.bigwig, HMEC_mH2A2_FE.bigwig, NHM_mH2A1_FE.bigwig, NHM_mH2A2_FE.bigwig obtained by running ChIPseq pipeline on data generated in our lab. 

# HepG2 mH2A1 and mH2A2 obtained from: 
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58175
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1402nnn/GSM1402782/suppl/GSM1402782_M1.HepG2.MACS.Bt2.bw -O HepG2_mH2A1.bigwig
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1402nnn/GSM1402783/suppl/GSM1402783_M2.HepG2.MACS.Bt2.bw -O HepG2_mH2A2.bigwig

################## All computeMatrix commands #####################

# # NHM
# computeMatrix scale-regions -S E059-H2BK12ac.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H2BK12ac \
                            # --outFileName NHM_H2BK12ac.tab.gz
# computeMatrix scale-regions -S E059-H2BK5ac.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H2BK5ac \
                            # --outFileName NHM_H2BK5ac.tab.gz
# computeMatrix scale-regions -S E059-H2BK120ac.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H2BK120ac \
                            # --outFileName NHM_H2BK120ac.tab.gz
# computeMatrix scale-regions -S E059-H3K27ac.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H3K27ac \
                            # --outFileName NHM_H3K27ac.tab.gz
# computeMatrix scale-regions -S E059-H3K27me3.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H3K27me3 \
                            # --outFileName NHM_H3K27me3.tab.gz
# computeMatrix scale-regions -S E059-H3K36me3.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H3K36me3 \
                            # --outFileName NHM_H3K36me3.tab.gz
# computeMatrix scale-regions -S E059-H3K79me2.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H3K79me2 \
                            # --outFileName NHM_H3K79me2.tab.gz
# computeMatrix scale-regions -S E059-H2A.Z.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H2AZ \
                            # --outFileName NHM_H2AZ.tab.gz
# computeMatrix scale-regions -S E059-H3K4me1.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H3K4me1 \
                            # --outFileName NHM_H3K4me1.tab.gz
# computeMatrix scale-regions -S E059-H3K4me3.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H3K4me3 \
                            # --outFileName NHM_H3K4me3.tab.gz
# computeMatrix scale-regions -S E059-H3K9me3.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H3K9me3 \
                            # --outFileName NHM_H3K9me3.tab.gz
# computeMatrix scale-regions -S E059-H3K9ac.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H3K9ac \
                            # --outFileName NHM_H3K9ac.tab.gz
# computeMatrix scale-regions -S E059-H4K8ac.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H4K8ac \
                            # --outFileName NHM_H4K8ac.tab.gz
# computeMatrix scale-regions -S E059-H4K20me1.imputed.pval.signal.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_H4K20me1 \
                            # --outFileName NHM_H4K20me1.tab.gz
# computeMatrix scale-regions -S NHM_mH2A1_FE.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_mH2A1 \
                            # --outFileName NHM_mH2A1.tab.gz
# computeMatrix scale-regions -S NHM_mH2A2_FE.bigwig \
                            # -R E059_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix NHM_mH2A2 \
                            # --outFileName NHM_mH2A2.tab.gz

# # HepG2
# computeMatrix scale-regions -S E118-H2BK12ac.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H2BK12ac \
                            # --outFileName HepG2_H2BK12ac.tab.gz
# computeMatrix scale-regions -S E118-H2BK5ac.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H2BK5ac \
                            # --outFileName HepG2_H2BK5ac.tab.gz
# computeMatrix scale-regions -S E118-H2BK120ac.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H2BK120ac \
                            # --outFileName HepG2_H2BK120ac.tab.gz
# computeMatrix scale-regions -S E118-H3K27ac.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H3K27ac \
                            # --outFileName HepG2_H3K27ac.tab.gz
# computeMatrix scale-regions -S E118-H3K27me3.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H3K27me3 \
                            # --outFileName HepG2_H3K27me3.tab.gz
# computeMatrix scale-regions -S E118-H3K36me3.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H3K36me3 \
                            # --outFileName HepG2_H3K36me3.tab.gz
# computeMatrix scale-regions -S E118-H3K79me2.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H3K79me2 \
                            # --outFileName HepG2_H3K79me2.tab.gz
# computeMatrix scale-regions -S E118-H2A.Z.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H2AZ \
                            # --outFileName HepG2_H2AZ.tab.gz
# computeMatrix scale-regions -S E118-H3K4me1.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H3K4me1 \
                            # --outFileName HepG2_H3K4me1.tab.gz
# computeMatrix scale-regions -S E118-H3K4me3.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H3K4me3 \
                            # --outFileName HepG2_H3K4me3.tab.gz
# computeMatrix scale-regions -S E118-H3K9me3.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H3K9me3 \
                            # --outFileName HepG2_H3K9me3.tab.gz
# computeMatrix scale-regions -S E118-H3K9ac.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H3K9ac \
                            # --outFileName HepG2_H3K9ac.tab.gz
# computeMatrix scale-regions -S E118-H4K8ac.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H4K8ac \
                            # --outFileName HepG2_H4K8ac.tab.gz
# computeMatrix scale-regions -S E118-H4K20me1.imputed.pval.signal.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_H4K20me1 \
                            # --outFileName HepG2_H4K20me1.tab.gz
# computeMatrix scale-regions -S HepG2_mH2A1.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_mH2A1 \
                            # --outFileName HepG2_mH2A1.tab.gz
# computeMatrix scale-regions -S HepG2_mH2A2.bigwig \
                            # -R E118_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HepG2_mH2A2 \
                            # --outFileName HepG2_mH2A2.tab.gz

# # HMEC
# computeMatrix scale-regions -S E119-H2BK12ac.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H2BK12ac \
                            # --outFileName HMEC_H2BK12ac.tab.gz
# computeMatrix scale-regions -S E119-H2BK5ac.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H2BK5ac \
                            # --outFileName HMEC_H2BK5ac.tab.gz
# computeMatrix scale-regions -S E119-H2BK120ac.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H2BK120ac \
                            # --outFileName HMEC_H2BK120ac.tab.gz
# computeMatrix scale-regions -S E119-H3K27ac.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H3K27ac \
                            # --outFileName HMEC_H3K27ac.tab.gz
# computeMatrix scale-regions -S E119-H3K27me3.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H3K27me3 \
                            # --outFileName HMEC_H3K27me3.tab.gz
# computeMatrix scale-regions -S E119-H3K36me3.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H3K36me3 \
                            # --outFileName HMEC_H3K36me3.tab.gz
# computeMatrix scale-regions -S E119-H3K79me2.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H3K79me2 \
                            # --outFileName HMEC_H3K79me2.tab.gz
# computeMatrix scale-regions -S E119-H2A.Z.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H2AZ \
                            # --outFileName HMEC_H2AZ.tab.gz
# computeMatrix scale-regions -S E119-H3K4me1.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H3K4me1 \
                            # --outFileName HMEC_H3K4me1.tab.gz
# computeMatrix scale-regions -S E119-H3K4me3.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H3K4me3 \
                            # --outFileName HMEC_H3K4me3.tab.gz
# computeMatrix scale-regions -S E119-H3K9me3.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H3K9me3 \
                            # --outFileName HMEC_H3K9me3.tab.gz
# computeMatrix scale-regions -S E119-H3K9ac.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H3K9ac \
                            # --outFileName HMEC_H3K9ac.tab.gz
# computeMatrix scale-regions -S E119-H4K8ac.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H4K8ac \
                            # --outFileName HMEC_H4K8ac.tab.gz
# computeMatrix scale-regions -S E119-H4K20me1.imputed.pval.signal.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_H4K20me1 \
                            # --outFileName HMEC_H4K20me1.tab.gz
# computeMatrix scale-regions -S HMEC_mH2A1_FE.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_mH2A1 \
                            # --outFileName HMEC_mH2A1.tab.gz
# computeMatrix scale-regions -S HMEC_mH2A2_FE.bigwig \
                            # -R E119_25_imputed12marks_stateno.bed \
                            # --binSize 50 -m 50 -p 16 \
                            # --outFileNameMatrix HMEC_mH2A2 \
                            # --outFileName HMEC_mH2A2.tab.gz

# Consolidate all data 
marks <- c("H2AZ", "H2BK120ac", "H2BK12ac", "H2BK5ac", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K79me2", "H3K9ac", "H3K9me3", "H4K20me1", "H4K8ac", "mH2A1", "mH2A2")

read_mark <- function(mark){
    tab <- read.table(gzfile(mark), skip=1, sep="\t")
    return(tab %>% select(V7))
}

HMEC1 <- lapply(paste0("HMEC_", marks, ".tab.gz"), read_mark)
HMEC1 <- bind_cols(HMEC1) %>% set_colnames(marks)
HMEC_bed <- read.table("E119_25_imputed12marks_stateno.bed", sep="\t")
HMEC1 <- HMEC1 %>% bind_cols(HMEC_bed %>% select(V4)) %>% set_colnames(c(marks, "state"))
write.table(HMEC1, file = "HMEC.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

NHM1 <- lapply(paste0("NHM_", marks, ".tab.gz"), read_mark)
NHM1 <- bind_cols(NHM1) %>% set_colnames(marks)
NHM_bed <- read.table("E059_25_imputed12marks_stateno.bed", sep="\t")
NHM1 <- NHM1 %>% bind_cols(NHM_bed %>% select(V4)) %>% set_colnames(c(marks, "state"))
write.table(NHM1, file = "NHM.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

HepG2_1 <- lapply(paste0("HepG2_", marks, ".tab.gz"), read_mark)
HepG2_1 <- bind_cols(HepG2_1) %>% set_colnames(marks)
HepG2_bed <- read.table("E118_25_imputed12marks_stateno.bed", sep="\t")
HepG2_1 <- HepG2_1 %>% bind_cols(HepG2_bed %>% select(V4)) %>% set_colnames(c(marks, "state"))
write.table(HepG2_1, file = "HepG2.tab", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Read consolidated data
marks <- c("H2AZ", "H2BK120ac", "H2BK12ac", "H2BK5ac", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K79me2", "H3K9ac", "H3K9me3", "H4K20me1", "H4K8ac", "mH2A1", "mH2A2")
HMEC_bed <- read.table("E119_25_imputed12marks_stateno.bed", sep="\t")
HMEC1 <- read.table(file = "HMEC.tab", sep = "\t", skip = 1) %>% set_colnames(c(marks, "state"))
HMEC1 <- HMEC1 %>% mutate(state = as.integer(state)) %>% mutate(state = factor(state, levels = 1:25))
NHM_bed <- read.table("E059_25_imputed12marks_stateno.bed", sep="\t")
NHM1 <- read.table(file = "NHM.tab", sep = "\t", skip = 1) %>% set_colnames(c(marks, "state"))
NHM1 <- NHM1 %>% mutate(state = as.integer(state)) %>% mutate(state = factor(state, levels = 1:25))
HepG2_bed <- read.table("E118_25_imputed12marks_stateno.bed", sep="\t")
HepG2_1 <- read.table(file = "HepG2.tab", sep = "\t", skip = 1) %>% set_colnames(c(marks, "state"))
HepG2_1 <- HepG2_1 %>% mutate(state = as.integer(state)) %>% mutate(state = factor(state, levels = 1:25))
HMEC_bed <- HMEC_bed %>% select(V1:V4) %>% set_colnames(c("chr", "start", "end", "state")) %>% mutate(len = end-start) %>% 
    mutate(state = as.integer(state)) %>% mutate(state = factor(state, levels = 1:25))
HMEC_states <- HMEC_bed %>% select(state, len) %>% group_by(state) %>% summarize_all(sum)
NHM_bed <- NHM_bed %>% select(V1:V4) %>% set_colnames(c("chr", "start", "end", "state")) %>% mutate(len = end-start) %>% 
    mutate(state = as.integer(state)) %>% mutate(state = factor(state, levels = 1:25))
NHM_states <- NHM_bed %>% select(state, len) %>% group_by(state) %>% summarize_all(sum)
HepG2_bed <- HepG2_bed %>% select(V1:V4) %>% set_colnames(c("chr", "start", "end", "state")) %>% mutate(len = end-start) %>% 
    mutate(state = as.integer(state)) %>% mutate(state = factor(state, levels = 1:25))
HepG2_states <- HepG2_bed %>% select(state, len) %>% group_by(state) %>% summarize_all(sum)

# Mann-Whitney test with Bonferroni correction to test for significance of enrichment
combos <- expand.grid(as.character(1:25), marks)
combos_split <- split(combos, seq(nrow(combos)))
get_tests <- function(combo){
    print(paste0(as.character(combo$Var1), " ", as.character(combo$Var2)))
    x <- curr %>% filter(state == as.integer(as.character(combo$Var1))) %>% select(combo$Var2) %>% unlist() %>% unname()
    y <- curr %>% filter(state != as.integer(as.character(combo$Var1))) %>% select(combo$Var2) %>% unlist() %>% unname()
    p <- wilcox.test(x = x, y = y)$p.value
    padj <- p * 400 # Bonferroni correction
    return(padj >= 0.05)
}
curr <- HMEC1
HMEC_insignif <- unname(unlist(lapply(combos_split, get_tests)))
curr <- NHM1
NHM_insignif <- unname(unlist(lapply(combos_split, get_tests)))
curr <- HepG2_1
HepG2_insignif <- unname(unlist(lapply(combos_split, get_tests)))

insignif <- bind_rows(list(HMEC = combos[HMEC_insignif,], NHM = combos[NHM_insignif,], HepG2 = combos[HepG2_insignif,]), .id = "cell")
write.table(insignif, file = "fig1a_stats.tab", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Plot heatmaps
notsigtab <- read.table("fig1a_stats.tab")
HMEC1_medians <- HMEC1 %>% group_by(state) %>% summarize_all(median, na.rm = TRUE)
pdf("fig1a-1.pdf", width=8, height=10)
notsig <- notsigtab %>% 
    set_colnames(c("CellLine", "state", "mark")) %>% 
    filter(CellLine == "HMEC") %>% 
    select(-CellLine) %>% 
    mutate(notsig = "o", state = factor(state, levels = 1:25))
df <- HMEC1_medians %>% pivot_longer(!state, names_to="mark", values_to="median") %>% left_join(notsig) %>% mutate(notsig = replace_na(notsig, ""))
m <- df %>% select(-notsig) %>% pivot_wider(names_from="mark", values_from="median") %>% select(-state) %>% as.matrix()
mp <- df %>% select(-median) %>% pivot_wider(names_from="mark", values_from="notsig") %>% select(-state) %>% as.matrix()
heatmap.2(m, col=bluered(256), scale='column',key=T,keysize=1, trace="none", cellnote = mp, notecol="black", cexRow=1, cexCol=1,symkey=T,symbreaks=T, Colv=TRUE, Rowv = FALSE, dendrogram = "column", margins = c(6, 10), labRow = paste0(1:25, ", len = ", prettyNum(HMEC_states$len, format="d", big.mark=",")))
dev.off()

NHM1_medians <- NHM1 %>% group_by(state) %>% summarize_all(median, na.rm = TRUE)
pdf("fig1a-2.pdf", width=8, height=10)
notsig <- notsigtab %>% 
    set_colnames(c("CellLine", "state", "mark")) %>% 
    filter(CellLine == "NHM") %>% 
    select(-CellLine) %>% 
    mutate(notsig = "o", state = factor(state, levels = 1:25))
df <- NHM1_medians %>% pivot_longer(!state, names_to="mark", values_to="median") %>% left_join(notsig) %>% mutate(notsig = replace_na(notsig, ""))
m <- df %>% select(-notsig) %>% pivot_wider(names_from="mark", values_from="median") %>% select(-state) %>% as.matrix()
mp <- df %>% select(-median) %>% pivot_wider(names_from="mark", values_from="notsig") %>% select(-state) %>% as.matrix()
heatmap.2(m, col=bluered(256), scale='column',key=T,keysize=1, trace="none", cellnote = mp, notecol="black", cexRow=1, cexCol=1,symkey=T,symbreaks=T, Colv=TRUE, Rowv = FALSE, dendrogram = "column", margins = c(6, 10), labRow = paste0(1:25, ", len = ", prettyNum(NHM_states$len, format="d", big.mark=",")))
dev.off()

HepG2_1_medians <- HepG2_1 %>% group_by(state) %>% summarize_all(median, na.rm = TRUE)
pdf("fig1a-3.pdf", width=8, height=10)
notsig <- notsigtab %>% 
    set_colnames(c("CellLine", "state", "mark")) %>% 
    filter(CellLine == "HepG2") %>% 
    select(-CellLine) %>% 
    mutate(notsig = "o", state = factor(state, levels = 1:25))
df <- HepG2_1_medians %>% pivot_longer(!state, names_to="mark", values_to="median") %>% left_join(notsig) %>% mutate(notsig = replace_na(notsig, ""))
m <- df %>% select(-notsig) %>% pivot_wider(names_from="mark", values_from="median") %>% select(-state) %>% as.matrix()
mp <- df %>% select(-median) %>% pivot_wider(names_from="mark", values_from="notsig") %>% select(-state) %>% as.matrix()
heatmap.2(m, col=bluered(256), scale='column',key=T,keysize=1, trace="none", cellnote = mp, notecol="black", cexRow=1, cexCol=1,symkey=T,symbreaks=T, Colv=TRUE, Rowv = FALSE, dendrogram = "column", margins = c(6, 10), labRow = paste0(1:25, ", len = ", prettyNum(HepG2_states$len, format="d", big.mark=",")))
dev.off()

