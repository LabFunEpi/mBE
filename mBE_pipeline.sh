. $1

bedtools intersect -a ${ATAC_peaks} -b ${H3K4me1_peaks} -u > ${outdir}/temp.bed
bedtools subtract -a ${outdir}/temp.bed -b ${blacklist} -A > ${outdir}/candidate_enhancers.bed
rm -f ${outdir}/temp.bed

awk '{print $1"\t" int($2+($3-$2)/2) "\t" int($2+($3-$2)/2)+1 "\t" $4}' ${outdir}/candidate_enhancers.bed > ${outdir}/candidate_enhancer_centers.bed
bedtools intersect -loj -a ${outdir}/candidate_enhancer_centers.bed -b ${ENCODE_cCRE} > ${outdir}/candidates_wENCODE.bed

ab=1000
cd ${outdir}

computeMatrix reference-point -S ${H3K27ac} \
                              -R ${outdir}/candidate_enhancer_centers.bed \
                              --referencePoint center \
                              -b ${ab} -a ${ab} -p 16 \
                              --outFileNameMatrix H3K27ac \
                              --outFileName H3K27ac.tab.gz
computeMatrix reference-point -S ${H3K4me1} \
                              -R ${outdir}/candidate_enhancer_centers.bed \
                              --referencePoint center \
                              -b ${ab} -a ${ab} -p 16 \
                              --outFileNameMatrix H3K4me1 \
                              --outFileName H3K4me1.tab.gz
computeMatrix reference-point -S ${H3K4me3} \
                              -R ${outdir}/candidate_enhancer_centers.bed \
                              --referencePoint center \
                              -b ${ab} -a ${ab} -p 16 \
                              --outFileNameMatrix H3K4me3 \
                              --outFileName H3K4me3.tab.gz
computeMatrix reference-point -S ${H2Az} \
                              -R ${outdir}/candidate_enhancer_centers.bed \
                              --referencePoint center \
                              -b ${ab} -a ${ab} -p 16 \
                              --outFileNameMatrix H2Az \
                              --outFileName H2Az.tab.gz
computeMatrix reference-point -S ${H3K27me3} \
                              -R ${outdir}/candidate_enhancer_centers.bed \
                              --referencePoint center \
                              -b ${ab} -a ${ab} -p 16 \
                              --outFileNameMatrix H3K27me3 \
                              --outFileName H3K27me3.tab.gz
computeMatrix reference-point -S ${mH2A1} \
                              -R ${outdir}/candidate_enhancer_centers.bed \
                              --referencePoint center \
                              -b ${ab} -a ${ab} -p 16 \
                              --outFileNameMatrix mH2A1 \
                              --outFileName mH2A1.tab.gz
computeMatrix reference-point -S ${mH2A2} \
                              -R ${outdir}/candidate_enhancer_centers.bed \
                              --referencePoint center \
                              -b ${ab} -a ${ab} -p 16 \
                              --outFileNameMatrix mH2A2 \
                              --outFileName mH2A2.tab.gz
computeMatrix reference-point -S ${CTCF} \
                              -R ${outdir}/candidate_enhancer_centers.bed \
                              --referencePoint center \
                              -b ${ab} -a ${ab} -p 16 \
                              --outFileNameMatrix CTCF \
                              --outFileName CTCF.tab.gz

Rscript cluster.R ${outdir}/candidate_enhancers.bed ${outdir}/candidates_wENCODE.bed ${outdir} ${k}




