COMP = c("E12.5_CTX_WT_MUT","E12.5_HC_WT_MUT","E15.5_HC_WT_MUT","E15.5_CTX_WT_MUT","E12.5_E15.5_CTX","E12.5_E15.5_HC","E12.5_CTX_E12.5_HC","E15.5_CTX_E15.5_HC")

et_tags = matrix(ncol = 1,nrow = 8)
for(i in 1:5){
  print(tags_index = topTags(COMP[[i]],n=NULL))
  
  print(file_location =  paste("/home/darwin/NEW_LHX2_RNA_SEQ/DEG", COMP[[i]], sep = "/") )
  print(txt_file_name = paste(file_location, ".txt", sep = ""))
  
  print(write.table(up_deg, file =  txt_file_name ,row.names = FALSE, sep = "\t"))
  
}

print(up_deg = tags_index[[i]][tags_index$table$FDR<0.01,])



PEAK1_GR = GRanges(E12_PEAK1$chr, IRanges(E12_PEAK1$start,E12_PEAK1$end),strand = "*")
PEAK2_GR = GRanges(E12_PEAK2$chr, IRanges(E12_PEAK2$start,E12_PEAK2$end),strand = "*")

COV_PLOT1 = covplot(PEAK1_GR)
COV_PLOT2 = covplot(PEAK2_GR)


PROMTOER1 = getPromoters(TxDb = txbd , upstream = 3000, downstream = 3000)
TAGMATRIX1 = getTagMatrix(PEAK1_GR, windows = PROMTOER1)
PROMTOER2 = getPromoters(TxDb = txbd , upstream = 3000, downstream = 3000)
TAGMATRIX2 = getTagMatrix(PEAK2_GR, windows = PROMTOER1)

tagheatMAP1 = tagHeatmap(TAGMATRIX1, xlim=c(-3000, 3000), color="red")
tagheatMAP2 = tagHeatmap(TAGMATRIX2, xlim=c(-3000, 3000), color="red")

plotAvgProf(TAGMATRIX1, xlim=c(-3000, 3000))

ANNO_PEAK1 = annotatePeak(PEAK1_GR, tssRegion=c(-3000, 3000),TxDb = txbd, annoDb="org.Mm.eg.db")
ANNO_PEAK2 = annotatePeak(PEAK2_GR, tssRegion=c(-3000, 3000),TxDb = txbd, annoDb="org.Mm.eg.db")

PEAK1_pie = plotAnnoPie(ANNO_PEAK1)
PEAK2_pie = plotAnnoPie(ANNO_PEAK2)

library(ReactomePA)
PATHWAY1 = enrichPathway(as.data.frame(ANNO_PEAK1@anno)$geneId)


E12.5_CTX = dba.peakset(DBA = NULL , peaks = "/home/darwin/NEW_CHIP/E12.5_CTX_PEAKS_REP1"
                                                , peak.caller = "raw", sampID = "E12.5_CTX",replicate = 1)
E12.5_CTX = dba.peakset(DBA = E12.5_CTX , peaks = "/home/darwin/NEW_CHIP/E12.5_CTX_PEAKS_REP2"
                        , peak.caller = "raw", sampID = "E12.5_CTX",replicate = 2)




plot(E12.5_CTX)
test = dba.count(E12.5_CTX)


peakfiles = 
  peaksfiles

datadir = "/home/darwin/NEW_CHIP/"
peakfiles = dir(file.path(getwd(), datadir), pattern = "*.txt", full.names = TRUE)
names(peakfiles) = gsub(".txt","",basename(txtfiles))



##DIFFFBIND

  
dba_test = dba.peakset(DBA = NULL, peaks = "/home/darwin/NEW_CHIP/DIFFBIND_TEST/E10.5_1.txt"
                       ,sampID = "e10.5_1",factor = "LHX2",replicate = 1,peak.format = "raw"
                       , bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E10.5/E10.5_LHX2_IP_REP1/E10.5_LHX2_IP_REP1.bam"
                       , bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E10.5/E10.5_LHX2_IP_REP1/E10.5_LHX2_IP_REP1.bam",condition = "E10.5")

dba_test = dba.peakset(DBA = dba_test , peaks = "/home/darwin/NEW_CHIP/DIFFBIND_TEST/E10.5_2.txt"
                       ,sampID = "e10.5_2",factor = "LHX2",replicate = 2,peak.format = "raw"
                       , bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E10.5/E10.5_LHX2_IP_REP2/E10.5_LHX2_IP_REP2.bam"
                       , bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E10.5/E10.5_LHX2_IP_REP2/E10.5_LHX2_IP_REP2.bam", condition = "E10.5")

dba_test = dba.peakset(DBA = dba_test ,peaks = "/home/darwin/NEW_CHIP/DIFFBIND_TEST/E12.5_CTX_1.txt"
                       , sampID = "e12.5_ctx_1", factor = "LHX2", replicate = 1, peak.format = "raw"
                       ,bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/CTX/E12.5_CTX_LHX2_REP1/E12.5_CTX_LHX2_REP1.bam"
                       ,bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/CTX/E12.5_CTX_INPUT_REP1/E12.5_CTX_INPUT_REP1.bam",condition = "E12.5_ctx")

dba_test = dba.peakset(DBA = dba_test, peaks = "/home/darwin/NEW_CHIP/DIFFBIND_TEST/E12.5_CTX_2.txt"
                       ,sampID = "e12.5_ctx_2",factor = "LHX2",replicate = 2, peak.format = "raw", 
                       bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/CTX/E12.5_CTX_LHX2_REP2/E12.5_CTX_LHX2_REP2.bam"
                       ,bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/CTX/E12.5_CTX_INPUT_REP2/E12.5_CTX_INPUT_REP2.bam",condition = "E12.5_ctx" )

dba_test = dba.peakset(DBA = dba_test,peaks = "/home/darwin/NEW_CHIP/DIFFBIND_TEST/E12.5_HC_1.txt"
                       ,sampID = "e12.5_HC_1",factor = "LHX2",replicate = 1,peak.format = "raw"
                       , bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/HC/E12.5_HC_LHX2_REP1/E12.5_HC_LHX2_REP1.bam"
                       ,bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/HC/E12.5_HC_INPUT_REP1/E12.5_HC_INPUT_REP1.bam", condition = "E12.5_hc")

dba_test = dba.peakset(DBA = dba_test,peaks = "/home/darwin/NEW_CHIP/DIFFBIND_TEST/E12.5_HC_2.txt"
                       ,sampID = "e12.5_HC_2",factor = "LHX2",replicate = 2,peak.format = "raw"
                       , bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/HC/E12.5_HC_LHX2_REP2/E12.5_HC_LHX2_REP2.bam"
                       ,bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/HC/E12.5_HC_INPUT_REP2/E12.5_HC_INPUT_REP2.bam", condition = "E12.5_hc")

dba_test = dba.peakset(DBA = dba_test,peaks = "/home/darwin/NEW_CHIP/DIFFBIND_TEST/E15.5_CTX_1.txt"
                       ,sampID = "e15.5_ctx_1",factor = "LHX2", replicate = 1, peak.format = "raw"
                       ,bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/CTX/E15.5_CTX_LHX2_REP1_repeated/E15.5_CTX_LHX2_REP1.bam"
                       , bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/CTX/E15.5_CTX_INPUT_REP1_repeated/E15.5_CTX_INPUT_REP1.bam",condition = "E15.5_ctx")

dba_test = dba.peakset(DBA = dba_test,peaks = "/home/darwin/NEW_CHIP/DIFFBIND_TEST/E15.5_CTX_2.txt"
                       ,sampID = "e15.5_ctx_2",factor = "LHX2",replicate = 2 ,peak.format = "raw"
                       ,bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/CTX/E15.5_CTX_LHX2_REP2/E15.5_CTX_LHX2_REP2.bam"
                       ,bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/CTX/E15.5_CTX_INPUT_REP2/E15.5_CTX_INPUT_REP2.bam",condition = "E15.5_ctx")

dba_test = dba.peakset(DBA = dba_test,peaks = "/home/darwin/NEW_CHIP/DIFFBIND_TEST/E15.5_HC_1.txt"
                       ,sampID = "e15.5_hc_1",factor = "LHX2",replicate = 1,peak.format = "raw"
                       ,bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/HC/E15.5_HC_LHX2_REP1/E15.5_HC_LHX2_REP1.bam"
                       ,bamControl ="/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/HC/E15.5_HC_INPUT_REP1/E15.5_HC_INPUT_REP1.bam", condition = "E15.5_hc")

dba_test = dba.peakset(DBA = dba_test, peaks = "/home/darwin/NEW_CHIP/DIFFBIND_TEST/E15.5_HC_2.txt"
                       ,sampID = "e15.5_hc_2", factor = "LHX2",replicate = 2, peak.format = "raw"
                       ,bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/HC/e15.5-hc-peaks-rep2/E15.5_HC_INPUT_LHX2_REP2_peaks.bed"
                       , bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/HC/E15.5_HC_INPUT_REP2/E15.5_HC_INPUT_REP2.bam" ,condition = "E15.5_hc")


test = dba.count(dba_test)
test2 = dba.contrast(test, categories = DBA_CONDITION)
test = dba.analyze(test)
















































