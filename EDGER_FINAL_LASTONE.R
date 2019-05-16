library(readr)
E12_5_E15_5_WT_HC <- read_delim("~/Lhx2_RNASeq/COUNTS/E15.5_WT_HC.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)

y <- DGEList(counts=E15_5_WT_Lhx2_mut_HC[,7:10], genes=E15_5_WT_Lhx2_mut_HC[,1:2],group = factor(c("E15_5_WT_HC","E15_5_WT_HC","E15_5_Lhx2_mut_HC","E15_5_Lhx2_mut_HC")))

E15_5_WT_Lhx2_mut_ctx
library(org.Mm.eg.db)
idfound <- y$genes$Geneid %in% mappedRkeys(org.Mm.egSYMBOL)
y <- y[idfound,]

countsPerMillion <- cpm(y)
summary(countsPerMillion)
countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
y1 <- y[keep,]
summary(cpm(y1))

y2 <- calcNormFactors(y1, method="TMM")
d1 <- estimateCommonDisp(y2, verbose=T)
d1 <- estimateTagwiseDisp(d1)


design.mat <- model.matrix(~ 0 + y2$samples$group)
colnames(design.mat) <- levels(y2$samples$group)
d2 <- estimateGLMCommonDisp(y2,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

et = exactTest(d1, pair=c("E15_5_WT_HC", "E15_5_Lhx2_mut_HC"))
tTags = topTags(et,n=NULL)

is.de = decideTestsDGE(et, p.value = 0.05)
summary(is.de)
plotMD(et)
abline(h=c(-1, 1), col="blue")

sig_DE=tTags[tTags$table$FDR<0.01,]
write.table(sig_DE, file='sig_DE', sep='\t', quote=F, row.names=F)

up_sig_DE=sig_DE[sig_DE$table$logFC>0.58,]
write.table(up_sig_DE, file='up_sig_DE', sep='\t', quote=F, row.names=F)

down_sig_DE=sig_DE[sig_DE$table$logFC<(-0.58),]
write.table(down_sig_DE, file='down_sig_DE', sep='\t', quote=F, row.names=F)

up_dge_list = up_sig_DE_E12.5_HC_WT_MUT$table$Geneid
down_dge_list = down_sig_DE_E12.5_HC_WT_MUT$table$Geneid
DGE = c(up_dge_list,down_dge_list)

newdata <- subset(sig_DE, logFC >= 5 | logFC < -5)

eg <- bitr(DGE, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)  
eg
gene = as.vector(eg$ENTREZID)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Mm.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

ego <- enrichGO(gene          = gene,
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

dotplot(ego)
enrichMap(ego)
plotGOgraph(ego)

library(topGO)
go = goana(et)
topGO(go,ontology ="BP",sort = "Up")






dba_test = dba.peakset(DBA = NULL, peaks = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E10.5/peaks-rep1/E10.5_INPUT_LHX2_MUT_REP1_peaks.bed",sampID = "e10.5_1",factor = "LHX2",replicate = 1,peak.format = "bed", bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E10.5/E10.5_LHX2_IP_REP1/E10.5_LHX2_IP_REP1.bam", bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E10.5/E10.5_LHX2_IP_REP1/E10.5_LHX2_IP_REP1.bam",condition = "E10.5")

dba_test = dba.peakset(DBA = dba_test, peaks = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E10.5/peaks-rep2/E10.5_INPUT_LHX2_MUT_REP2_peaks.bed",sampID = "e10.5_2",factor = "LHX2",replicate = 2,peak.format = "bed", bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E10.5/E10.5_LHX2_IP_REP2/E10.5_LHX2_IP_REP2.bam", bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E10.5/E10.5_LHX2_IP_REP2/E10.5_LHX2_IP_REP2.bam", condition = "E10.5")

dba_test = dba.peakset(DBA = dba_test,peaks = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/CTX/E12.5_CTX_peaks1/E12.5_CTX_INPUT_LHX2_REP1_peaks.bed", sampID = "e12.5_1", factor = "LHX2", replicate = 1, peak.format = "bed",bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/CTX/E12.5_CTX_LHX2_REP1/E12.5_CTX_LHX2_REP1.bam",bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/CTX/E12.5_CTX_INPUT_REP1/E12.5_CTX_INPUT_REP1.bam",condition = "E12.5_ctx")

dba_test = dba.peakset(DBA = dba_test, peaks = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/CTX/e12.5_cortes_peak2/E12.5_CTX_INPUT_LHX2_REP2_peaks.bed",sampID = "e12.5_2",factor = "LHX2",replicate = 2, peak.format = "bed", bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/CTX/E12.5_CTX_LHX2_REP2/E12.5_CTX_LHX2_REP2.bam",bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/CTX/E12.5_CTX_INPUT_REP2/E12.5_CTX_INPUT_REP2.bam",condition = "E12.5_ctx" )

dba_test = dba.peakset(DBA = dba_test,peaks = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/HC/e12.5-hc-peak-rep1/E12.5_HC_INPUT_LHX2_REP1_peaks.bed",sampID = "e12.5_HC_1",factor = "LHX2",replicate = 1,peak.format = "bed", bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/HC/E12.5_HC_LHX2_REP1/E12.5_HC_LHX2_REP1.bam",bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/HC/E12.5_HC_INPUT_REP1/E12.5_HC_INPUT_REP1.bam", condition = "E12.5_hc")
dba_test = dba.peakset(DBA = dba_test,peaks = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/HC/e12.5-hc-peak-rep2/E12.5_HC_INPUT_LHX2_REP2_peaks.bed",sampID = "e12.5_HC_2",factor = "LHX2",replicate = 2,peak.format = "bed", bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/HC/E12.5_HC_LHX2_REP2/E12.5_HC_LHX2_REP2.bam",bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E12.5/HC/E12.5_HC_INPUT_REP2/E12.5_HC_INPUT_REP2.bam", condition = "E12.5_hc")

dba_test = dba.peakset(DBA = dba_test,peaks = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/CTX/e15.5-ctx-peaks-rep1/E15.5_CTX_INPUT_LHX2_REP1_peaks.bed",sampID = "e15.5_ctx_1",factor = "LHX2", replicate = 1, peak.format = "bed",bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/CTX/E15.5_CTX_LHX2_REP1/E15.5_CTX_LHX2_REP1.bam", bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/CTX/E15.5_CTX_INPUT_REP1/E15.5_CTX_INPUT_REP1.bam",condition = "e15.5_ctx")
dba_test = dba.peakset(DBA = dba_test,peaks = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/CTX/e15.5-ctx-peaks-rep2/E15.5_CTX_INPUT_LHX2_REP2_peaks.bed",sampID = "e15.5ctx_2",factor = "LHX2",replicate = 2 ,peak.format = "bed",bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/CTX/E15.5_CTX_LHX2_REP2/E15.5_CTX_LHX2_REP2.bam",bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/CTX/E15.5_CTX_INPUT_REP2/E15.5_CTX_INPUT_REP2.bam",condition = "e15.5_ctx")

dba_test = dba.peakset(DBA = dba_test,peaks = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/HC/e15.5-hc-peaks-rep1/E15_5_HC_INPUT_LHX2_REP1_peaks.bed",sampID = "e15.5_hc_1",factor = "LHX2",replicate = 1,peak.format = "bed",bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/HC/E15.5_HC_LHX2_REP1/E15.5_HC_LHX2_REP1.bam",bamControl ="/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/HC/E15.5_HC_INPUT_REP1/E15.5_HC_INPUT_REP1.bam", condition = "e15.5_hc")
dba_test = dba.peakset(DBA = dba_test, peaks = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/HC/e15.5-hc-peaks-rep2/E15.5_HC_INPUT_LHX2_REP2_peaks.bed",sampID = "e15.5_hc_2", factor = "LHX2",replicate = 2, peak.format = "bed", bamReads = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/HC/e15.5-hc-peaks-rep2/E15.5_HC_INPUT_LHX2_REP2_peaks.bed", bamControl = "/media/darwin/STOLE_HD1/1704KHP-0068_Shubha_Saurabh_Mauli_Satya_Rutika/report/Lhx2_ChIPseq/E15.5/HC/E15.5_HC_INPUT_REP2/E15.5_HC_INPUT_REP2.bam" ,condition = "e15.5_hc")
















