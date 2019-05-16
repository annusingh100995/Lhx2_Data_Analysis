y <- DGEList(counts=LHX2_RNA_COUNT[,7:ncol(LHX2_RNA_COUNT)], genes=LHX2_RNA_COUNT[,1:2],group = factor(c( "E12_5_CTX_LHX2_MUT", "E12_5_CTX_LHX2_MUT", "E12_5_CTX_LHX2_MUT",
                                                                                                          "E12_5_CTX_WT", "E12_5_CTX_WT", "E12_5_CTX_WT",
                                                                                                          "E12_5_HC_LHX2_MUT", "E12_5_HC_LHX2_MUT", "E12_5_HC_LHX2_MUT",
                                                                                                          "E12_5_HC_WT", "E12_5_HC_WT", "E12_5_HC_WT",
                                                                                                        "E15_5_CTX_LHX2_MUT","E15_5_CTX_LHX2_MUT","E15_5_CTX_LHX2_MUT",
                                                                                                         "E15_5_CTX_WT","E15_5_CTX_WT","E15_5_CTX_WT",
                                                                                                         "E15_5_HC_LHX2_MUT","E15_5_HC_LHX2_MUT","E15_5_HC_LHX2_MUT",
                                                                                                         "E15_5_HC_WT","E15_5_HC_WT","E15_5_HC_WT")))


countsPerMillion <- cpm(y)
summary(countsPerMillion)
countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
y1_cpm <- y[keep,]
summary(cpm(y1_cpm))

Count_LHX2_RNA_2 =  LHX2_RNA_COUNT[,-c(13,17)]

Y_LHX2 = DGEList(counts=Count_LHX2_RNA_2[,7:ncol(Count_LHX2_RNA_2)], genes=Count_LHX2_RNA_2[,1:2],group = factor(c( "E12_5_CTX_LHX2_MUT", "E12_5_CTX_LHX2_MUT", "E12_5_CTX_LHX2_MUT",
                                                                                                             "E12_5_CTX_WT", "E12_5_CTX_WT", "E12_5_CTX_WT",
                                                                                                             "E12_5_HC_LHX2_MUT", "E12_5_HC_LHX2_MUT",
                                                                                                             "E12_5_HC_WT", "E12_5_HC_WT",
                                                                                                             "E15_5_CTX_LHX2_MUT","E15_5_CTX_LHX2_MUT","E15_5_CTX_LHX2_MUT",
                                                                                                             "E15_5_CTX_WT","E15_5_CTX_WT","E15_5_CTX_WT",
                                                                                                             "E15_5_HC_LHX2_MUT","E15_5_HC_LHX2_MUT","E15_5_HC_LHX2_MUT",
                                                                                                             "E15_5_HC_WT","E15_5_HC_WT","E15_5_HC_WT")))

plotMDS(Y_LHX2)
countsPerMillion <- cpm(Y_LHX2)
summary(countsPerMillion)
countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
Y_LHX2_cpm <- Y_LHX2[keep,]
summary(cpm(Y_LHX2_cpm))

#########################################################################################################

y2 <- calcNormFactors(Y_LHX2_cpm, method="TMM")
d1 <- estimateCommonDisp(y2, verbose=T)
d1 <- estimateTagwiseDisp(d1)


design.mat <- model.matrix(~ 0 + y2$samples$group)
colnames(design.mat) <- levels(y2$samples$group)
d2 <- estimateGLMCommonDisp(y2,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

############################################################################################################

#### 1 E12.5_CTX_WT_MUT

et_E12.5_CTX_WT_MUT = exactTest(d1, pair=c("E12_5_CTX_WT","E12_5_CTX_LHX2_MUT"))
tTags_E12.5_CTX_WT_MUT = topTags(et_E12.5_CTX_WT_MUT,n=NULL)

is.de_E12.5_CTX_WT_MUT = decideTestsDGE(et_E12.5_CTX_WT_MUT, p.value = 0.05)
summary(is.de_E12.5_CTX_WT_MUT)
a = plotMD(et_E12.5_CTX_WT_MUT)
abline(h=c(-1, 1), col="blue")

sig_DE_E12.5_CTX_WT_MUT=tTags_E12.5_CTX_WT_MUT[tTags_E12.5_CTX_WT_MUT$table$FDR<0.01,]
write.table(sig_DE_E12.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/sig_DE_E12.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)

up_sig_DE_E12.5_CTX_WT_MUT=sig_DE_E12.5_CTX_WT_MUT[sig_DE_E12.5_CTX_WT_MUT$table$logFC>0.58,]
write.table(up_sig_DE_E12.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/up_sig_DE_E12.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)

down_sig_DE_E12.5_CTX_WT_MUT=sig_DE_E12.5_CTX_WT_MUT[sig_DE_E12.5_CTX_WT_MUT$table$logFC<(-0.58),]
write.table(down_sig_DE_E12.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/down_sig_DE_E12.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)

#################################################################################################
 ######  2 E12.5_HC_WT_MUT

et_E12.5_HC_WT_MUT = exactTest(d1, pair=c( "E12_5_HC_WT","E12_5_HC_LHX2_MUT"))
tTags_E12.5_HC_WT_MUT = topTags(et_E12.5_HC_WT_MUT,n=NULL)

is.de_E12.5_HC_WT_MUT = decideTestsDGE(et_E12.5_HC_WT_MUT, p.value = 0.05)
summary(is.de_E12.5_HC_WT_MUT)
plotMD(et_E12.5_HC_WT_MUT)
abline(h=c(-1, 1), col="blue")

sig_DE_E12.5_HC_WT_MUT=tTags_E12.5_HC_WT_MUT[tTags_E12.5_HC_WT_MUT$table$FDR<0.01,]
write.table(sig_DE_E12.5_HC_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/sig_DE_E12.5_HC_WT_MUT', sep='\t', quote=F, row.names=F)

up_sig_DE_E12.5_HC_WT_MUT=sig_DE_E12.5_HC_WT_MUT[sig_DE_E12.5_HC_WT_MUT$table$logFC>0.58,]
write.table(up_sig_DE_E12.5_HC_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/up_sig_DE_E12.5_HC_WT_MUT', sep='\t', quote=F, row.names=F)

down_sig_DE_E12.5_HC_WT_MUT=sig_DE_E12.5_HC_WT_MUT[sig_DE_E12.5_HC_WT_MUT$table$logFC<(-0.58),]
write.table(down_sig_DE_E12.5_HC_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/down_sig_DE_E12.5_HC_WT_MUT', sep='\t', quote=F, row.names=F)
##############################################################################################
###################3  E15.5_CTX_WT_MUT

et_E15.5_CTX_WT_MUT = exactTest(d1, pair=c("E15_5_CTX_WT","E15_5_CTX_LHX2_MUT"))
tTags_E15.5_CTX_WT_MUT = topTags(et_E15.5_CTX_WT_MUT,n=NULL)

is.de_E15.5_CTX_WT_MUT = decideTestsDGE(et_E15.5_CTX_WT_MUT, p.value = 0.05)
summary(is.de_E15.5_CTX_WT_MUT)
plotMD(et_E15.5_CTX_WT_MUT)
abline(h=c(-1, 1), col="blue")


sig_DE_E15.5_CTX_WT_MUT=tTags_E15.5_CTX_WT_MUT[tTags_E15.5_CTX_WT_MUT$table$FDR<0.01,]
write.table(sig_DE_E15.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/sig_DE_E15.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)

up_sig_DE_E15.5_CTX_WT_MUT=sig_DE_E15.5_CTX_WT_MUT[sig_DE_E15.5_CTX_WT_MUT$table$logFC>0.58,]
write.table(up_sig_DE_E15.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/up_sig_DE_E15.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)

down_sig_DE_E15.5_CTX_WT_MUT=sig_DE_E15.5_CTX_WT_MUT[sig_DE_E15.5_CTX_WT_MUT$table$logFC<(-0.58),]
write.table(down_sig_DE_E15.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/down_sig_DE_E15.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)
####################################################################################################
##################  4 E15.5_HC_LHX2_MUT

et_E15.5_HC_WT_MUT = exactTest(d1, pair=c("E15_5_HC_WT","E15_5_HC_LHX2_MUT"))
tTags_E15.5_HC_WT_MUT = topTags(et_E15.5_HC_WT_MUT,n=NULL)

is.de_E15.5_HC_WT_MUT = decideTestsDGE(et_E15.5_HC_WT_MUT, p.value = 0.05)
summary(is.de_E15.5_HC_WT_MUT)
plotMD(et_E15.5_HC_WT_MUT)
abline(h=c(-1, 1), col="blue")


sig_DE_E15.5_HC_WT_MUT=tTags_E15.5_HC_WT_MUT[tTags_E15.5_HC_WT_MUT$table$FDR<0.01,]
write.table(sig_DE_E15.5_HC_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/sig_DE_E15.5_HC_WT_MUT', sep='\t', quote=F, row.names=F)

up_sig_DE_E15.5_HC_WT_MUT=sig_DE_E15.5_HC_WT_MUT[sig_DE_E15.5_HC_WT_MUT$table$logFC>0.58,]
write.table(up_sig_DE_E15.5_HC_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/up_sig_DE_E15.5_HC_WT_MUT', sep='\t', quote=F, row.names=F)

down_sig_DE_E15.5_HC_WT_MUT=sig_DE_E15.5_HC_WT_MUT[sig_DE_E15.5_HC_WT_MUT$table$logFC<(-0.58),]
write.table(down_sig_DE_E15.5_HC_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/down_sig_DE_E15.5_HC_WT_MUT', sep='\t', quote=F, row.names=F)
#####################################################################################################
########### 5 E12.5_E15.5_CTX

et_E12.5_E15.5_CTX = exactTest(d1, pair=c("E12_5_CTX_LHX2_MUT", "E15_5_CTX_LHX2_MUT"))
tTags_E12.5_E15.5_CTX = topTags(et_E12.5_E15.5_CTX,n=NULL)

is.de_E12.5_E15.5_CTX = decideTestsDGE(et_E12.5_E15.5_CTX, p.value = 0.05)
summary(is.de_E12.5_E15.5_CTX)
plotMD(et_E12.5_E15.5_CTX)
abline(h=c(-1, 1), col="blue")


sig_DE_E12.5_E15.5_CTX=tTags_E12.5_E15.5_CTX[tTags_E12.5_E15.5_CTX$table$FDR<0.01,]
write.table(sig_DE_E12.5_E15.5_CTX, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/sig_DE_E12.5_E15.5_CTX', sep='\t', quote=F, row.names=F)

up_sig_DE_E12.5_E15.5_CTX=sig_DE_E12.5_E15.5_CTX[sig_DE_E12.5_E15.5_CTX$table$logFC>0.58,]
write.table(up_sig_DE_E12.5_E15.5_CTX, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/up_sig_E12.5_E15.5_CTX', sep='\t', quote=F, row.names=F)

down_sig_DE_E12.5_E15.5_CTX=sig_DE_E12.5_E15.5_CTX[sig_DE_E12.5_E15.5_CTX$table$logFC<(-0.58),]
write.table(down_sig_DE_E12.5_E15.5_CTX, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/down_sig_DE_E12.5_E15.5_CTX', sep='\t', quote=F, row.names=F)
################################################################################################
#########6 E12.5_E15.5_HC

et_E12.5_E15.5_HC = exactTest(d1, pair=c("E12_5_HC_LHX2_MUT", "E15_5_HC_LHX2_MUT"))
tTags_E12.5_E15.5_HC = topTags(et_E12.5_E15.5_HC,n=NULL)

is.de_E12.5_E15.5_HC = decideTestsDGE(et_E12.5_E15.5_HC, p.value = 0.05)
summary(is.de_E12.5_E15.5_HC)
plotMD(et_E12.5_E15.5_HC)
abline(h=c(-1, 1), col="blue")


sig_DE_E12.5_E15.5_HC=tTags_E12.5_E15.5_HC[tTags_E12.5_E15.5_HC$table$FDR<0.01,]
write.table(sig_DE_E12.5_E15.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/sig_DE_E12.5_E15.5_HC', sep='\t', quote=F, row.names=F)

up_sig_DE_E12.5_E15.5_HC=sig_DE_E12.5_E15.5_HC[sig_DE_E12.5_E15.5_HC$table$logFC>0.58,]
write.table(up_sig_DE_E12.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/up_sig_DE_E12.5_E15.5_HC', sep='\t', quote=F, row.names=F)

down_sig_DE_E12.5_E15.5_HC=sig_DE_E12.5_E15.5_HC[sig_DE_E12.5_E15.5_HC$table$logFC<(-0.58),]
write.table(down_sig_DE_E12.5_E15.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/down_sig_DE_E12.5_E15.5_HC', sep='\t', quote=F, row.names=F)

##################################################################################################
################7  E12.5_CTX_E12.5_HC

et_E12.5_CTX_E12.5_HC = exactTest(d1, pair=c("E12_5_CTX_LHX2_MUT", "E12_5_HC_LHX2_MUT"))
tTags_E12.5_CTX_E12.5_HC = topTags(et_E12.5_CTX_E12.5_HC,n=NULL)

is.de_E12.5_CTX_E12.5_HC = decideTestsDGE(et_E12.5_CTX_E12.5_HC, p.value = 0.05)
summary(is.de_E12.5_CTX_E12.5_HC)
plotMD(et_E12.5_CTX_E12.5_HC)
abline(h=c(-1, 1), col="blue")


sig_DE_E12.5_CTX_E12.5_HC=tTags_E12.5_CTX_E12.5_HC[tTags_E12.5_CTX_E12.5_HC$table$FDR<0.01,]
write.table(sig_DE_E12.5_CTX_E12.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/sig_DE_E12.5_CTX_E12.5_HC', sep='\t', quote=F, row.names=F)

up_sig_DE_E12.5_CTX_E12.5_HC=sig_DE_E12.5_CTX_E12.5_HC[sig_DE_E12.5_CTX_E12.5_HC$table$logFC>0.58,]
write.table(up_sig_DE_E12.5_CTX_E12.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/up_sig_DE_E12.5_CTX_E12.5_HC', sep='\t', quote=F, row.names=F)

down_sig_DE_E12.5_CTX_E12.5_HC=sig_DE_E12.5_CTX_E12.5_HC[sig_DE_E12.5_CTX_E12.5_HC$table$logFC<(-0.58),]
write.table(down_sig_DE_E12.5_CTX_E12.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/down_sig_DE_E12.5_CTX_E12.5_HC', sep='\t', quote=F, row.names=F)

###########################################################################################################
################# 8 E15.5_CTX_E15.5_HC
et_E15.5_CTX_E15.5_HC = exactTest(d1, pair=c("E15_5_CTX_LHX2_MUT", "E15_5_HC_LHX2_MUT"))
tTags_E15.5_CTX_E15.5_HC = topTags(et_E15.5_CTX_E15.5_HC,n=NULL)

is.de_E15.5_CTX_E15.5_HC = decideTestsDGE(et_E15.5_CTX_E15.5_HC, p.value = 0.05)
summary(is.de_E15.5_CTX_E15.5_HC)
a = plotMD(et_E15.5_CTX_E15.5_HC)
abline(h=c(-1, 1), col="blue")


sig_DE_E15.5_CTX_E15.5_HC=tTags_E15.5_CTX_E15.5_HC[tTags_E15.5_CTX_E15.5_HC$table$FDR<0.01,]
write.table(sig_DE_E15.5_CTX_E15.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/sig_DE_E15.5_CTX_E15.5_HC', sep='\t', quote=F, row.names=F)

up_sig_DE_E15.5_CTX_E15.5_HC=sig_DE_E15.5_CTX_E15.5_HC[sig_DE_E15.5_CTX_E15.5_HC$table$logFC>0.58,]
write.table(up_sig_DE_E15.5_CTX_E15.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/up_sig_DE_E15.5_CTX_E15.5_HC', sep='\t', quote=F, row.names=F)

down_sig_DE_E15.5_CTX_E15.5_HC=sig_DE_E15.5_CTX_E15.5_HC[sig_DE_E15.5_CTX_E15.5_HC$table$logFC<(-0.58),]
write.table(down_sig_DE_E15.5_CTX_E15.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG_P_0.05/down_sig_DE_E15.5_CTX_E15.5_HC', sep='\t', quote=F, row.names=F)

############################################################################################################################
p == 0.1

#### 1 E12.5_CTX_WT_MUT

et_E12.5_CTX_WT_MUT = exactTest(d1, pair=c("E12_5_CTX_WT","E12_5_CTX_LHX2_MUT"))
tTags_E12.5_CTX_WT_MUT = topTags(et_E12.5_CTX_WT_MUT,n=NULL)

is.de_E12.5_CTX_WT_MUT = decideTestsDGE(et_E12.5_CTX_WT_MUT, p.value = 0.1)
summary(is.de_E12.5_CTX_WT_MUT)
a = plotMD(et_E12.5_CTX_WT_MUT)
abline(h=c(-1, 1), col="blue")

sig_DE_E12.5_CTX_WT_MUT=tTags_E12.5_CTX_WT_MUT[tTags_E12.5_CTX_WT_MUT$table$FDR<0.01,]
write.table(sig_DE_E12.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/sig_DE_E12.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)

up_sig_DE_E12.5_CTX_WT_MUT=sig_DE_E12.5_CTX_WT_MUT[sig_DE_E12.5_CTX_WT_MUT$table$logFC>0.58,]
write.table(up_sig_DE_E12.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/up_sig_DE_E12.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)

down_sig_DE_E12.5_CTX_WT_MUT=sig_DE_E12.5_CTX_WT_MUT[sig_DE_E12.5_CTX_WT_MUT$table$logFC<(-0.58),]
write.table(down_sig_DE_E12.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/down_sig_DE_E12.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)

#################################################################################################
######  2 E12.5_HC_WT_MUT

et_E12.5_HC_WT_MUT = exactTest(d1, pair=c( "E12_5_HC_WT","E12_5_HC_LHX2_MUT"))
tTags_E12.5_HC_WT_MUT = topTags(et_E12.5_HC_WT_MUT,n=NULL)

is.de_E12.5_HC_WT_MUT = decideTestsDGE(et_E12.5_HC_WT_MUT, p.value = 0.1)
summary(is.de_E12.5_HC_WT_MUT)
plotMD(et_E12.5_HC_WT_MUT)
abline(h=c(-1, 1), col="blue")

sig_DE_E12.5_HC_WT_MUT=tTags_E12.5_HC_WT_MUT[tTags_E12.5_HC_WT_MUT$table$FDR<0.01,]
write.table(sig_DE_E12.5_HC_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/sig_DE_E12.5_HC_WT_MUT', sep='\t', quote=F, row.names=F)

up_sig_DE_E12.5_HC_WT_MUT=sig_DE_E12.5_HC_WT_MUT[sig_DE_E12.5_HC_WT_MUT$table$logFC>0.58,]
write.table(up_sig_DE_E12.5_HC_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/up_sig_DE_E12.5_HC_WT_MUT', sep='\t', quote=F, row.names=F)

down_sig_DE_E12.5_HC_WT_MUT=sig_DE_E12.5_HC_WT_MUT[sig_DE_E12.5_HC_WT_MUT$table$logFC<(-0.58),]
write.table(down_sig_DE_E12.5_HC_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/down_sig_DE_E12.5_HC_WT_MUT', sep='\t', quote=F, row.names=F)
##############################################################################################
###################3  E15.5_CTX_WT_MUT

et_E15.5_CTX_WT_MUT = exactTest(d1, pair=c("E15_5_CTX_WT","E15_5_CTX_LHX2_MUT"))
tTags_E15.5_CTX_WT_MUT = topTags(et_E15.5_CTX_WT_MUT,n=NULL)

is.de_E15.5_CTX_WT_MUT = decideTestsDGE(et_E15.5_CTX_WT_MUT, p.value = 0.1)
summary(is.de_E15.5_CTX_WT_MUT)
plotMD(et_E15.5_CTX_WT_MUT)
abline(h=c(-1, 1), col="blue")


sig_DE_E15.5_CTX_WT_MUT=tTags_E15.5_CTX_WT_MUT[tTags_E15.5_CTX_WT_MUT$table$FDR<0.01,]
write.table(sig_DE_E15.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/sig_DE_E15.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)

up_sig_DE_E15.5_CTX_WT_MUT=sig_DE_E15.5_CTX_WT_MUT[sig_DE_E15.5_CTX_WT_MUT$table$logFC>0.58,]
write.table(up_sig_DE_E15.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/up_sig_DE_E15.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)

down_sig_DE_E15.5_CTX_WT_MUT=sig_DE_E15.5_CTX_WT_MUT[sig_DE_E15.5_CTX_WT_MUT$table$logFC<(-0.58),]
write.table(down_sig_DE_E15.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/down_sig_DE_E15.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)
####################################################################################################
##################  4 E15.5_HC_LHX2_MUT

et_E15.5_HC_WT_MUT = exactTest(d1, pair=c("E15_5_HC_WT","E15_5_HC_LHX2_MUT"))
tTags_E15.5_HC_WT_MUT = topTags(et_E15.5_HC_WT_MUT,n=NULL)

is.de_E15.5_HC_WT_MUT = decideTestsDGE(et_E15.5_HC_WT_MUT, p.value = 0.1)
summary(is.de_E15.5_HC_WT_MUT)
plotMD(et_E15.5_HC_WT_MUT)
abline(h=c(-1, 1), col="blue")


sig_DE_E15.5_HC_WT_MUT=tTags_E15.5_HC_WT_MUT[tTags_E15.5_HC_WT_MUT$table$FDR<0.01,]
write.table(sig_DE_E15.5_HC_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/sig_DE_E15.5_HC_WT_MUT', sep='\t', quote=F, row.names=F)

up_sig_DE_E15.5_HC_WT_MUT=sig_DE_E15.5_HC_WT_MUT[sig_DE_E15.5_HC_WT_MUT$table$logFC>0.58,]
write.table(up_sig_DE_E15.5_HC_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/up_sig_DE_E15.5_HC_WT_MUT', sep='\t', quote=F, row.names=F)

down_sig_DE_E15.5_HC_WT_MUT=sig_DE_E15.5_HC_WT_MUT[sig_DE_E15.5_HC_WT_MUT$table$logFC<(-0.58),]
write.table(down_sig_DE_E15.5_HC_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/down_sig_DE_E15.5_HC_WT_MUT', sep='\t', quote=F, row.names=F)
#####################################################################################################
########### 5 E12.5_E15.5_CTX

et_E12.5_E15.5_CTX = exactTest(d1, pair=c("E12_5_CTX_LHX2_MUT", "E15_5_CTX_LHX2_MUT"))
tTags_E12.5_E15.5_CTX = topTags(et_E12.5_E15.5_CTX,n=NULL)

is.de_E12.5_E15.5_CTX = decideTestsDGE(et_E12.5_E15.5_CTX, p.value = 0.1)
summary(is.de_E12.5_E15.5_CTX)
plotMD(et_E12.5_E15.5_CTX)
abline(h=c(-1, 1), col="blue")


sig_DE_E12.5_E15.5_CTX=tTags_E12.5_E15.5_CTX[tTags_E12.5_E15.5_CTX$table$FDR<0.01,]
write.table(sig_DE_E12.5_E15.5_CTX, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/sig_DE_E12.5_E15.5_CTX', sep='\t', quote=F, row.names=F)

up_sig_DE_E12.5_E15.5_CTX=sig_DE_E12.5_E15.5_CTX[sig_DE_E12.5_E15.5_CTX$table$logFC>0.58,]
write.table(up_sig_DE_E12.5_E15.5_CTX, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/up_sig_E12.5_E15.5_CTX', sep='\t', quote=F, row.names=F)

down_sig_DE_E12.5_E15.5_CTX=sig_DE_E12.5_E15.5_CTX[sig_DE_E12.5_E15.5_CTX$table$logFC<(-0.58),]
write.table(down_sig_DE_E12.5_E15.5_CTX, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/down_sig_DE_E12.5_E15.5_CTX', sep='\t', quote=F, row.names=F)
################################################################################################
#########6 E12.5_E15.5_HC

et_E12.5_E15.5_HC = exactTest(d1, pair=c("E12_5_HC_LHX2_MUT", "E15_5_HC_LHX2_MUT"))
tTags_E12.5_E15.5_HC = topTags(et_E12.5_E15.5_HC,n=NULL)

is.de_E12.5_E15.5_HC = decideTestsDGE(et_E12.5_E15.5_HC, p.value = 0.1)
summary(is.de_E12.5_E15.5_HC)
plotMD(et_E12.5_E15.5_HC)
abline(h=c(-1, 1), col="blue")


sig_DE_E12.5_E15.5_HC=tTags_E12.5_E15.5_HC[tTags_E12.5_E15.5_HC$table$FDR<0.01,]
write.table(sig_DE_E12.5_E15.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/sig_DE_E12.5_E15.5_HC', sep='\t', quote=F, row.names=F)

up_sig_DE_E12.5_E15.5_HC=sig_DE_E12.5_E15.5_HC[sig_DE_E12.5_E15.5_HC$table$logFC>0.58,]
write.table(up_sig_DE_E12.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/up_sig_DE_E12.5_E15.5_HC', sep='\t', quote=F, row.names=F)

down_sig_DE_E12.5_E15.5_HC=sig_DE_E12.5_E15.5_HC[sig_DE_E12.5_E15.5_HC$table$logFC<(-0.58),]
write.table(down_sig_DE_E12.5_E15.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/down_sig_DE_E12.5_E15.5_HC', sep='\t', quote=F, row.names=F)

##################################################################################################
################7  E12.5_CTX_E12.5_HC

et_E12.5_CTX_E12.5_HC = exactTest(d1, pair=c("E12_5_CTX_LHX2_MUT", "E12_5_HC_LHX2_MUT"))
tTags_E12.5_CTX_E12.5_HC = topTags(et_E12.5_CTX_E12.5_HC,n=NULL)

is.de_E12.5_CTX_E12.5_HC = decideTestsDGE(et_E12.5_CTX_E12.5_HC, p.value = 0.1)
summary(is.de_E12.5_CTX_E12.5_HC)
plotMD(et_E12.5_CTX_E12.5_HC)
abline(h=c(-1, 1), col="blue")


sig_DE_E12.5_CTX_E12.5_HC=tTags_E12.5_CTX_E12.5_HC[tTags_E12.5_CTX_E12.5_HC$table$FDR<0.01,]
write.table(sig_DE_E12.5_CTX_E12.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/sig_DE_E12.5_CTX_E12.5_HC', sep='\t', quote=F, row.names=F)

up_sig_DE_E12.5_CTX_E12.5_HC=sig_DE_E12.5_CTX_E12.5_HC[sig_DE_E12.5_CTX_E12.5_HC$table$logFC>0.58,]
write.table(up_sig_DE_E12.5_CTX_E12.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/up_sig_DE_E12.5_CTX_E12.5_HC', sep='\t', quote=F, row.names=F)

down_sig_DE_E12.5_CTX_E12.5_HC=sig_DE_E12.5_CTX_E12.5_HC[sig_DE_E12.5_CTX_E12.5_HC$table$logFC<(-0.58),]
write.table(down_sig_DE_E12.5_CTX_E12.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/down_sig_DE_E12.5_CTX_E12.5_HC', sep='\t', quote=F, row.names=F)

###########################################################################################################
################# 8 E15.5_CTX_E15.5_HC
et_E15.5_CTX_E15.5_HC = exactTest(d1, pair=c("E15_5_CTX_LHX2_MUT", "E15_5_HC_LHX2_MUT"))
tTags_E15.5_CTX_E15.5_HC = topTags(et_E15.5_CTX_E15.5_HC,n=NULL)

is.de_E15.5_CTX_E15.5_HC = decideTestsDGE(et_E15.5_CTX_E15.5_HC, p.value = 0.1)
summary(is.de_E15.5_CTX_E15.5_HC)
a = plotMD(et_E15.5_CTX_E15.5_HC)
abline(h=c(-1, 1), col="blue")


sig_DE_E15.5_CTX_E15.5_HC=tTags_E15.5_CTX_E15.5_HC[tTags_E15.5_CTX_E15.5_HC$table$FDR<0.01,]
write.table(sig_DE_E15.5_CTX_E15.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/sig_DE_E15.5_CTX_E15.5_HC', sep='\t', quote=F, row.names=F)

up_sig_DE_E15.5_CTX_E15.5_HC=sig_DE_E15.5_CTX_E15.5_HC[sig_DE_E15.5_CTX_E15.5_HC$table$logFC>0.58,]
write.table(up_sig_DE_E15.5_CTX_E15.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/up_sig_DE_E15.5_CTX_E15.5_HC', sep='\t', quote=F, row.names=F)

down_sig_DE_E15.5_CTX_E15.5_HC=sig_DE_E15.5_CTX_E15.5_HC[sig_DE_E15.5_CTX_E15.5_HC$table$logFC<(-0.58),]
write.table(down_sig_DE_E15.5_CTX_E15.5_HC, file='/home/darwin/NEW_LHX2_RNA_SEQ/DEG/down_sig_DE_E15.5_CTX_E15.5_HC', sep='\t', quote=F, row.names=F)

############################################################################################################################





#### 1 E12.5_CTX_WT_MUT

et_E12.5_CTX_WT_MUT = exactTest(d1, pair=c("E12_5_CTX_WT","E12_5_CTX_LHX2_MUT"))
tTags_E12.5_CTX_WT_MUT = topTags(et_E12.5_CTX_WT_MUT,n=NULL)

is.de_E12.5_CTX_WT_MUT = decideTestsDGE(et_E12.5_CTX_WT_MUT, p.value = 0.1)
summary(is.de_E12.5_CTX_WT_MUT)
a = plotMD(et_E12.5_CTX_WT_MUT)
abline(h=c(-1, 1), col="blue")


sig_DE_E12.5_CTX_WT_MUT=tTags_E12.5_CTX_WT_MUT[tTags_E12.5_CTX_WT_MUT$table$PValue<0.1,]
write.table(sig_DE_E12.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/test/sig_DE_E12.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)

up_sig_DE_E12.5_CTX_WT_MUT=sig_DE_E12.5_CTX_WT_MUT[sig_DE_E12.5_CTX_WT_MUT$table$logFC>0.58,]
write.table(up_sig_DE_E12.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/test/up_sig_DE_E12.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)

down_sig_DE_E12.5_CTX_WT_MUT=sig_DE_E12.5_CTX_WT_MUT[sig_DE_E12.5_CTX_WT_MUT$table$logFC<(-0.58),]
write.table(down_sig_DE_E12.5_CTX_WT_MUT, file='/home/darwin/NEW_LHX2_RNA_SEQ/test/down_sig_DE_E12.5_CTX_WT_MUT', sep='\t', quote=F, row.names=F)




















