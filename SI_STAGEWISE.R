#lLOADING SAMPLE
library(readr)
A <- read_delim("~/zebra_stage_RNA/lineage/RAW_DATA/SHIELD_BUD_DE_unique_normCPM_PC.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

SI_GENE = A[grep("si:", A$Gene_ID),]

SI_GENES_BUD = SI_GENE[,colnames(SI_GENE)[grep("Gene_ID|_B", colnames(SI_GENE))]]

#MEAN
EC_B_MEAN = rowMeans(SI_GENES_BUD[,c("EC_B.1","EC_B.2")])
NC_B_MEAN = rowMeans(SI_GENES_BUD[,c("NC_B.1","NC_B.2")])
PPL_B_MEAN = rowMeans(SI_GENES_BUD[,c("PPL_B.1","PPL_B.2")])
ME_B_MEAN = rowMeans(SI_GENES_BUD[,c("ME_B.1","ME_B.2")])
EN_B_MEAN = rowMeans(SI_GENES_BUD[,c("EN_B.1","EN_B.2")])

AVG_BUD = cbind(EC_B_MEAN,NC_B_MEAN,PPL_B_MEAN,ME_B_MEAN,EN_B_MEAN)

#STAGE SPECIFIC GENE

SampleB = c("STAGE_EC" , "STAGE_NC" , "STAGE_PPL", "STAGE_ME" , "STAGE_EN" )

out = matrix(ncol = 5,nrow = 1384)

for( i in 1:5 ){
  for(j in 1:5) {
    if (i == j){ next } 
    out[,j] <- AVG_BUD[,i]/AVG_BUD[,j]
  }
  sam_mat   = out[,-i]
  flag = cbind(SI_GENES_BUD$Gene_ID,sam_mat)
  test_ABA = flag[apply(flag[,2:ncol(flag)],1,function(x) all(x>1.5)),]
  
  switch(i,
         
         "i = 1"  =  (colnames(test_ABA) <- c("GENE_ID","EC/NC","EC/PPL","EC/ME","EC/EN")),
         "i = 2"  =  colnames(test_ABA) <- c("GENE_ID","NC/EC","NC/PPL","NC/ME","NC/EN"),
         "i = 3"  =  colnames(test_ABA) <- c("GENE_ID","PPL/EC","PPL/NC","PPL/ME","PPL/EN"),
         "i = 4"  =  colnames(test_ABA) <- c("GENE_ID","ME/EC","ME/NC","ME/PPL","ME/EN"),
         "i = 5"  =  (colnames(test_ABA) <- c("GENE_ID","EN/EC","EN/NC","EN/PPL","EN/ME")),
  )
  
  assign(SampleB[[i]],test_ABA)
}

#########################################################################################
# WRITING TXT FILE AND HEAT MAP

SampleB = c("STAGE_EC" , "STAGE_NC" , "STAGE_PPL", "STAGE_ME" , "STAGE_EN" ) 
MATRICES = list(STAGE_EC,STAGE_NC,STAGE_PPL,STAGE_ME,STAGE_EN)

SI_BUD_COMMON_MATRIC_NAME = c("STAGE_BUD_EC" , "STAGE_BUD_NC" , "STAGE_BUD_PPL", "STAGE_BUD_ME" , "STAGE_BUD_EN")
MATRICES2 = list(STAGE_BUD_EC,STAGE_BUD_NC,STAGE_BUD_PPL,STAGE_BUD_ME,STAGE_BUD_EN)
for(i in 1:5){
  common_SI_genes = merge(SI_GENES_BUD, MATRICES[i], by.x = "Gene_ID", by.y = "GENE_ID") 
  COMMOM_STAGEWISE = common_SI_genes[,-c(12,13,14,15)] 
  
  file_location =  paste("/home/darwin/zebra_stage_RNA/lineage/SI_BUD", SampleB[[i]], sep = "/") 
  txt_file_name = paste(file_location, ".txt", sep = "")
  
  
  print(file_location)
  print(txt_file_name)
  write.table(COMMOM_STAGEWISE, file =  txt_file_name ,row.names = FALSE, sep = "\t")
  
  
  tabel_read = read.table(txt_file_name,"\t", header=T,row.names=1)
  
  assign(SI_BUD_COMMON_MATRIC_NAME[[i]], tabel_read)
  
  converted_matrix = data.matrix(tabel_read,1:ncol(tabel_read))
  #converted_matrix_log = log(converted_matrix,1:ncol(converted_matrix))
  
  my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
  heatmap.2(converted_matrix, Colv = NA, Rowv = NA , trace = "none",col = my_palette,
            key = TRUE, key.xlab = "Log(COUNTS)", key.ylab = NULL , main = SampleB[[i]],
            xlab = "STAGES", ylab = "Log(COUNTS)")
}




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


COMPARISON = c("et_E12.5_CTX_WT_MUT","et_E12.5_HC_WT_MUT","et_E15.5_HC_WT_MUT","et_E15.5_CTX_WT_MUT","et_E12.5_E15.5_CTX","et_E12.5_E15.5_HC","et_E12.5_CTX_E12.5_HC","et_E15.5_CTX_E15.5_HC")

COMP = c("E12.5_CTX_WT_MUT","E12.5_HC_WT_MUT","E15.5_HC_WT_MUT","E15.5_CTX_WT_MUT","E12.5_E15.5_CTX","E12.5_E15.5_HC","E12.5_CTX_E12.5_HC","E15.5_CTX_E15.5_HC")

et_tags = matrix(ncol = 1,nrow = 8)
 for(i in 1:5){
	print(tags_index = topTags(COMP[[i]],n=NULL))
	
	print(file_location =  paste("/home/darwin/NEW_LHX2_RNA_SEQ/DEG", COMP[[i]], sep = "/") )
  	print(txt_file_name = paste(file_location, ".txt", sep = ""))
  	print(up_deg = tags_index[[i]][tags_index$table$FDR<0.01,])
	print(write.table(up_deg, file =  txt_file_name ,row.names = FALSE, sep = "\t"))

}



















