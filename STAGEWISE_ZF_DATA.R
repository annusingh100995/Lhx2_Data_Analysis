p <- vector("numeric", length = nrow(AVG_BUD))
p = NULL
for(i in colnames(AVG_BUD)){
  for(j in 1:5) {
    if (i == j){ next } 
  p[[i]] = AVG_BUD[[i]]/AVG_BUD[[j]]
  }
  
  #print(mydata[[j]])
  #test[i] = AVG_BUD[,i]/AVG_BUD[,j]
  #assign(paste("x",j,sep=""), AVG_BUD[,i]/AVG_BUD[,j])
  #assign(paste("BC",j,sep="_"), cbind(y[j]))
  #cbind(xj)
  #assign(paste("x",i,sep="_"), cbind(SI_GENES_BUD$Gene_ID,AVG_BUD[,i]/AVG_BUD[,j]))
  #cbind(x[i])  

  mydata = list()
# mydf = data.frame()
  for(i in 2){
    
    for(j in 1:5) {
      
      if (i == j){ next } 
     
      mydata[[j]] = AVG_BUD[,i]/AVG_BUD[,j]
      #mydf[[i]] = data.frame(mydata = unlist(mydata[[j]]))
    }
    mydf <- data.frame(mydata)
  }
  #assign(paste("t",i,sep="_"),cbind.data.frame(assign(paste("x",j,sep="_"), (AVG_BUD[,i]/AVG_BUD[,j]))))
  
  
  
  
  
  
  ###############################################################################################################
  Sample = c("Z" , "B" , "C", "D" , "E" )
   
 out = matrix(ncol = 5,nrow = 1384)
  for(k in Sample){
      for( i in 1:2 ){
          for(j in 1:5) {
               if (i == j){ next } 
               out[,j] <- AVG_BUD[,i]/AVG_BUD[,j]
                        }
     sam_mat   = out[,-i]
                                                                       #ABA = out[,-i]
                   }
  test_ABA = sam_mat
      assign(Sample[[i]],test_ABA)
  }
 
 colnames(EC2) = c("GENE_ID","NC","PPL","ME","EN")
 
##########################################################################################################
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
   
 
###################################################################################################################

 common44_ppl = merge(SI_GENES_BUD, STAGE_PPL, by.x = "Gene_ID", by.y = "GENE_ID") 
 COMMOM_BUD_PPL = common_ppl[,-c(12,13,14,15)]

 write.table(COMMOM_BUD_PPL, file = "/home/darwin/zebra_stage_RNA/lineage/SI_BUD/COMMON_SI_BUD_PPL.txt",row.names = FALSE, sep = "\t")
 BUD_PPL = read.table("/home/darwin/zebra_stage_RNA/lineage/SI_BUD/COMMON_SI_BUD_PPL.txt","\t", header=T,row.names=1)
 DATA_BUD_PPL = data.matrix(BUD_PPL[,1:ncol(BUD_PPL)])
 DATA_BUD_PPL_LOG = log(DATA_BUD_PPL[,1:ncol(DATA_BUD_PPL)])
 my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
 heatmap.2(DATA_BUD_PPL_LOG, Colv = NA, Rowv = NA , trace = "none",col = my_palette,
           key = TRUE, key.xlab = "Log(COUNTS)", key.ylab = NULL , main = "PPL",
           xlab = "STAGES", ylab = "Log(COUNTS)")
 
 
 SampleB = c("STAGE_EC" , "STAGE_NC" , "STAGE_PPL", "STAGE_ME" , "STAGE_EN" ) 
 
 FILENAME = c("COMMON_SI_BUD_EC.txt","COMMON_SI_BUD_NC.txt","COMMON_SI_BUD_PPL.txt","COMMON_SI_BUD_ME.txt","COMMON_SI_BUD_EN.txt")
 output_column = matrix(ncol = 15,nrow = 1384)
 remove_FC = matrix(ncol = 11,nrow = 1384)
 load_matrix =  matrix(ncol = 15,nrow = 1384)
for(i in 1:5){
  print( as.name(SampleB[[i]]))
  foo =  as.name(SampleB[[i]])
  load_matrix = as.matrix(foo)
 out_common = merge(SI_GENES_BUD ,load_matrix, by.x = "Gene_ID", by.y = "GENE_ID")
  #remove_FC = as.matrix(out_common[,-c(12,13,14,15)])
  #write.table(remove_FC,file = "/home/darwin/zebra_stage_RNA/lineage/SI_BUD/SampleB[[i]]",row.names = FALSE, sep = "\t")
}
 
 
 
 
 
 
 
 my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
 heatmap.2(COMMOM_BUD_PPL[,2:ncol(COMMOM_BUD_PPL)], Colv = NA, Rowv = NA , trace = "none",col = my_palette,
           key = TRUE,
           key.xlab = "COUNTS", key.ylab = NULL , main = "PPL",
           xlab = "STAGES", ylab = "COUNTS")
 
 
 
 
 ANO = read.table("~/zebra_stage_RNA/RNASeq/ANO.txt","\t", header=T,row.names=1)
 DATA_ANO = data.matrix(ANO[,2:ncol(ANO)])
 my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
 heatmap.2(DATA_ANO, Colv = NA, Rowv = NA , trace = "none",col = my_palette,
           key = TRUE,
           key.xlab = "Log(CPM)", key.ylab = NULL , main = "Anoctamin",
           xlab = "STAGES", ylab = "Log(CPM)")
 
  
 
 DATA_COMMOM_PPL = data.matrix(common_ppl[,2:11]) 
 my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
 heatmap.2(DATA_COMMOM_PPL, Colv = NA, Rowv = NA , trace = "none",col = my_palette,
           key = TRUE,
           key.xlab = "COUNTS", key.ylab = NULL , main = "PPL",
           xlab = "STAGES", ylab = "COUNTS")
 
 
############################################################################################################################################ 
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
   converted_matrix_log = log(tabel_read,1:ncol(tabel_read))
   
   my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
   heatmap.2(converted_matrix_log, Colv = NA, Rowv = NA , trace = "none",col = my_palette,
             key = TRUE, key.xlab = "Log(COUNTS)", key.ylab = NULL , main = SampleB[[i]],
             xlab = "STAGES", ylab = "Log(COUNTS)")
 }
 
 
 common44_ppl = merge(SI_GENES_BUD, STAGE_PPL, by.x = "Gene_ID", by.y = "GENE_ID") 
 COMMOM_BUD_PPL = common_ppl[,-c(12,13,14,15)]
  write.table(COMMOM_BUD_PPL, file = "/home/darwin/zebra_stage_RNA/lineage/SI_BUD/COMMON_SI_BUD_PPL.txt",row.names = FALSE, sep = "\t")
 BUD_PPL = read.table("/home/darwin/zebra_stage_RNA/lineage/SI_BUD/COMMON_SI_BUD_PPL.txt","\t", header=T,row.names=1)
 DATA_BUD_PPL = data.matrix(BUD_PPL[,1:ncol(BUD_PPL)])
 DATA_BUD_PPL_LOG = log(DATA_BUD_PPL[,1:ncol(DATA_BUD_PPL)])
 my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
 heatmap.2(DATA_BUD_PPL_LOG, Colv = NA, Rowv = NA , trace = "none",col = my_palette,
           key = TRUE, key.xlab = "Log(COUNTS)", key.ylab = NULL , main = "PPL",
           xlab = "STAGES", ylab = "Log(COUNTS)")
 
 ########################################################################################
 go_id_chr_remod = GOID( GOTERM[ Term(GOTERM) == "chromatin remodeling"])
 allegs = get(go_id_chr_remod,org.Dr.egGO2ALLEGS)
 genes_chromatin_remodeling = unlist(mget(allegs,org.Dr.egSYMBOL))
 
 GENES_BUD = A[!(A$Gene_ID %in% SI_GENE$Gene_ID),]
 
Matched_gene2 = GENES_BUD[GENES_BUD$Gene_ID %in% genes_chromatin_remodeling ,] 


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
