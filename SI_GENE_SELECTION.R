AVG_BUD = cbind(EC_B,NC_B,PPL_B,ME_B,EN_B)
max = apply(AVG_BUD[,1:ncol(AVG_BUD)],1,max)
AVG_BUD = cbind(AVG_BUD,max)
AVG_BUD = as.data.frame(AVG_BUD)
test = sweep(AVG_BUD[,1:5],1,AVG_BUD[,6],"/")
test2 = apply(test,2,function(x) x^-1)

col_max = colnames(AVG_BUD[,apply(AVG_BUD,1,which.max)])
test3 = cbind(SI_GENES_BUD$Gene_ID,test2,col_max)
test3 = as.data.frame(test3)
test4 = test3[grep("PPL",test3$col_max),]

test5 = test4[apply(test4[,-c(4,7)],1, function(x) all(x>1.5)),]
test6 = test5[,-7]
 
test7 = intersect(SI_GENES_BUD$Gene_ID , test6$V1)
test8 = merge(SI_GENES_BUD, test6 , by.x = "Gene_ID", by.y = "V1")

AVG_BUD =  cbind(SI_GENES_BUD$Gene_ID, AVG_BUD)
 AVG_BUD =  as.numeric(AVG_BUD[,2:6])
 
 
 
PPL = cbind(SI_GENES_BUD$Gene_ID , AVG_BUD[,3]/AVG_BUD[,1],AVG_BUD[,3]/AVG_BUD[,2],AVG_BUD[,3]/AVG_BUD[,4],AVG_BUD[,3]/AVG_BUD[,5])
PPL2 = PPL[apply(PPL[,2:ncol(PPL)],1,function(x) all(x>1.5)),]
colnames(PPL2) = c("GENE_ID","PPL/EC","PPL/NC","PPL/ME","PPL/EN")

EC = cbind(SI_GENES_BUD$Gene_ID,AVG_BUD[,1]/AVG_BUD[,2],AVG_BUD[,1]/AVG_BUD[,3],AVG_BUD[,1]/AVG_BUD[,4],AVG_BUD[,1]/AVG_BUD[,5])
EC2 = EC[apply(EC[,2:ncol(EC)],1,function(x) all(x>1.5)),]
colnames(EC2) = c("GENE_ID","NC","PPL","ME","EN")




Sample = c("EC_B","NC_B","PPL_B","ME_B","EN_B")

for (i in Sample){
  assign(i, AVG_BUD[[i]])
}

y = numeric(length(EC_B))
for (i in c(1,2)){
 for(k in seq_along(EC_B)){
  y[k] = Sample[[i]]/Sample[[(i+1)]]
   paste(i,'RATIO',"_") = y[k]
   
    }
}
rm(x_2,x_3,x_4,x_5)
x <- vector("list",27) 
for (i in 1:27) {
  x[[i]] <- cbind(rep(1,aux[i]), rnorm(aux[i]))
}

x = vector("y",5)
v <- c()
for(i in 1){
  
  for(j in 1:5) {
    if (i == j){ next } 
    print(j)
    
    #test[i] = AVG_BUD[,i]/AVG_BUD[,j]
    
     assign(paste("x",j,sep=""), AVG_BUD[,i]/AVG_BUD[,j])
     #assign(paste("BC",j,sep="_"), cbind(y[j]))
    cbind(xj)
    #assign(paste("x",i,sep="_"), cbind(SI_GENES_BUD$Gene_ID,AVG_BUD[,i]/AVG_BUD[,j]))
    #cbind(x[i])
  }
}


result <- vector("list",5)
for(i in 5){
  for(j in 1:5) {
    if (i == j){ next } 
   assign(paste("x",j,sep=""), AVG_BUD[,i]/AVG_BUD[,j])
  }
}
rm(x1,x2,x3,x4,x5)

