
E12.5_CTX_WT_MUT

up_dge_list = up_sig_DE_E12.5_CTX_WT_MUT$table$Geneid
down_dge_list = down_sig_DE_E12.5_CTX_WT_MUT$table$Geneid
DGE_E12.5_CTX_WT_MUT = c(up_dge_list,down_dge_list)
length(DGE_E12.5_CTX_WT_MUT)
eg_E12.5_CTX_WT_MUT <- bitr(DGE_E12.5_CTX_WT_MUT, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)  

eg_E12.5_CTX_WT_MUT
gene_E12.5_CTX_WT_MUT = as.vector(eg_E12.5_CTX_WT_MUT$ENTREZID)
ggo_E12.5_CTX_WT_MUT<- groupGO(gene     = gene,
                                 OrgDb    = org.Mm.eg.db,
                                 level    = 3,
                                 ont      = "CC",
                                 readable = TRUE)

ego_E12.5_CTX_WT_MUT <- enrichGO(gene          = gene_E12.5_CTX_WT_MUT,
                                   OrgDb         = org.Mm.eg.db,
                                   ont           = "CC",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
head(ego_E12.5_CTX_WT_MUT)

dotplot(ego_E12.5_CTX_WT_MUT)
enrichMap(ego_E12.5_CTX_WT_MUT)
plotGOgraph(ego_E12.5_CTX_WT_MUT)

E12.5_CTX_WT_MUT_gene_ratio = as.vector(ego_E12.5_CTX_WT_MUT@result$GeneRatio)
E12.5_CTX_WT_MUT_description = as.vector(ego_E12.5_CTX_WT_MUT@result$Description)
E12.5_CTX_WT_MUT_go_ids = as.vector(ego_E12.5_CTX_WT_MUT@result$ID)
E12.5_CTX_WT_MUT_gene_ids = as.vector(ego_E12.5_CTX_WT_MUT@result$geneID)
E12.5_CTX_WT_MUT = cbind(E12.5_CTX_WT_MUT_gene_ratio,E12.5_CTX_WT_MUT_description,E12.5_CTX_WT_MUT_go_ids,E12.5_CTX_WT_MUT_gene_ids)


E12.5_E12.5_HC_WT_MUT

up_dge_list = up_sig_DE_E12.5_HC_WT_MUT$table$Geneid
down_dge_list = down_sig_DE_E12.5_HC_WT_MUT$table$Geneid
DGE_E12.5_HC_WT_MUT = c(up_dge_list,down_dge_list)
length(DGE_E12.5_HC_WT_MUT)
eg_E12.5_HC_WT_MUT <- bitr(DGE_E12.5_HC_WT_MUT, fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)  

eg_E12.5_HC_WT_MUT
gene_E12.5_HC_WT_MUT = as.vector(eg_E12.5_HC_WT_MUT$ENTREZID)
ggo_E12.5_HC_WT_MUT<- groupGO(gene     = gene_E12.5_HC_WT_MUT,
                               OrgDb    = org.Mm.eg.db,
                               level    = 3,
                               ont      = "CC",
                               readable = TRUE)

ego_E12.5_HC_WT_MUT<- enrichGO(gene          = gene_E12.5_HC_WT_MUT,
                                 OrgDb         = org.Mm.eg.db,
                                 ont           = "CC",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE)
head(ego_E12.5_HC_WT_MUT)

dotplot(ego_E12.5_HC_WT_MUT)
enrichMap(ego_E12.5_HC_WT_MUT)
plotGOgraph(ego_E12.5_HC_WT_MUT)

E12.5_HC_WT_MUT_gene_ratio = as.vector(ego_E12.5_HC_WT_MUT@result$GeneRatio)
E12.5_HC_WT_MUT_description = as.vector(ego_E12.5_HC_WT_MUT@result$Description)
E12.5_HC_WT_MUT_go_ids = as.vector(ego_E12.5_HC_WT_MUT@result$ID)
E12.5_HC_WT_MUT_gene_ids = as.vector(ego_E12.5_HC_WT_MUT@result$geneID)
E12.5_HC_WT_MUT = cbind(E12.5_HC_WT_MUT_gene_ratio,E12.5_HC_WT_MUT_description,E12.5_HC_WT_MUT_go_ids,E12.5_HC_WT_MUT_gene_ids)
