source("http://bioconductor.org/biocLite.R")
library("org.Hs.eg.db")
library(clusterProfiler)
library(mygene)
library(topGO)
library("pathview")
library("Rgraphviz")
setwd("D:\\工作\\大项目\\毓璜顶医院\\部分结果\\mir表达\\mir靶基因\\GCM_DCM")


Mutation_genes<-c('GUCA1C','NT5C1A','PDK4','SMCO1')
Mutation_genes<-read.table("GCM_DCM_靶基因.txt")

Mutation_genes_genes<-queryMany(Mutation_genes$V1,scopes="symbol", fields=c("uniprot", "entrezgene","ensembl"), species="human",returnall=T)
Mutation_genes_ggo_cc <- groupGO(gene = as.character(Mutation_genes_genes$response$entrezgene),OrgDb    = org.Hs.eg.db, ont = "CC",level    = 3,readable = TRUE,keytype="ENTREZID")
Mutation_genes_ggo_bp <- groupGO(gene = as.character(Mutation_genes_genes$response$entrezgene),OrgDb    = org.Hs.eg.db, ont = "BP",level    = 3,readable = TRUE,keytype="ENTREZID")
Mutation_genes_ggo_mf <- groupGO(gene = as.character(Mutation_genes_genes$response$entrezgene),OrgDb    = org.Hs.eg.db, ont = "MF",level    = 3,readable = TRUE,keytype="ENTREZID")
Mutation_genes_ggo_mf_sort<-Mutation_genes_ggo_mf[order(Mutation_genes_ggo_mf[,3],decreasing=T),]
Mutation_genes_ggo_cc_sort<-Mutation_genes_ggo_cc[order(Mutation_genes_ggo_cc[,3],decreasing=T),]
Mutation_genes_ggo_bp_sort<-Mutation_genes_ggo_bp[order(Mutation_genes_ggo_bp[,3],decreasing=T),]
####barplot
pdf("GCMvsDCM_targets_Group_GO_MF.barplot.pdf")
barplot(Mutation_genes_ggo_mf,drop=TRUE, showCategory=25,font.size = 10,width=2,title="Molecular functions of target genes of DE miRNAs of DCM vs Normal")
dev.off()
png("GCMvsDCM_targets_Group_GO_CC.barplot.png")
barplot(Mutation_genes_ggo_cc,drop=TRUE, showCategory=25,font.size = 10,title="Cellular components of target genes of DE miRNAs of DCM vs Normal")
dev.off()
png("GCMvsDCM_targets_Group_GO_BP.barplot.png")
barplot(Mutation_genes_ggo_bp,drop=TRUE, showCategory=25,font.size = 10,title="Biological process of target genes of DE miRNAs of DCM vs Normal")
dev.off()


write.table(Mutation_genes_ggo_mf_sort,"GCMvsDCM_targets_Group_GO_MF.txt",sep="\t")
write.table(Mutation_genes_ggo_cc_sort,"GCMvsDCM_targets_Group_GO_CC.txt",sep="\t")
write.table(Mutation_genes_ggo_bp_sort,"GCMvsDCM_targets_Group_GO_BP.txt",sep="\t")


Mutation_genes_ego_CC <- enrichGO(gene=as.character(Mutation_genes_genes$response$entrezgene),OrgDb=org.Hs.eg.db,keytype="ENTREZID",ont="CC",pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,minGSSize=1,readable=T)
Mutation_genes_ego_BP <- enrichGO(gene=as.character(Mutation_genes_genes$response$entrezgene),OrgDb=org.Hs.eg.db,keytype="ENTREZID",ont="BP",pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,minGSSize=1,readable=T)
Mutation_genes_ego_MF <- enrichGO(gene=as.character(Mutation_genes_genes$response$entrezgene),OrgDb=org.Hs.eg.db,keytype="ENTREZID",ont="MF",pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,minGSSize=1,readable=T)

png("GCMvsDCM_targets_Enrich_GO_MF.png")
barplot(Mutation_genes_ego_MF,drop=TRUE, showCategory=25,font.size = 10)
dev.off()
png("GCMvsDCM_targets_Enrich_GO_CC.png")
barplot(Mutation_genes_ego_CC,drop=TRUE, showCategory=25,font.size = 10)
dev.off()
png("GCMvsDCM_targets_Enrich_GO_BP.png")
barplot(Mutation_genes_ego_BP,drop=TRUE, showCategory=25,font.size = 10)
dev.off()

####dotplot
png("GCMvsDCM_targets_Group_GO_MF.dotplot.png")
dotplot(Mutation_genes_ego_MF, showCategory=25,font.size = 10)
dev.off()
png("GCMvsDCM_targets_Group_GO_CC.dotplot.png")
dotplot(Mutation_genes_ego_CC, showCategory=15,font.size = 10)
dev.off()
png("GCMvsDCM_targets_Group_GO_BP.dotplot.png")
dotplot(Mutation_genes_ego_BP, showCategory=25,font.size = 10)
dev.off()

write.table(Mutation_genes_ego_CC,"GCMvsDCM_targets_Enrich_GO_CC.txt",sep="\t")
write.table(Mutation_genes_ego_MF,"GCMvsDCM_targets_Enrich_GO_MF.txt",sep="\t")
write.table(Mutation_genes_ego_BP,"GCMvsDCM_targets_Enrich_GO_BP.txt",sep="\t")

####plotGOgraph
pdf("GCMvsDCM_targets_Enrich_plotGOgraph_GO_BP.dotplot.pdf")
plotGOgraph(Mutation_genes_ego_BP,firstSigNodes = 10)
dev.off()
pdf("GCMvsDCM_targets_Enrich_plotGOgraph_GO_CC.dotplot.pdf")
plotGOgraph(Mutation_genes_ego_CC,firstSigNodes = 20)
dev.off()
pdf("GCMvsDCM_targets_Enrich_plotGOgraph_GO_MF.dotplot.pdf")
plotGOgraph(Mutation_genes_ego_MF,firstSigNodes = 20)
dev.off()

#####################
#####kegg enrichment
#####################
Mutation_genes_kk_Mutation_genes_KEGG<-enrichKEGG(gene=Mutation_genes_genes$response$entrezgene,organism='hsa',pvalueCutoff=1,qvalueCutoff=1,keyType="ncbi-geneid",pAdjustMethod="BH",minGSSize=1)

write.table(Mutation_genes_kk_Mutation_genes_KEGG,"GCMvsDCM_targets_Enrich_KEGG.txt",sep="\t")



