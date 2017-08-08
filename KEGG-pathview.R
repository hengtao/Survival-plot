#!/usr/local/bin/R
library(pathview)
##�л�Ŀ¼
setwd("D:\\����\\����Ŀ\\ع諶�ҽԺ\\���ֽ��\\���Ի���-����\\GCM_DCM\\KEGG\\")
setwd("D:/����/����Ŀ/ع諶�ҽԺ/���ֽ��/���Ի���/DCM_Normal/pathview_kegg/DCM_Normal_Cardi_kegg")
##���غ���Ѫ����ص�keggͨ·
for (i in c("00230","00240","00260","00270","00350","00450","00500","00532","00534","00561","00564","00600","00603","00670","00740","00760","00770","01100","02010","03320","04010","04020","04060","04066","04080","04142","04144","04145","04150","04260","04270","04330","04350","04370","04380","04390","04510","04512","04520","04530","04540","04810","04910","04913","04920","04960","04975","04976","05410","05412","05414","05416")){
  download.kegg(pathway.id = i, species = "hsa", kegg.dir = ".",file.type=c("xml", "png")) 
}
GCM_Normal_kegg<-read.table("D:\\����\\����Ŀ\\ع諶�ҽԺ\\���ֽ��\\���Ի���-����\\GCM_DCM\\KEGG\\����ת¼��GCM_DCM.txt",header=T,row.names=1,sep="\t")
DCM_Normal_kegg<-read.table("D:/����/����Ŀ/ع諶�ҽԺ/���ֽ��/���Ի���/DCM_Normal/pathview_kegg/DCMvsNormal_genes_for_kegg_rld.txt",header=T,row.names=1,sep="\t")
Mutation_WGS_kegg<-read.table("D:/����/����Ŀ/ع諶�ҽԺ/���ֽ��/��Ҫ����/GO_KEGG/Variants_genes.txt",header=T,row.names=1,sep="\t")
class(Mutation_WGS_kegg) ##data frame��ʽ
Mutation_WGS_kegg_log<-log2(Mutation_WGS_kegg+1)
Mutation_WGS_kegg$GCM<-rowMeans(Mutation_WGS_kegg[,1:3])
Mutation_WGS_kegg$Normal<-rowMeans(Mutation_WGS_kegg[,4:8])
Mutation_WGS_kegg$DCM<-rowMeans(Mutation_WGS_kegg[,9:17])
Mutation_WGS_kegg_log$GCM_Normal<-Mutation_WGS_kegg_log[,2]-Mutation_WGS_kegg_log[,1]
##��Mutation_WGS_keggתΪmatrix,Ϊ���ܹ����к����ķ���
Mutation_WGS_kegg_log_matrix<-as.matrix(Mutation_WGS_kegg_log)
Mutation_WGS_kegg_log_matrix[,3]
##�Ժ���Ѫ����ص�����ͨ·���з���
for (i in c("00230","00240","00260","00270","00350","00450","00500","00532","00534","00561","00564","00600","00603","00670","00740","00760","00770","01100","02010","03320","04010","04020","04060","04066","04080","04142","04144","04145","04150","04260","04270","04330","04350","04370","04380","04390","04510","04512","04520","04530","04540","04810","04910","04913","04920","04960","04975","04976","05410","05412","05414","05416")){
  pv.out<-pathview(gene.data=Mutation_WGS_kegg_log_matrix[,3],pathway.id=i,species="hsa",out.suffix="GCM_DCM_kegg",gene.idtype="SYMBOL",kegg.native=T,same.layer=F,kegg.dir=".")
  
}

####����DCM��normal���Ի����kegg
setwd("D:/����/����Ŀ/ع諶�ҽԺ/���ֽ��/���Ի���/DCM_Normal/pathview_kegg/DCM_Normal_Cardi_kegg")
setwd("D:/����/����Ŀ/ع諶�ҽԺ/���ֽ��/��Ҫ����/GO_KEGG")