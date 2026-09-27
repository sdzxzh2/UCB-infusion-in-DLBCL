
library(R.utils)
library(GenomicFeatures)
library(ChIPseeker)
library(ggplot2)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

con_sam="T1"
test_sam=c("T2","T3","T4")
con_sam_all_position_dir=paste0("./",con_sam,"/results/Step4_VariantCalling/",con_sam,".calling.step1.tsv")
con_sam_mut_position_dir=paste0("./",con_sam,"/results/Step4_VariantCalling/",con_sam,".calling.step2.tsv")
cnv <- fread(con_sam_all_position_dir)
cnv <- as.data.frame(cnv)
con_sam_all_position=cnv
rm(cnv)
gc()
cnv <- fread(con_sam_mut_position_dir)
cnv <- as.data.frame(cnv)
con_sam_mut_position=cnv
rm(cnv)
gc()
con_sam_mut_position[,3]=paste0(con_sam_mut_position[,1],"_",con_sam_mut_position[,2])
list_mut_sam_position=list()

for (i in test_sam) {
  mut_sam_mut_position_dir=paste0("./",i,"/results/Step4_VariantCalling/",i,".calling.step2.tsv")
  cnv <- fread(mut_sam_mut_position_dir)
  cnv <- as.data.frame(cnv)
  mut_sam_mut_position=cnv
  list_mut_sam_position[[i]]=mut_sam_mut_position
  
}

rm(cnv)
gc()

for (i in test_sam ) {
  tmp=list_mut_sam_position[[i]]
  tmp[,3]=paste0(tmp[,1],"_",tmp[,2])
  list_mut_sam_position[[i]]=tmp
}

position=c()
for (i in test_sam ) {
  tmp=list_mut_sam_position[[i]]
  tmp2=as.vector(tmp[3])[c(1:nrow(tmp)),1]
  position=c(position,tmp2)
}

position=table(position)
position=as.data.frame(t(position))
position_mut_fin=position[position$Freq==length(test_sam),]$position
position_con_fin=as.vector(con_sam_mut_position[3])[c(1:nrow(con_sam_mut_position)),1]
consam=position_mut_fin%in%position_con_fin
position_fin=position_mut_fin[which(consam==FALSE)]
position_fin=as.vector(position_fin)
fin=data.frame(position_fin)
for (i in test_sam ) {
  tmp=list_mut_sam_position[[i]]
  tmp=tmp[tmp[,3]%in%position_fin,]
  fin$position_fin=tmp[,3]
  fin[,paste0("REF","_",i)]=tmp$REF
  fin[,paste0("ALT","_",i)]=tmp$ALT
  fin[,paste0("Cell_types","_",i)]=tmp$Cell_types
  fin[,paste0("FILTER","_",i)]=tmp$FILTER
  fin[,paste0("Cell_type_Filter","_",i)]=tmp$Cell_type_Filter
  fin[,paste0("Cell_types_min_BC","_",i)]=tmp$Cell_types_min_BC
  fin[,paste0("Cell_types_min_CC","_",i)]=tmp$Cell_types_min_CC
  fin[,paste0("Nc","_",i)]=tmp$Nc
  fin[,paste0("Cc","_",i)]=tmp$Cc
  fin[,paste0("DC","_",i)]=tmp$DC
  fin[,paste0("Mono","_",i)]=tmp$Mono
  fin[,paste0("NK","_",i)]=tmp$NK
  fin[,paste0("platelet","_",i)]=tmp$platelet
  fin[,paste0("T","_",i)]=tmp$T
  fin[,paste0("Tprolif","_",i)]=tmp$Tprolif
  fin[,paste0("VAF","_",i)]=tmp$VAF
  fin[,paste0("CCF","_",i)]=tmp$CCF
  fin[,paste0("BCp","_",i)]=tmp$BCp
  fin[,paste0("CCp","_",i)]=tmp$CCp
  fin[,paste0("Rest_BC","_",i)]=tmp$Rest_BC
  fin[,paste0("Rest_CC","_",i)]=tmp$Rest_CC
}

tmp=strsplit(position_fin,"_")
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
tmp=as.vector(unique(tmp$V2))
con_sam_all_position_sub=con_sam_all_position[con_sam_all_position[,2]%in%tmp,]
con_sam_all_position_sub[,3]=paste0(con_sam_all_position_sub[,1],"_",con_sam_all_position_sub[,2])
con_sam_all_position_fin=con_sam_all_position_sub[con_sam_all_position_sub[,3]%in%position_fin,]
rownames(con_sam_all_position_fin)=con_sam_all_position_fin[,3]
rownames(fin)=fin$position_fin


library("doParallel")     
library("foreach")       
cl<-makeCluster(14)
registerDoParallel(cl)

mydata <- foreach(i=rownames(fin),.combine = rbind) %dopar% {
  tmp=con_sam_all_position_fin[rownames(con_sam_all_position_fin)%in%i,]
  if (nrow(tmp)==0) {
    tmp[1,]=NA
  }
  
  sam="CON"
  dat=c(tmp$REF,tmp$ALT,tmp$Cell_types,tmp$FILTER,tmp$Cell_type_Filter,tmp$Cell_types_min_BC,
        tmp$Cell_types_min_CC,tmp$Nc,tmp$Cc,tmp$DC,tmp$Mono,tmp$NK,tmp$platelet,tmp$T,
        tmp$Tprolif,tmp$VAF,tmp$CCF,tmp$BCp,tmp$CCp,tmp$Rest_BC,tmp$Rest_CC)
  
  return(dat)
  
  
}

stopCluster(cl)

sam="CON"
mydata=as.data.frame(mydata)
rownames(mydata)=rownames(fin)
colnames(mydata)=c(paste0("REF","_",sam),paste0("ALT","_",sam),paste0("Cell_types","_",sam),
                   paste0("FILTER","_",sam),paste0("Cell_type_Filter","_",sam),paste0("Cell_types_min_BC","_",sam),
                   paste0("Cell_types_min_CC","_",sam),paste0("Nc","_",sam),paste0("Cc","_",sam),
                   paste0("DC","_",sam),
                   paste0("Mono","_",sam),paste0("NK","_",sam),paste0("platelet","_",sam),paste0("T","_",sam),
                   paste0("Tprolif","_",sam),paste0("VAF","_",sam),paste0("CCF","_",sam),paste0("BCp","_",sam),
                   paste0("CCp","_",sam),paste0("Rest_BC","_",sam),paste0("Rest_CC","_",sam)
)

fin=cbind(fin,mydata)
dat=data.frame()
for (sam in test_sam) {
  cnv <- fread(paste0("cell_genotype_mut_read_3_",sam,".csv"))
  cnv <- as.data.frame(cnv)
  tmp=cnv
  dat=rbind(dat,tmp)
}
dat= dat[,-1]
write.csv(dat,"cell_genotype_mut_read_3_all_sam.csv")
pos=data.frame(dat[,1],as.numeric(dat[,2]))
colnames(pos)=c("#CHROM","POS")
pos=as.data.frame(pos)
library(dplyr)
pos=distinct(pos)
write.table(pos, "vcf.txt", sep = "\t", row.names = FALSE)

ann <- read.delim("./ann.exonic_variant_function", header=FALSE, row.names=1)
ann$chr_pos_mut=paste0(ann$V4,"_",ann$V5,"_",ann$V8)
rownames(ann)=ann$chr_pos_mut
mut_info=strsplit(ann$V3,":")
gene=lapply(mut_info,FUN = function(x){x[1]})
gene=unlist(gene)
protein_mut=lapply(mut_info,FUN = function(x){x[5]})
protein_mut=unlist(protein_mut)
protein_mut=strsplit(protein_mut,",")
protein_mut=lapply(protein_mut,FUN = function(x){x[1]})
protein_mut=unlist(protein_mut)
ann$SYMBOL=gene
ann$Protein_Change=protein_mut
write.csv(ann,"cell_genotype_mut_read_3_all_sam_ann_position_annovar.csv")
pos_ann=read.csv("./cell_genotype_mut_read_3_all_sam_ann_position_annovar.csv", row.names=1)
dat_fin <- fread(paste0("cell_genotype_mut_read_3_all_sam.csv"))
dat_fin <- as.data.frame(dat_fin)
dat_fin<-dat_fin[,-1]
dat_fin$chr_pos_mut=paste0(dat_fin$chr_pos,"_",dat_fin$Base_observed)
tmp=pos_ann[dat_fin$chr_pos_mut,]
dat_fin=cbind(dat_fin,tmp)
write.csv(dat_fin,"cell_genotype_mut_read_3_all_sam_fin_ann_annovar.csv")
library("data.table")
dat_fin <- fread(paste0("maf_all_annovar.csv"))
dat_fin <- as.data.frame(dat_fin)
dat_fin=dat_fin[,-1]
meta <- data.combined@meta.data
meta$Tumor_Sample_Barcode=meta$name
maf_test<-dat_fin
maf_test=maf_test[-which(is.na(maf_test$Hugo_Symbol)),]
maf_test_a=maf_test[maf_test$tme_origin=="t1_a",]
maf_test_u=maf_test[maf_test$origin=="u",]
library(RColorBrewer)
mutcol=c("#1F78B4","#B2DF8A", "#A6CEE3" ,"#CAB2D6" ,"#FF7F00","black")
names(mutcol) <-c("nonsynonymous SNV","stopgain","synonymous SNV","unknown","stoploss","Multi_Hit")
library(maftools)
maf_a<-read.maf(maf_test_a,isTCGA = F,vc_nonSyn=unique(maf_test$Variant_Classification))
maf_u<-read.maf(maf_test_u,isTCGA = F,vc_nonSyn=unique(maf_test$Variant_Classification))
pdf("Summary_a.pdf",8,7)
plotmafSummary(maf = maf_a, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,color = mutcol)
dev.off()
pdf("Summary_u.pdf",8,7)
plotmafSummary(maf = maf_u, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,color = mutcol)
dev.off()
pdf("oncoplot_a.pdf",8,5.5)
oncoplot(maf_a,borderCol=NULL,
         top = 20,
         bgCol = "grey98",
         annotationDat = meta, 
         clinicalFeatures = "fincell",
        colors = mutcol,
        annotationColor = list(fincell=annocolors),
        sortByAnnotation = T
         )
dev.off()
pdf("oncoplot_a2.pdf",8,5.5)
oncoplot(maf_a,borderCol=NULL,
         top = 20,
         bgCol = "grey98", 
         annotationDat = meta,
         clinicalFeatures = "fincell",
         colors = mutcol,
         annotationColor = list(fincell=annocolors),
         sortByAnnotation = F
)
dev.off()

pdf("oncoplot_u.pdf",8,5.5)
oncoplot(maf_u,borderCol=NULL,
         top = 20,
         bgCol = "grey98", 
         annotationDat = meta, 
         clinicalFeatures = "fincell",
         colors = mutcol,
         annotationColor = list(fincell=annocolors),
         sortByAnnotation = T
)
dev.off()
pdf("oncoplot_u2.pdf",8,5.5)
oncoplot(maf_u,borderCol=NULL,
         top = 20,
         bgCol = "grey98", 
         annotationDat = meta,
         clinicalFeatures = "fincell",
         colors = mutcol,
         annotationColor = list(fincell=annocolors),
         sortByAnnotation = F
)
dev.off()

dat.titv = titv(maf = maf_a, plot = FALSE, useSyn = TRUE)
pdf("titv_a.pdf",7,6)
plotTiTv(res = dat.titv)
dev.off()
dat.titv = titv(maf = maf_u, plot = FALSE, useSyn = TRUE)
pdf("titv_u.pdf",7,6)
plotTiTv(res = dat.titv)
dev.off()

pdf("lollipopPlot_HLA-C.pdf",6,4.5)
lollipopPlot(
  maf = maf_a,
  gene = 'HLA-C',
  AACol = 'Protein_Change',
  showMutationRate = TRUE,
  colors = mutcol
)
dev.off()




