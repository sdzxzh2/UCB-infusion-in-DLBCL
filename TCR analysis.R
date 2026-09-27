
library(Seurat)
library(scRepertoire)
slot(seurat1, "meta.data")$cloneType <- factor(slot(seurat1, "meta.data")$cloneType, 
                                               levels = c("Hyperexpanded (100 < X <= 500)", 
                                                          "Large (20 < X <= 100)", 
                                                          "Medium (5 < X <= 20)", 
                                                          "Small (1 < X <= 5)", 
                                                          "Single (0 < X <= 1)", NA))
pdf("cloneType.pdf",6,5)
DimPlot(seurat1, group.by = "cloneType") +
scale_color_manual( values =  colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6","#0348A6","#0348A6"))(5),na.value="grey")
dev.off()

cell.info=meta
cell.info$celltype=paste0(cell.info$fincell,"_",cell.info$origin)
cell.info$celltype=cell.info$tme_origin
celltypes <- celltype
celltypes <- unique(cell.info$celltype)
count <- lapply(celltypes, function(A){
  sapply(celltypes, function(B){
    length(intersect(
      cell.info[["clone_type_all_sam"]][cell.info[["celltype"]] == A],
      cell.info[["clone_type_all_sam"]][cell.info[["celltype"]] == B]
    ))
  })
})
count <- do.call(rbind, count)
rownames(count) = colnames(count) = celltypes
write.table(count, "output_count.txt", row.names = T, col.names = T, sep = "\t", quote = F)
library(pheatmap)
library(ComplexHeatmap)
upper <- 6 
plot.order <-celltype[1:11]
count <- read.table("output_count.txt", header = T, sep = "\t", check.names = F)
count <- log(1 + count)
count[count > upper] = upper
count <- count[plot.order, plot.order]
hm<-pheatmap(as.matrix(count), 
               color = colorRampPalette(c("#440052", "#259290", "#F8E725"))(100),
               border_color = NA,
               show_colnames = T, 
               show_rownames = T,
               cellwidth = 12, 
               cellheight = 12,
               name = "log(1+Number of\nTCR clonaltypes)")

pdf(file = "count_heatmap.pdf", width = 8, height = 8)
draw(hm)
invisible(dev.off())

meta=seurat1@meta.data
meta=meta[-which(is.na(meta$CTstrict)),]
meta_fin=data.frame()
for (i in unique(meta$orig.ident)) {
  meta_sub=meta[meta$orig.ident==i,]
  clone_tab=as.data.frame(table(meta_sub$CTstrict))
  clone_tab$clone="no_clone"
  if (length(clone_tab[clone_tab$Freq>1,]$clone)>0) {
    clone_tab[clone_tab$Freq>1,]$clone<-"clone"
    clone_vdj=clone_tab[clone_tab$clone=="clone",]
    meta_sub$clone="no_clone"
    meta_sub[meta_sub$CTstrict%in%clone_vdj$Var1,]$clone<-"clone"
  }else{ meta_sub$clone="no_clone"}
  meta_sub$clone_count=1
  for (t in 1:nrow(meta_sub)) {
    meta_sub2=meta_sub[t,]
    count<-clone_tab[clone_tab$Var1==meta_sub2$CTstrict,]$Freq
    meta_sub[t,"clone_count"]<-count
  }
  clone_tab=clone_tab[rev(order(clone_tab$Freq)),]
  clone_tab$clone_type_each_sam=paste0(i,"_clonetype_",1:nrow(clone_tab))
  meta_sub$clone_type_each_sam="clonetype"
  for (t in 1:nrow(meta_sub)) {
    meta_sub3=meta_sub[t,]
    clonetype<-clone_tab[clone_tab$Var1==meta_sub3$CTstrict,]$clone_type_each_sam
    meta_sub[t,"clone_type_each_sam"]<-clonetype
  }
  meta_fin=rbind(meta_fin,meta_sub)
}
meta=meta_fin
clonetype_tab=as.data.frame(table(meta$CTstrict))
clonetype_tab=clonetype_tab[rev(order(clonetype_tab$Freq)),]
clonetype_tab$clone_type_all_sam=paste0("all_clonetype_",1:nrow(clonetype_tab))
meta$clone_type_all_sam="clonetype"
for (i in 1:nrow(meta)) {
  meta_sub=meta[i,]
  type<-clonetype_tab[clonetype_tab$Var1==meta_sub$CTstrict,]$clone_type_all_sam
  meta[i,"clone_type_all_sam"]<-type
}
write.csv(meta,"meta_fin1.csv")
meta_fin2=data.frame()
for (i in unique(meta$patient)) {
  meta_sub=meta[meta$patient==i,]
  clone_tab= as.data.frame(table(meta_sub$CTstrict))
  clone_tab=clone_tab[rev(order(clone_tab$Freq)),]
  clone_tab$clone_type_each_patient=paste0(i,"_clonetype_",1:nrow(clone_tab))
  meta_sub$clone_type_each_patient="clonetype"
  for (t in 1:nrow(meta_sub)) {
    meta_sub3=meta_sub[t,]
    clonetype<-clone_tab[clone_tab$Var1==meta_sub3$CTstrict,]$clone_type_each_patient
    meta_sub[t,"clone_type_each_patient"]<-clonetype
  }
   meta_fin2=rbind(meta_fin2,meta_sub)
 }
write.csv(meta_fin2,"meta_fin2.csv")

library(dplyr)
meta_data=data.combined@meta.data
sam=intersect(rownames(meta_data),rownames(meta_fin2))
meta_fin2=meta_fin2[rownames(meta_fin2)%in%sam,]
sam=as.vector(sam)
meta_data=meta_data[!rownames(meta_data)%in%sam,]
meta=bind_rows(meta_data,meta_fin2)
meta=meta[rownames(data.combined@meta.data),]
data.combined@meta.data=meta

cell.info=data.combined@meta.data
cell.info$clone_group=paste0(cell.info$clone,"_",cell.info$group)
cell.info$UMAP_1=data.combined@reductions$umap@cell.embeddings[,1]
cell.info$UMAP_2=data.combined@reductions$umap@cell.embeddings[,2]
cell.info$fincell=as.vector(cell.info$fincell)
col2=c("grey90")
names(col2)="Others"
celltype_col_raw=celltype_col
celltype_col=c(celltype_col,col2)
umapPlot.tissue <- function(data, tissue.use, topn) {
  clusters.use <- table(subset(data, Tissue == tissue.use)$ClusterID)
  if(topn > length(clusters.use)) {
    topn <- length(clusters.use)
  }
  clusters.use <- names(head(sort(clusters.use, decreasing = T), n=topn))
  data$new.cluster <- ifelse(data$Tissue == tissue.use & data$ClusterID %in% clusters.use, data$CellType, "Others")
  data$new.cluster <- factor(data$new.cluster, levels = c(setdiff(names(table(data$new.cluster)), "Others"), "Others"))
  data$pt.size <- ifelse(data$new.cluster == "Others", 0.4, 0.4)
  data1=data[data$new.cluster=="Others",]
  data2=data[data$new.cluster!="Others",]
  ggplot() + 
    geom_point(data = data1,
               aes(UMAP_1, UMAP_2, color=new.cluster),
               size=data1$pt.size,
               alpha=0.4) + 
    geom_point(data =data2,
               aes(UMAP_1, UMAP_2, color=new.cluster),
               size=data2$pt.size,
               alpha=0.5)+
    scale_color_manual(values =celltype_col) + 
    guides(colour = guide_legend(keyheight = .7, 
                                 keywidth = .1, 
                                 override.aes = list(size=3))) + 
    theme_bw(base_size = 12) + 
    ggtitle(tissue.use) + 
    theme(#legend.justification = c(0,0), 
      legend.title = element_blank(),
      legend.key = element_rect(fill = alpha("white", 0)),
      legend.background = element_rect(fill=alpha('white', 0)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color="black"),
      plot.title = element_text(hjust = .5, face = "bold")
    )
  
}
cell.info$clone_sam=cell.info$clone_group
cell.info$Tissue=cell.info$clone_group
cell.info$ClusterID=cell.info$fincell
cell.info$CellType=cell.info$fincell
no_clone_group=paste0("no_clone","_",unique(cell.info$group))
clone_group=paste0("clone","_",unique(cell.info$group))
for (i in c(no_clone_group,clone_group)) {
  p<-umapPlot.tissue(cell.info, i, 25)
  ggsave(paste0(i,"_type1.pdf"),plot = p,width = 6.5,height = 5)
}
celltype_col=rep("#3A7BBF",length(celltype))
names(celltype_col)=celltype
celltype_col=c(celltype_col,col2)
for (i in c(no_clone_group)) {
  p<-umapPlot.tissue(cell.info, i, 25)
  ggsave(paste0(i,"_type2_A_no_clone.pdf"),plot = p,width = 6.5,height = 5)
}

celltype_col=rep("#E6315E",length(celltype))
names(celltype_col)=celltype
celltype_col=c(celltype_col,col2)
for (i in c(clone_group)) {
  p<-umapPlot.tissue(cell.info, i, 25)
  ggsave(paste0(i,"_type2_B_clone.pdf"),plot = p,width = 6.5,height = 5)
}
meta_tcr=meta_tcr_fin
meta_tcr$Frequency=meta_tcr$Freq_sam_origin
meta_tcr=meta[!is.na(meta$clone),]
meta_post=meta_tcr[meta_tcr$treat=="post",]
meta_post_clone=meta_post[meta_post$clone=="clone",]
clone_type=unique(meta_post_clone$clone_type_each_patient)
clonetype_emergent_in_treat=c()
for (i in clone_type) {
  meta_sub=meta_tcr[meta_tcr$clone_type_each_patient==i,]
  if (length(unique(meta_sub$treat))==1) {
    clonetype_emergent_in_treat=c(clonetype_emergent_in_treat,i)
  }else{
    dat=data.frame(meta_sub$treat,meta_sub$Frequency)
    colnames(dat)<-c("group","freq")
    dat=aggregate(freq ~ group, data = dat, FUN = mean, na.action = na.omit)
    if (dat[dat$group=="post",]$freq>dat[dat$group=="before",]$freq) {
      clonetype_emergent_in_treat=c(clonetype_emergent_in_treat,i)
    }
  }
  
}

meta$clone_emergent="no"
meta[which(meta$clone=="clone"),]$clone_emergent<-"no_emergent"
meta[meta$clone_type_each_patient%in%clonetype_emergent_in_treat,]$clone_emergent<-"emergent"
write.csv(meta,"meta_emergent.csv")
cell.info=meta
cell.info$clone_emergent_patient=paste0(cell.info$clone_emergent,"_",cell.info$patient)
cell.info$UMAP_1=data.combined@reductions$umap@cell.embeddings[,1]
cell.info$UMAP_2=data.combined@reductions$umap@cell.embeddings[,2]
cell.info$fincell=as.vector(cell.info$fincell)
cell.info$clone_sam=cell.info$clone_emergent
cell.info$Tissue=cell.info$clone_emergent
cell.info_3sam=cell.info[cell.info$patient%in%c("LXD","XL","XYH"),]
cell.info$ClusterID=cell.info$fincell
cell.info$CellType=cell.info$fincell
celltype_col=rep("#3A7BBF",length(celltype))
names(celltype_col)=celltype
celltype_col=c(celltype_col,col2)
p<-umapPlot.tissue(cell.info_3sam, "emergent", 25)
ggsave(paste0("clonetype_emergent.pdf"),plot = p,width = 6.5,height = 5)

my36colors <- c('#53A85F', '#E63863', '#57C3F3', '#68A180', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

cell.info_3sam$ClusterID=cell.info_3sam$patient
cell.info_3sam$CellType=cell.info_3sam$patient

celltype_col=my36colors[1:length(unique(cell.info_3sam$patient))]
names(celltype_col)=unique(cell.info_3sam$patient)
celltype_col=c(celltype_col,col2)

p<-umapPlot.tissue(cell.info_3sam, "emergent", 25)
ggsave(paste0("clonetype_emergent_time.pdf"),plot = p,width = 5.5,height = 5)


library(ggplot2)
library(jjPlot)
df.long <- reshape2::melt(pie_stat,id.vars = c('group','fincell','group_fincell'),
                          variable.name = 'type',value.name = 'per')
df=df.long
df=df[df$fincell%in%celltype,]
df$fincell=factor(df$fincell,levels = celltype)
p<-ggplot(df,aes(x = group,y = fincell,group = group_fincell)) +
  geom_jjPointPie(aes(pievar = per,
                      fill = type),
                  width = 1.2,
                  color = 'white',
                  line.size = 1) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,color =sample_color, hjust = 1),
        axis.text.y = element_text(color =celltype_col, hjust = 1)) +
  ggsci::scale_fill_npg()+ 
  geom_vline(xintercept=c(1.5,3.5,5.5), linetype="dotted",size=1.2)

seurat1$group2=paste0(seurat1$group,"_",seurat1$fincell)
seurat1$treat_response=paste0(seurat1$treat,"_",seurat1$response)
xx=StartracDiversity(seurat1, type = "treat_response", 
                     sample = "treat_response", 
                     by = "overall")#+scale_fill_manual(values =celltype_col)

write.csv(xx$data,"StartracDiversity_test5.csv")

library(ggpubr)
plot=plot_pair
plot$ID=plot$aid
ggplot(plot, aes(x = group, y = value, color = group)) +
  geom_line(aes(group = id), color = "grey80",alpha=0.5, size = 1) +
  geom_boxplot(size=1,alpha=0.5,width =0.5) +
  geom_point(size = 2) +
  stat_compare_means() +
  scale_fill_manual(values=sample_color)+
  scale_color_manual(values=sample_color)+
  theme_classic2(base_size = 15) +
  xlab('group')







