library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(pheatmap)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

#fig4b
#data.combined is the seuratObj is this study
celltype=c("T_cell","B_cell",
           "Mono_Marco","Mast_cell",
           "Fib",
           "Endo",
           "Epi_Malignant")

celltype_col=c("#BBD4E0",
      "#60AFB3",
      "#AC9ABE",
      "#462672",
      "#D9A833",
      "#79992F",
      "#962120")
names(celltype_col)= celltype

pdf("fig5d.pdf",5,5)
DimPlot(data.combined, reduction = 'umap', label=T,raster=FALSE,cols =celltype_col )
dev.off()

#fig4c
heatmap_gene <- c("CD3D","CD3E","TRBC2",
                  "CD79A","JCHAIN","IGKC",
                  "CD68","C1QA", "C1QB",
                  "CPA3","GATA2","TPSAB1",
                  "DCN","COL1A1","COL1A2",
                  "VWF","CDH5","SELE",
                  "KRT5","KRT13","KRT4"
                  
                  
)
data.combined@meta.data$cluster <- factor(data.combined@meta.data$cluster, levels = cluster_sort_by_celltype)
Idents(data.combined)="celltype"
heatmap_AveE <- AverageExpression(data.combined, assays = "RNA", features = heatmap_gene,verbose = TRUE) %>% .$RNA
gene_num <- c(3,3,3,3,3,3,3)
gaps_row <- cumsum(gene_num)
cluster_num <- c(1,1,1,1,3,2,3)
gaps_col <- cumsum(cluster_num)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
annotation_row <- data.frame(row.names = rownames(heatmap_AveE),
                             `GeneType` = rep(as.factor(c(1:length(gene_num))),gene_num))
annotation_col <- data.frame(row.names = colnames(heatmap_AveE),
                             `CellType` = colnames(heatmap_AveE))
annotation_colors = list(`CellType` = celltype_col)
pdf("fig5c_1.pdf",8,5)
pheatmap(heatmap_AveE,cluster_cols = F,cluster_rows = F,show_colnames=F,show_rownames=T,
                border=F,#border_color = "white",
                color = c(colorRampPalette(colors = c("#2166ac","#f7fbff"))(length(bk)/2),
                          colorRampPalette(colors = c("#f7fbff","#b2182b"))(length(bk)/2)),
                breaks=bk,scale="row",legend_breaks=seq(-2,2,2),
                gaps_row = gaps_row,
                gaps_col = gaps_col,
                #annotation_row = annotation_row,
                annotation_col = annotation_col,
                annotation_colors = annotation_colors,
                annotation_names_row = F,annotation_names_col = T)
dev.off()

library(reshape2)
library(ggpubr)
pB2_df <- table(data.combined@meta.data$cluster,data.combined@meta.data$group) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
pB2_df$Cluster <- factor(pB2_df$Cluster,levels = cluster_sort_by_celltype)
for (i in unique(pB2_df$Sample)) {
  pB2_df[pB2_df$Sample==i,]$Number=pB2_df[pB2_df$Sample==i,]$Number/sum(pB2_df[pB2_df$Sample==i,]$Number)
}
sample_color <- c("#1CC5FE","#FB7D80")
pB2 <- ggplot(data = pB2_df, aes(x = Cluster, y = Number, fill = Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  rotate_x_text()+
  #scale_y_continuous(position = "right",labels = percent)+
  scale_y_continuous(position = "right")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))
stat_raw=data.frame(data.combined$celltype,data.combined$time)
table1=table(stat_raw$data.combined.celltype,stat_raw$data.combined.time)
table2=table1[,-5]
num_neot=sum(table1[,5])
num_other= sum(table2)  
ratio_neo=table1[,5]/num_neot
ratio_other=table2/num_other
fin_ratio=ratio_neo/ratio_other
AveExpression=data.frame(fin_ratio)
AveExpression$Gene <- "fin_ratio"
AveExpression$Cluster=rownames(AveExpression)
Ave_df=AveExpression
colnames(Ave_df) <- c("Expression", "Gene", "Cluster")
Ave_df$Group <- paste(Ave_df$Cluster,Ave_df$Gene,sep = "_")
Ave_df$Expression[which(Ave_df$Expression>3)] <- 3 Ave_df$Cluster=factor(Ave_df$Cluster,levels = celltype)
pB1 <- ggplot(Ave_df,aes(x=Gene, y=Expression))+
  geom_hline(yintercept = seq(0, 2, 0.5),linetype = 2, color = "lightgray",size=1)+
  geom_line()+
  geom_segment(aes(x=Gene,xend=Gene,y=0,yend=Expression),color="lightgray",size = 1.5)+
  geom_point(size=4,aes(color=Gene))+
  scale_color_manual(values=c("#FB7D80","#AF1E1F")) +
 # scale_color_continuous(values=c("#00AFBB", "#E7B800", "#FC4E07", "#41ab5d"))+
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="ratio")
pB1 <- facet(pB1, facet.by = "Cluster",ncol = length(unique(Ave_df$Cluster)),panel.labs.font = list(size = 12),panel.labs.background = list(fill = "#a6cee3"))
pB1 <- pB1 + scale_y_continuous(position = "right")+ 
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_blank())+ 
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "right",
        panel.border = element_blank(),
        axis.ticks.x = element_line(color =  NA))
require(ggplotify)
pB_1_2 <- pB1 + pB2 + plot_layout(ncol = 1, heights = c(1.2, 2))
ggsave(paste0("fig5c_2.pdf"), plot = pB_1_2, width = 8, height =5)

#fig4d
library(ggplot2)
library(dplyr)
library(reshape2)
bar_plot_input$N=-bar_plot_input$N
x=melt(bar_plot_input,id.vars = "cell")
x$ylab=x$cell
x$ylab=factor(x$ylab,levels = bar_plot_input[order(bar_plot_input$T),]$cell)
x$variable=factor(x$variable,levels = c("T","N"))
x$variable=factor(x$variable,levels = c("T","N"))
col <- c("#FB7D80","#1CC5FE")
ggplot(x,aes(x=ylab,y=value,fill=variable))+
  geom_col(position=position_dodge(0),width = 1)+
  coord_flip()+
  ylim(-1.65,1.65)+
  labs(x=NULL,y=NULL)+
  theme_bw()+
  scale_fill_manual(guide = guide_legend(title = NULL),values = col)+
  theme(legend.position ="top", 
        legend.key.size=unit(0.8,"line"),
        axis.text.y = element_text(color =celltype_col[levels(x$ylab)],size = 12))
 
#fig4e
sample_color <- c("#1CC5FE","#FB7D80","#39398A","#AF1E1F")
library(stringr)
data.combined$sam_group=paste0(data.combined$sam,":",data.combined$group)
cellratio <- prop.table(table(data.combined$celltype, data.combined$sam_group), margin = 2)%>%
  as.data.frame()%>%
  `colnames<-`(c('celltype','orig.ident','Freq'))
write.csv(cellratio,"cellratio.csv")
group=strsplit(as.vector(cellratio$orig.ident),":")
group=lapply(group,FUN = function(x){x[2]})
group=unlist(group)
cellratio$group <-group
cellratio1 <- cellratio%>%
  group_by(celltype, group) %>%
  summarise(avg.freq = mean(Freq),
            sd.freq = sd(Freq),
            se.freq = sd(Freq)/sqrt(n()))
cellratio1$celltype=factor(cellratio1$celltype,levels = celltype)
cellratio$celltype=factor(cellratio$celltype,levels = celltype)

library(ggpubr)
for (i in unique(cellratio$celltype)) {
  cellratio_sub=cellratio[cellratio$celltype==i,]
  cellratio1_sub <- cellratio_sub%>%
    group_by(celltype, group) %>%
    summarise(avg.freq = mean(Freq),
              sd.freq = sd(Freq),
              se.freq = sd(Freq)/sqrt(n()))  
  ebtop<-function(x){
    return(mean(x)+sd(x)/sqrt(length(x)))
  }
  ebbottom<-function(x){
    return(mean(x)-sd(x)/sqrt(length(x)))
  }
  
  c <- ggplot(data=cellratio_sub,aes(x=group,y=Freq,fill=group))+
    stat_summary(geom = "bar",
                 fun = mean,
                 )+
    stat_summary(geom = "errorbar",
                 fun.min = ebbottom,
                 fun.max = ebtop,
                 width=0.2)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10))  + 
    scale_fill_manual(values=sample_color)+
    scale_color_manual(values=sample_color)+
    geom_jitter(aes(fill=group),width = 0.2,size=1,shape=21,color="black") 
  ggsave(plot = c,paste0("box1_",i,".pdf"), width=3.8, height=4.3)
  p <- ggplot(cellratio_sub,
              aes(x=group, y=Freq, 
                  color=group
              )) + 
    geom_boxplot(position = position_dodge(0.8)) +
    geom_point(position = position_jitterdodge()) +
    scale_color_manual(values = sample_color) +
    theme_bw()+#+stat_compare_means() 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10))  + 
    stat_compare_means(comparisons =my_comparisons, simplify =FALSE)
  ggsave(plot = p,paste0("box2_",i,".pdf"), width=3.8, height=4)
}

#fig4f
#data.combined_tmp is the seuratObj of a specific group in this study
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpointdensity) 
sce_cca=data.combined_tmp
data <- cbind(Embeddings(object=sce_cca[['umap']]),
              FetchData(sce_cca,c("celltype","group")))
#plot
ggplot(data = data, mapping = aes(x = UMAP_1,y = UMAP_2)) + 
  geom_pointdensity(adjust=6) +
  geom_density_2d(bins = 5, colour="black")+
   theme_bw()+
  theme(
       axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    legend.text = element_text(size =10),
    aspect.ratio=1,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size=12, color = "black",
                                vjust = 0.5,margin = margin(b = 3,t=3)), 
    strip.background = element_blank(),
    plot.margin=unit(c(1, 1, 1, 1),'cm'),
    legend.margin = margin(-0.2,-0.2,0,0,'cm'),
    legend.key.height = unit(1,'cm'),
    legend.key.width = unit(0.4,'cm'))+
  scale_color_gradientn(colours = c('#5749a0', '#0f7ab0', '#00bbb1',
                                    '#bef0b0', '#fdf4af', '#f9b64b',
                                    '#ec840e', '#ca443d', '#a51a49'), 
                        guide = guide_colorbar(frame.colour = "black", 
                                               ticks.colour = NA),
                        name = "Density",
                        labels = c('low',"high"),
                        breaks = c(100,700)) 
ggsave(“fig5e_1.pdf”, width = 4, height = 4)

cell.info=data.combined@meta.data
cell.info$UMAP_1=data.combined@reductions$umap@cell.embeddings[,1]
cell.info$UMAP_2=data.combined@reductions$umap@cell.embeddings[,2]
cell.info$fincell=as.vector(cell.info$fincell)

col2=c("grey90")
names(col2)="Others"
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

cell.info$Tissue=cell.info$group
cell.info$ClusterID=cell.info$celltype
cell.info$CellType=cell.info$celltype
p<-umapPlot.tissue(cell.info, "selected group", 25)
ggsave("selected group.pdf",plot = p,width = 6,height = 4.6)

library(Seurat)
library(ggplot2)
library(limma)
degs <- FindMarkers(object = data.combined, 
                       only.pos = TRUE,
                       ident.1="B",
                       ident.2="A",
                       logfc.threshold = 0.1) 

diff=degs
diff$P.Value=diff$p_val
diff$logFC=diff$avg_log2FC
write.csv(diff,"diff_B_vs_A.csv")
diff_sub=diff[diff$P.Value<0.05,]
diff_sub=diff_sub[diff_sub$logFC>0.05,]
degs_all=unique(rownames(diff_sub))
sigScores <- as.matrix(data.combined@assays$RNA@counts[degs_all,])
sigScores_degs=sigScores[degs_all,]
sigScores_degs[which(sigScores_degs>0)]<-1
degs_express=colSums(sigScores_degs)
meta=data.combined@meta.data
meta=meta[names(degs_express),]
meta$degs_express=degs_express
library("viridis")
UMap = data.frame(data.combined@reductions$umap@cell.embeddings)
UMap$TranscriptCount = meta$degs_express
UMap$Cluster = meta$fincell
max(UMap$TranscriptCount)
UMap$TranscriptCount[which(UMap$TranscriptCount>435)]<-435
ggplot(UMap[order(UMap$TranscriptCount),],aes(x=UMAP_1, y=UMAP_2, colour=TranscriptCount)) + 
  geom_point(alpha=0.8, size=0.9) +
  xlab("UMAP 1") + 
  ylab("UMAP 2") + 
  theme_bw()+
  theme(panel.border = element_rect(fill=NA),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        aspect.ratio=1,  # width : height =1 
        legend.position = 'right',
        legend.key = element_rect(size=1),# legend size
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
  ) + 
  labs(colour="Number of genes")+
  scale_color_viridis(option = "C", direction = 1,limits = c(0, 435))

all_fin=data.frame()
for (sam in unique(meta$sample)) {
  meta_sub1=meta[which(meta$sample==sam),]
  for (time in c("all")) {
    meta_sub2=meta_sub1[which(meta_sub1$time==time),]
    
    
    roe_celltype_u=c()
    roe_celltype_a=c()
    p=c()
    
    for (i in unique(meta$celltype)) {
      tmp=meta_sub2
      tmp_sub1=tmp[which(tmp$celltype==i),]
      tmp_sub1$chi_id="in_cell"
      tmp_sub2=tmp[which(tmp$celltype!=i),]
      tmp_sub2$chi_id="out_cell"
      tmp_fin=rbind(tmp_sub1,tmp_sub2)
      ka<- xtabs(~tmp_fin$chi_id+tmp_fin$origin,data=tmp_fin)
      chi=chisq.test(ka)
      roe=chi$observed/chi$expected
      roe=as.data.frame(roe)
      p=c(p,chi$p.value)
      roe_celltype_u=c(roe_celltype_u,roe[3,3])
      roe_celltype_a=c(roe_celltype_a,roe[1,3])
    }
    names(roe_celltype_u)=unique(meta$celltype)
    names(roe_celltype_a)=unique(meta$celltype)
    names(p)=unique(meta$celltype)
    fin=data.frame(roe_celltype_u,roe_celltype_a,p)
    fin$time=time
    fin$sample=sam
    fin$cell=unique(meta$fincell)
    all_fin=rbind(all_fin,fin)
  }
}

data=all_fin[all_fin$time=="all",]
data$pstar <- ifelse(data$p < 0.05,
                     ifelse(data$p < 0.01,"**","*"),
                     "")
data$pstar[1:20]
data$id=factor(data$cell,levels = rev(celltype))
celltype_col=rev(celltype_col)
data$variable=factor(data$sample)
data$col<-cut(data$roe_celltype_u,
              breaks = c(0,0.5,1,1.5,3))
data_sub=data

ggplot(data_sub, aes(variable,id,fill=col)) + 
  geom_tile(colour = "white",size=1.5)+

  scale_fill_manual(breaks = levels(data$col),
                    values =c("white","#CAE4ED","#E88860","#C6483D","#C6483D"))+
  geom_text(aes(label=pstar),col ="black",size = 4)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
        axis.text.y = element_text(size = 10,color =celltype_col))+
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","roe"))









