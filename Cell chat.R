library(CellChat)
library(patchwork)

mycellchat=cellchat
LR <- mycellchat@LR$GROUP$LRsig
write.csv(LR, file = 'LR.csv')
Target=LR_sec$receptor
Target=strsplit(Target,"_")
target_fin=c()
for (i in 1:length(Target)) {
  tmp=Target[[i]]
  target_fin=c(target_fin,tmp[1])
}

target_raw=LR_sec$receptor
target_fin=target_fin

target_data=data.frame(target_raw,target_fin)
receptor_dot <- DotPlot(data.combined, features = unique(target_fin),assay = "RNA")
target_fin_fin<-unique(receptor_dot$data$features.plot)
target_data=target_data[target_data$target_fin%in%target_fin_fin,]
tile=LR_sec
tile=tile[tile$receptor%in%target_data$target_raw,]
tile_raw=tile
Target=tile_raw$ligand
Target=strsplit(Target,"_")
target_fin=c()
for (i in 1:length(Target)) {
  tmp=Target[[i]]
  target_fin=c(target_fin,tmp[1])
}
target_raw=tile_raw$ligand
target_fin=target_fin
target_data=data.frame(target_raw,target_fin)
receptor_dot <- DotPlot(data.combined, features = unique(target_fin),assay = "RNA")
target_fin_fin<-unique(receptor_dot$data$features.plot)
target_data=target_data[target_data$target_fin%in%target_fin_fin,]
tile=tile[tile$ligand%in%target_data$target_raw,]

LR_sec=tile
colnames(LR_sec)
LR_sec <- LR_sec[,c("ligand","receptor","pathway_name")]
receptor <- unique(LR_sec$receptor)
receptor
length(receptor)
LR_sec$Target <- LR_sec$receptor
for (i in 1:length(receptor)){
  LR_sec[, "receptor"][LR_sec[, "receptor"]==receptor[i]] = i
  
}
ligand <- unique(LR_sec$ligand)
ligand
length(ligand)
LR_sec$source <- LR_sec$ligand

for (i in 1:length(ligand)){
  LR_sec[, "ligand"][LR_sec[, "ligand"]==ligand[i]] = i
  
}
library(GGally)
library(dittoSeq)
library(ggplot2)
LR_sec$ligand <- factor(LR_sec$ligand, levels = c(1:length(unique(LR_sec$ligand))))
LR_sec$receptor <- factor(LR_sec$receptor, levels = c(1:length(unique(LR_sec$receptor))))
ylab_right <- c(unique(LR_sec$Target),rep("",length(unique(LR_sec$source))-length(unique(LR_sec$Target))))
ylab_right <- c(unique(LR_sec$source),rep("",length(unique(LR_sec$Target))-length(unique(LR_sec$source))))

ylab_left_color <- LR_sec[!duplicated(LR_sec$source), ]
tab=table(ylab_left_color$pathway_name)
tab=data.frame(tab)
rownames(tab)=tab[,1]
path_name=unique(ylab_left_color$pathway_name)
col_all <- c("#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2","#D55E00",
             '#E5D2DD', '#53A85F', '#F1BB72', '#E95C59', '#D6E7A3', '#57C3F3', '#476D87',
             '#F3B1A0', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
             '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
             '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
             '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
             '#968175')
tab=tab[path_name,]
col=col_all[1:nrow(tab)]
col=rep(col,tab[,2])


ylab_left_color$col=col
col_unique=col_all[1:nrow(tab)]
names(col_unique)=rownames(tab)
col_ligand=col
names(col_ligand)=ylab_left_color$ligand

if (length(unique(LR_sec$source))-length(unique(LR_sec$Target))>0) {
  p11 <- LR_sec %>% 
    ggparcoord(
      columns = 1:2, 
      groupColumn = 3, 
      showPoints = F, 
      alphaLines = 0.6,scale="globalminmax",
      splineFactor=T)+
    scale_y_continuous(breaks = c(1:length(unique(LR_sec$source))),
                       labels =unique(LR_sec$source),
                       sec.axis = sec_axis(~./1, 
                                           breaks=c(1:length(unique(LR_sec$source))),
                                           labels=ylab_right))+
    theme(legend.position = 'none')+
    theme(plot.margin = margin(t = 0,  
                               r = 0,  
                               b = 0,  
                               l = 0, 
                               unit = "cm"))+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          legend.position = 'none',
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y.left = element_text(colour =col))+
    scale_color_manual(values = col_unique)  
}else{
  p11 <- LR_sec %>% 
    ggparcoord(
      columns = 1:2, 
      groupColumn = 3,
      showPoints = F, 
      alphaLines = 0.6,scale="globalminmax",
      splineFactor=T)+
    scale_y_continuous(breaks = c(1:length(unique(LR_sec$Target))),
                       labels =ylab_right,
                       sec.axis = sec_axis(~./1, 
                                           breaks=c(1:length(unique(LR_sec$Target))),
                                           labels=unique(LR_sec$Target)))+
    theme(legend.position = 'none')+
    theme(plot.margin = margin(t = 0,  
                               r = 0,  
                               b = 0,  
                               l = 0,  
                               unit = "cm"))+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          legend.position = 'none',
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y.left = element_text(colour =col))+
    scale_color_manual(values =col_unique)
  
}

library(Seurat)
data.combined=data.combined
Idents(data.combined)<-"group"
HC <-data.combined
ligand_dot <- DotPlot(HC, features = ligand,assay = "RNA")
ligand_exp <- ligand_dot$data
unique(ligand_exp$id)
group_level=unique(as.vector(sort(as.vector(ligand_exp$id))))
ligand_exp$id <- factor(ligand_exp$id, levels = celltype)
ligand_exp$id <- factor(ligand_exp$id, levels = group_level)
write.csv(ligand_exp,"ligand_exp_dc_cell.csv")
p2 <- ggplot(ligand_exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  theme_bw()+
  scale_x_discrete(position = "top")+
  scale_y_discrete(position = "right")+
  theme(axis.text.x=element_text(angle=90,hjust = 0,vjust=1,colour = 'black',size = 9),
        axis.text.y = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.3,'cm'),
        legend.key.width = unit(0.3,'cm'),
        legend.title = element_text(size=8,vjust = 1),
        legend.text = element_text(size = 5))+
  theme(plot.margin = margin(t = 0, 
                             r = 0,  
                             b = 0,  
                             l = 0,  
                             unit = "cm"))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#1A5592','white',"#B83D3D"),
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"),
                        name = "Scaled expr")+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(title = "",order = 1,nrow = 2))+
  theme(plot.margin=unit(c(1, 1, 1, 1),'cm'))+
  coord_cartesian(clip = 'off') 
p4 <- p2
p4

tile[tile$prob>0.01,]$prob<-0.01
p_tile<-ggplot(data=tile,aes(x=L,y= R)) + 
  geom_tile(data=tile,aes(x=L,y= R,  fill=prob),color="white",size=0)+
    scale_fill_gradient2(high="#9733CF", low="white")+#对应修改有填充的点颜色
    theme_bw()+
  theme(legend.position = "none",
   plot.margin = margin(t = 0,  
                         r = 0,  
                         b = 0,  
                         l = 0,  
                         unit = "cm"),
    panel.grid = element_blank(),
    axis.text.x=element_text(angle=90,hjust = 0,vjust=1),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
  coord_cartesian(clip = 'off')
p_tile
library(foreach)
library(doParallel)
library(data.table)
cl <- makeCluster(detectCores()-3)
registerDoParallel(cl)
sam=c("Untitled.txt")
data=data.frame()
for (i in sam) {
  tmp <- fread(paste0("./",i),check.names = F,stringsAsFactors = F,header = T,data.table = F) 
  study=gsub(".txt","",i)
  tmp$study=study
  tmp[8]=tmp[8]-min(tmp[8])
  tmp[9]=tmp[9]-min(tmp[9])
  data=rbind(data,tmp)
  }
colnames(data)[8]<-"Centroid_X"
colnames(data)[9]<-"Centroid_y"
type_marker<-foreach(i=1:nrow(data),.combine = "rbind") %dopar% {
  data_sub=data[i,]
  Class_raw=data_sub$Classification
  Class_raw=gsub(" ","",Class_raw)
  Class_raw=strsplit(Class_raw,":")
  Class_raw=Class_raw[[1]]
  type_marker=c(rep(0,6))
colnames(type_marker)<-c("CXCR5","CD8","CD4","CD20","CD21","CXCL13")
rownames(type_marker)<-rownames(data)
data=cbind(data,type_marker)
save.image("all.RData")
cols=c("CD8","CD4","FDC","B_cell","CD8_CXCR5","B_cell_CXCR5","CD4_CXCL13","Tumor_Stroma")
data_fin=data[,c(cols,"Centroid_X","Centroid_y","study")]
celltype<-foreach(i=1:nrow(data_fin),.combine = "c") %dopar%  {
  data_sub=data_fin[i,]
  celltype="a_Tumor_Stroma"
  if (sum(data_sub[,1:6])==0) {
    celltype<-"a_Tumor_Stroma"
  }
  
  if (sum(data_sub[,1])>0) {
    celltype<-"b_CD8"
  } 
  
  if (sum(data_sub[,2])>0) {
    celltype<-"c_CD4"
  } 
  
  if (sum(data_sub[,3])>0) {
    celltype<-"e_FDC"
  } 
  
  if (sum(data_sub[,4])>0) {
    celltype<-"f_B_cell"
  } 
  
  if (sum(data_sub[,5])>0) {
    celltype<-"d_CD8_CXCR5"
  } 
  
  if (sum(data_sub[,6])>0) {
    celltype<-"g_B_cell_CXCR5"
  } 
  
  if (sum(data_sub[,7])>0) {
    celltype<-"h_CD4_CXCL13"
  } 
  celltype
}
data_fin$celltype=celltype
library(ggvoronoi)
library(ggplot2)
library(dplyr)
library(ggforce)
col=c("#C6E0C0","#BBD4E0","#8CB4D2","#6B68CB","#DE4C00")
celltype=sort(unique(data_fin$celltype))
names(col)= celltype
celltype_col=col
sam=gsub(".txt","",sam)
data_fin_raw=data_fin
data_fin=data_fin_raw
data_fin_1=data_fin[data_fin$celltype=="a_Tumor_Stroma",]
data_fin_2=data_fin[data_fin$celltype!="a_Tumor_Stroma",]
data_fin_1=data_fin[sample(rownames(data_fin_1),10000),]
data_fin=data_fin_2
for (i in sam) {
   data_fin_sub=data_fin[data_fin$study==i,]
  data_fin_sub=data_fin_sub[order(data_fin_sub$celltype),]
  
  point<-ggplot(data=data_fin_sub,aes(x=Centroid_X,y=Centroid_y))+ geom_point(
    aes(colour=celltype, 
        fill=celltype), size = 1,alpha=1)+
    scale_color_manual(values=celltype_col) +
    scale_fill_manual(values=celltype_col)+
       theme()+theme(panel.background = element_rect(fill = 'white'),
                  panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  ggsave(paste0("point1_",i,".pdf"),plot =point, width = 80, height = 81,limitsize = FALSE)
  
  
  point<-ggplot(data=data_fin_sub,aes(x=Centroid_X,y=Centroid_y))+ geom_point(
    aes(colour=celltype, 
        fill=celltype), size = 2,alpha=1)+
    scale_color_manual(values=celltype_col) +
    scale_fill_manual(values=celltype_col)+
     theme()+theme(panel.background = element_rect(fill = 'black'),
                  panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  ggsave(paste0("point2_",i,".pdf"),plot =point, width = 80, height = 81,limitsize = FALSE)

  data_fin3=data_fin_sub
  
  map<-  ggplot(data_fin3, aes(Centroid_X, Centroid_y, group = -1L)) +
    geom_voronoi_tile(aes(fill = celltype,colour = celltype), normalize = F,expand = unit(-1.5, 'pt'), radius = unit(0, 'pt'), max.radius = 20)+
    scale_color_manual(values=celltype_col) +
    scale_fill_manual(values=celltype_col)+
    theme()+
    theme(panel.background = element_rect(fill = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) 
  ggsave(paste0("voronoi1_",i,".pdf"),plot =map, width = 250, height = 251,limitsize = FALSE)
  map<-  ggplot(data_fin3, aes(Centroid_X, Centroid_y, group = -1L)) +
    geom_voronoi_tile(aes(fill = celltype,colour = celltype), normalize = F,expand = unit(-1.5, 'pt'), radius = unit(0, 'pt'), max.radius = 20)+
    scale_color_manual(values=celltype_col) +
    scale_fill_manual(values=celltype_col)+
    theme()+
    theme(panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) 
  ggsave(paste0("voronoi2_",i,".pdf"),plot =map, width = 250, height = 251,limitsize = FALSE)  
}
stopCluster(cl) 
data_fin_for_tile=data_fin
data_fin_for_tile=data_fin_for_tile[which(data_fin_for_tile$celltype!="b_Tumor_Stroma"),]
data_fin_for_tile=data_fin_for_tile[which(data_fin_for_tile$celltype!="a_multi"),]
data_fin_for_tile=data_fin_for_tile[which(data_fin_for_tile$celltype!="TFAM_HLA_DC"),]
distance=100
min_cell=1
celltype=sort(unique(data_fin$celltype))
windows_cell_freq_all=data.frame()
for (i in sam) {
  data_fin_sub=data_fin_for_tile[data_fin_for_tile$study==i,]
  X_max=max(data_fin_sub$Centroid_X)
  Y_max=max(data_fin_sub$Centroid_y)
  x_tile=floor(X_max/distance)
  y_tile=floor(Y_max/distance)
  windows_cell_freq=data.frame()
  for (x in 1:x_tile) {
    img_sub_x=data_fin_sub[distance*(x-1)<data_fin_sub$Centroid_X & data_fin_sub$Centroid_X<distance*x,]
    for (y in 1:y_tile) {
      windows_name=paste0(i,"_",x,"_",y)
      img_sub_x_y=img_sub_x[distance*(y-1)<img_sub_x$Centroid_y & img_sub_x$Centroid_y<distance*y,]
      if (nrow(img_sub_x_y)>min_cell) {
        cell_freq=prop.table(table(img_sub_x_y$celltype))
        windows_cell_freq[windows_name,celltype[1]]=cell_freq[celltype[1]]
        windows_cell_freq[windows_name,celltype[2]]=cell_freq[celltype[2]]
        windows_cell_freq[windows_name,celltype[3]]=cell_freq[celltype[3]]
        windows_cell_freq[windows_name,celltype[4]]=cell_freq[celltype[4]]
        windows_cell_freq[windows_name,celltype[5]]=cell_freq[celltype[5]]
       # windows_cell_freq[windows_name,celltype[6]]=cell_freq[celltype[6]]
        #return(windows_cell_freq)
      }
    }
    
  } 
  
  windows_cell_freq_all=rbind(windows_cell_freq_all,windows_cell_freq)
  
}
windows_cell_freq_all[is.na(windows_cell_freq_all)]=0
windows_cell_freq_all[windows_cell_freq_all>0.2]<-1
windows_cell_freq_all[windows_cell_freq_all<0.2]<-0
windows_cell_freq_all[windows_cell_freq_all==0.2]<-0
pheatmap(windows_cell_freq_all,show_rownames = F) 
kmeans.result <- kmeans(windows_cell_freq_all, 7)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #合并
kmeans_df_s=arrange(kmeans_df,kmeans_class) #排序
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:7排序
head(kmeans_df_s)
library(ComplexHeatmap)
Heatmap(windows_cell_freq_all[rownames(kmeans_df_s),],cluster_rows = F,show_row_names = F,
        left_annotation= rowAnnotation(df = kmeans_df_s))

write.csv(kmeans_df_s,"kmeans_df_s.csv")      
data_fin_for_cell=data_fin
distance=distance
min_cell=min_cell
celltype=sort(unique(data_fin$celltype))
all_cell_windows=data.frame()
for (i in sam) {
  data_fin_sub=data_fin_for_cell[data_fin_for_cell$study==i,]
  X_max=max(data_fin_sub$Centroid_X)
  Y_max=max(data_fin_sub$Centroid_y)
  x_tile=floor(X_max/distance)
  y_tile=floor(Y_max/distance)
  img_sub_each_sam=data.frame()
  for (x in 1:x_tile) {
    img_sub_x=data_fin_sub[distance*(x-1)<data_fin_sub$Centroid_X & data_fin_sub$Centroid_X<distance*x,]
    for (y in 1:y_tile) {
      windows_name=paste0(i,"_",x,"_",y)
      img_sub_x_y=img_sub_x[distance*(y-1)<img_sub_x$Centroid_y & img_sub_x$Centroid_y<distance*y,]
      if (nrow(img_sub_x_y)>min_cell) {
        img_sub_x_y[,"windows"]=windows_name}
      if (nrow(img_sub_x_y)==min_cell) {
        img_sub_x_y[,"windows"]="few_cell"}
      if (nrow(img_sub_x_y)<min_cell & nrow(img_sub_x_y)>0) {
        img_sub_x_y[,"windows"]="few_cell"}
      if ( nrow(img_sub_x_y)>0) {
        img_sub_each_sam=rbind(img_sub_each_sam,img_sub_x_y)
      }
      
    }
    
  } 
  
  all_cell_windows=rbind(all_cell_windows,img_sub_each_sam)
  
}

input=all_cell_windows

input_few_cell=input[input$windows=="few_cell",]
input_few_cell$windows_type="few_cell"
input_cells=input[which(input$windows!="few_cell"),]
input_cells_imm=input_cells[input_cells$windows%in%rownames(kmeans_df_s),]
input_cells_imm$windows_type=kmeans_df_s[input_cells_imm$windows,"kmeans_class"]
input_cells_tumor=input[-which(input_cells$windows%in%rownames(kmeans_df_s)),]
input_cells_tumor$windows_type="tumor"
input_windows_type=rbind(input_few_cell,input_cells_imm,input_cells_tumor)
data_fin_sub=input_windows_type[input_windows_type$study=="Untitled",]
data_fin_sub=data_fin_sub[order(data_fin_sub$windows_type),]

library(tidyr)
library(data.table)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
require(ggpubr)
library(Seurat)
library(NMF)
library(ggradar)
library(corrr)
library(patchwork)
windows_cell_freq_all$windows=rownames(windows_cell_freq_all)
df_long <- melt(windows_cell_freq_all, id.vars = "windows")
colnames(df_long)<-c("sam_x_y","fincell","Freq")
nmf_input_long = df_long %>% 
  mutate(pt_stg_2 = sam_x_y) %>% separate(pt_stg_2, into = c('sam','x',"y"), sep = '_')# %>% 

nmf_input=reshape2::dcast(nmf_input_long,sam_x_y~fincell,value.var = "Freq")
rownames(nmf_input)=nmf_input[,1]
nmf_input_fin=nmf_input[2:ncol(nmf_input)]
nmf_input_fin[is.na(nmf_input_fin)]<-0
library(ggraph)
library(tidyverse)
library(ggraph)
library(Hmisc)
library(igraph)
group=celltype
nmf_input_fin2=nmf_input_fin
nmf_input_fin2$sam=rownames(nmf_input_fin2)
df_gene=melt(nmf_input_fin2)
colnames(df_gene)<-c("sam","cluster","value")
cluster=names(table(df_gene$cluster))
df_gene<-df_gene%>%group_by(cluster)%>%mutate(median=median(value))
df_gene$ratio<-"low"
df_gene[df_gene$value>df_gene$median,]$ratio<-"high"
df_gene=df_gene[df_gene$cluster%in%group,]
df_ja=matrix(nrow = length(group)-1,ncol = length(group))
for (i in 1:(length(group)-1)) {
  cell_i=group[i]
  #ja=c()
  for (j in (i+1):length(group)) {
    cell_j=group[j]
    a=df_gene[df_gene$cluster==cell_i,]
    a=a[a$ratio=="high",]
    a=a$sam
    b=df_gene[df_gene$cluster==cell_j,]
    b=b[b$ratio=="high",]
    b=b$sam
    jaccard=length(intersect(a,b))/length(union(a,b))
    
    #ja=c(ja,jaccard)
    df_ja[i,j]=jaccard
  }
  
}
df_ja=df_ja[,-1]
rownames(df_ja)=paste0(group[1:(length(group)-1)])
colnames(df_ja)=paste0(group[2:(length(group))])
edge_list <- df_ja %>%
  as_tibble(rownames = "from") %>%
  pivot_longer(cols = -from, names_to = "to", values_to = "weight") %>%
  filter(weight != 0, from != to) # 过滤掉权重为0的边和自环
edge_list=distinct(edge_list)
node=as.data.frame(group)
colnames(node)<-"name"
df_igraph <- graph_from_data_frame(edge_list, directed = FALSE)
df.weight <- E(df_igraph)$weight
edge_attributes <- tibble(weight = df.weight) %>%
  mutate(
    color = case_when(
      weight > 0 ~ "#E6956F",
      weight < 0 ~ "#788FCE",
      TRUE ~ "gray"),width = abs(weight) * 2.5)
df_ja[is.na(df_ja)]<-0
node_sizes_1 <- df_ja %>%
  rowSums() 
node_sizes_2 <- df_ja %>%
  colSums() 
node_sizes=c(node_sizes_1[1],node_sizes_1[2:length(node_sizes_1)]+node_sizes_2[1:(length(node_sizes_2)-1)],node_sizes_2[length(node_sizes_2)])
edge_list=cbind(edge_list,edge_attributes[,-1])
node$size=log10(abs(node_sizes)) * 2.5
df_igraph <- graph_from_data_frame(edge_list,node, directed = FALSE)
ggraph(df_igraph, layout = "linear",circular = TRUE)+
  geom_edge_fan(aes(color=color,edge_width=weight,edge_alpha = weight),show.legend = T)+
  geom_node_point(aes(fill=size,color=size),size=20,shape=21)+
  geom_node_text(aes(label=name, size = weight),
                 angle=0,hjust=0.1,size=3) +
  scale_edge_width_continuous(range = c(1,3),
                              limits = c(0.01,0.8),
                              #breaks = c(0.1),
                              guide = guide_legend(title = "Jaccard index",
                                                   order=1))+
  scale_alpha_continuous(range = c(0.4,0.7))+
  scale_color_gradientn( colours =colorRampPalette(c("#FFED99","#FF7600"))(10))+
  scale_fill_gradientn(#values = seq(0,1,0.2),
     colours =colorRampPalette(c("#FFED99","#FF7600"))(10),
    guide = guide_colorbar(ticks.colour = "black",
                           frame.colour = "black"))+
  scale_edge_colour_manual(values="grey")+
  theme_graph()+
  theme(legend.position = 'none')+
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

library(spatstat)
library(ggplot2)
celltype
cell1="A"
cell2="B"
for (i in unique(data_fin$study)) {
  
  data_fin_for_g_cross=data_fin[data_fin$study==i,]
  data_fin_for_g_cross=data_fin_for_g_cross[data_fin_for_g_cross$celltype%in%c(cell1,cell2),]
  data_fin_for_g_cross$Centroid_X=data_fin_for_g_cross$Centroid_X-min(data_fin_for_g_cross$Centroid_X)
  data_fin_for_g_cross$Centroid_y=data_fin_for_g_cross$Centroid_y-min(data_fin_for_g_cross$Centroid_y)
  window=owin(range(data_fin_for_g_cross$Centroid_X),
              range(data_fin_for_g_cross$Centroid_y))
  spab <- ppp(
    data_fin_for_g_cross$Centroid_X,
    data_fin_for_g_cross$Centroid_y,
    window = window,
    marks = as.factor(data_fin_for_g_cross$celltype)
  )
   df=spab
  G01 <- Gcross(df,r = 0:500)
  save(G01,file = paste0("Gcross",i,"_",cell1,"_",cell2,".RData"))
  pdf(paste0(i,"_",cell1,"_",cell2,"_.pdf"),20,20)
  plot(G01)
  dev.off()
  
}