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
fib.cell=rownames(meta_all)
sub_CT<- subset(data.combined, name %in% fib.cell)
data.combined=sub_CT
data.combined$fincell=Idents(data.combined)
meta_all=meta_all[rownames(data.combined@meta.data),]
meta_all2=meta_all[order(meta_all$fincell),]
meta_all2=meta_all2[order(meta_all2$celltype_main),]
data.combined@meta.data=meta_all
data.combined$fincell=factor(data.combined$fincell,levels = unique(meta_all2$fincell))
Idents(data.combined)<-"fincell"

pdf("Overview_Cluster.pdf",11,4.6)
DimPlot(data.combined, reduction = 'umap', label=F,raster=FALSE,cols =col ) 
dev.off()

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

ranks <- 3:8  
nrun <- 20   
meta = data.combined@meta.data
meta$Patient_origin=paste0(meta$stim,"-",meta$group)
sub_freq_df=data.frame()

for (i in unique(meta$Patient_origin)) {
  meta_sub=meta[meta$Patient_origin==i,]
  tmp=prop.table(table(meta_sub$Patient_origin,meta_sub$fincell),margin = 1)
  sub_freq_df=rbind(sub_freq_df,tmp)
}

colnames(sub_freq_df)<-c("Patient_origin","fincell","Freq")

nmf_input_long = sub_freq_df %>% 
  mutate(pt_stg_2 = Patient_time_origin) %>% separate(pt_stg_2, into = c('Patient','time_origin'), sep = '-')  
nmf_input=reshape2::dcast(nmf_input_long,Patient_time_origin~fincell,value.var = "Freq")
rownames(nmf_input)=nmf_input[,1]
nmf_input_fin=nmf_input[2:ncol(nmf_input)]
nmf_input_fin[is.na(nmf_input_fin)]<-0
nmf_input_fin<-nmf_input_fin[,colSums(nmf_input_fin)!=0]
write.csv(nmf_input_fin,"nmf_input_fin.csv")
cophenetic_coeffs <- numeric(length(ranks))

for (i in seq_along(ranks)) {
  rank <- ranks[i]
  nmf_results <- nmf(nmf_input_fin, rank = rank, nrun = nrun, .options = "v", seed = 123456)
  cophenetic_coeffs[i] <- cophcor(nmf_results)
}
cophenetic_df <- data.frame(FactorizationRank = ranks, CopheneticCoefficient = cophenetic_coeffs)
ggplot(cophenetic_df, aes(x = FactorizationRank, y = CopheneticCoefficient)) +
  geom_point(color = "purple", size = 3) +
  geom_line(color = "purple") +
  labs(title = "Cophenetic correlation survey", x = "Factorization rank", y = "Coefficient") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16))

nrank=4
mtx = nmf_input_fin[, 1:ncol(nmf_input_fin)]
nmf.res = nmf(t(nmf_input_fin[, 1:ncol(nmf_input_fin)]), rank=nrank, nrun=20, seed=123456)
w = nmf.res@fit@W
rownames(w) = colnames(mtx)
colnames(w) = paste0('NMF', seq(1,nrank))

w_ = w %>% t %>% scale
p=ComplexHeatmap::Heatmap(w_, width = 11, height = 4, column_km = nrank, 
                          cluster_rows = F, clustering_method_columns = 'single') 
heatmap = draw(p)
h = nmf.res@fit@H
rownames(h) = paste0('NMF', seq(1,nrank))
h_ = h %>% t %>% scale
p=ComplexHeatmap::Heatmap(h_) 
heatmap = draw(p)

#fig1c,figS1e,f
df_igraph <- graph_from_data_frame(edge_list,node, directed = FALSE)
ggraph(df_igraph, layout = "linear",circular = TRUE)+
  geom_edge_fan(aes(color=color,edge_width=weight,edge_alpha = weight),show.legend = T)+
  geom_node_point(aes(fill=size,color=size),size=20,shape=21)+
  geom_node_text(aes(label=name, size = weight),
                  angle=0,hjust=0.1,size=3) +
  scale_edge_width_continuous(range = c(1,3),
                              limits = c(0.01,1),
                   #breaks = c(0.1),
                   guide = guide_legend(title = "Jaccard index",
                                        order=1))+
  scale_alpha_continuous(range = c(0.4,0.7))+
  scale_color_gradientn( colours =colorRampPalette(c("#FFED99","#FF7600"))(10))+
  scale_fill_gradientn(  colours =colorRampPalette(c("#FFED99","#FF7600"))(10),
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"))+
   scale_edge_colour_manual(values="grey")+
  theme_graph()+
  theme(legend.position = 'none')+
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

ggsave('net1_group4.pdf', width = 6, height = 6)

ggraph(df_igraph, layout = "linear",circular = TRUE)+
  geom_edge_fan(aes(color=color,edge_width=weight,edge_alpha = weight),show.legend = T)+
  geom_node_point(aes(fill=name,color=name),size=20,shape=21,alpha=1)+
  geom_node_text(aes(label=name, size = weight),
                 angle=0,hjust=0.1,size=3) +
  scale_edge_width_continuous(range = c(1,3),
                              limits = c(0.01,1),
                              #breaks = c(0.1),
                              guide = guide_legend(title = "Jaccard index",
                                                   order=1))+
  scale_alpha_continuous(range = c(0.4,0.7))+
  scale_color_manual(values = cols_for_net)+
  scale_fill_manual(values = cols_for_net)+
    scale_edge_colour_manual(values="grey")+
  theme_graph()+
  theme(legend.position = 'none')+
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))
ggsave('net2_group4.pdf', width = 6, height =6)

#fig2d
library(scales)
library(ggradar)
library(tibble)
group="A"
radar_input <- apply(as.matrix(t(h)), 1, function(x) (x-min(x))/(max(x)-min(x))) %>% 
  data.frame(check.names = F)
radar_input=as.data.frame(t(radar_input))
radar_input_fin = radar_input %>% 
  rownames_to_column(var="id")%>% 
  separate(id, into = c('Patient','time_origin'), sep = '-')# %>% 
radar_input$group=radar_input_fin$time_origin
radar_input_sub=radar_input[radar_input$group==group,]
radar_input_sub=radar_input_sub[,1:nrank]
radar_input_sub <- radar_input_sub %>% 
  as_tibble(rownames = "group")
ggradar(
  radar_input_sub,
  values.radar = c("0", "0.5", "1"),
    group.line.width = 1, 
  group.point.size = 0, 
   group.colours = c(rep("#FB7D80",length(unique(radar_input_sub$group)))), 
  grid.line.width = 0.8, 
  background.circle.colour = "white",
  gridline.mid.colour = "grey", 
  legend.position = "bottom" 
)

ggsave(paste0('radar_',group,".pdf"), width = 5, height = 5)

#fig2e
box_input=as.data.frame(t(h))
box_input$id=rownames(box_input)
box_input = box_input %>% 
  separate(id, into = c('Patient','group'), sep = '-') 
box_input=box_input
box_input$group=factor(box_input$group,levels = sort(unique(box_input$group)))
box_input$id=box_input$Patient
write.csv(box_input,"box_input.csv")
box_input_long=melt(box_input)
c <- ggplot(box_input_long,
            aes(x=group, y=value, 
               
                color = group
            )) + 
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", 
               outlier.size = 0.65) +
  facet_wrap(~ variable, scales = 'free_x', nrow = 2) +
  geom_point()+
  scale_color_manual(values = sample_color) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")
p <- c + stat_compare_means(label = "p.signif",comparisons = my_comparisons)
ggsave('box1.pdf',plot = p, width = 6, height = 8)


#fig2f

library(vistributions)
library(copykat)
library(infercnv)
library(ggplot2)
library(ggpubr)
library(ggsci)

dir_main_all="./"
for (group_type in unique(data.combined$group)) {
  
  dir_main_set=paste0(dir_main_all,group_type)
  setwd(dir_main_set)
   infercnv_obj = readRDS("./run.final.infercnv_obj")
  expr <- infercnv_obj@expr.data
  normal_loc <- infercnv_obj@reference_grouped_cell_indices
  normal_loc <- normal_loc$cCON
  test_loc <- infercnv_obj@observation_grouped_cell_indices
  test_loc<-unlist(test_loc)
  pvalue_signif_cutoff=0.01
  expr=expr-1  
  dat_fin=data.frame()
  
  for (i in 1:nrow(expr_tumor)) {
    bin_normal=as.numeric(expr_normal[i,8:ncol(expr_normal)])
    bin_tumor=as.numeric(expr_tumor[i,8:ncol(expr_tumor)])
    
    normal_mean = mean(bin_normal)
    normal_sd = sd(bin_normal)
    
    cut_top=qnorm((1-pvalue_signif_cutoff),mean=normal_mean, sd=normal_sd)
    cut_floor=qnorm((pvalue_signif_cutoff),mean=normal_mean, sd=normal_sd)
    Amp_pos=bin_tumor>cut_top
    Del_pos=bin_tumor<cut_floor

    
  }
  
  colnames(dat_fin)=colnames(expr_tumor[i,8:ncol(expr_tumor)])
  
  Amp_freq= apply(dat_fin, 1, function(x) sum(x=="Amp"))/ncol(dat_fin)
  Del_freq= -apply(dat_fin, 1, function(x)sum(x=="Del"))/ncol(dat_fin)
   chrom=paste0("Chr",data_all$chromosome_name)
  data_all$chrom=chrom
  data_all$chrom=factor(data_all$chrom,levels = paste0("Chr",1:22))
  
  
  for (chorm in unique(data_all$chrom)) {
    data_all[data_all$chrom==chorm,]$abspos=data_all[data_all$chrom==chorm,]$abspos-min(data_all[data_all$chrom==chorm,]$abspos)
  }

  p<-ggplot(data_all) +     
    geom_point(mapping=aes(x=abspos, y=Amp_freq), colour="red", size=0.1) + 
    geom_line(mapping=aes(x=abspos, y=Amp_freq), colour="red") +
    
    geom_point(mapping=aes(x=abspos, y=Del_freq), colour="blue", size=0.1) + 
    geom_line(mapping=aes(x=abspos, y=Del_freq), colour="blue") +     
    geom_hline(yintercept=0, colour="black") +    
    geom_vline(data=data.frame(f=2, x=1), mapping=aes(xintercept=x),
               xintercept=0, 
               linetype=2, size=0.5, col="black") + 
    
    xlab("Chromosome / Position") + #ylab(ylabel) +   
    theme_bw() +    
    theme(panel.spacing = unit(0, "mm"),          
          panel.border = element_blank(),           
          panel.background = element_blank(),           
          panel.grid = element_blank(),          
          strip.background = element_blank(),          
          strip.text.x = element_text(colour="black",size=8,face="plain", angle=90),          
          axis.text.x=element_blank(),          
          axis.text.y = element_text(colour="black",size=8,face="plain"),          
          axis.ticks.x = element_blank(),          
          axis.line.x=element_blank(),          
          axis.title = element_text(colour="black",size=8,face="plain"),          
          plot.title = element_text(colour="black",size=8,face="plain")) #,hjust = 0.5
  
  
  ggsave(paste0("cnv_type1_",group_type,".pdf"), plot = p, width = 12, height = 2.5)
  
}






