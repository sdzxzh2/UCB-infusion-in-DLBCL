
library(VISION)
library(msigdbr)

for (i in gmt) {
counts <-as.matrix(data.combined@assays$RNA@counts)
n.umi <- colSums(counts)
scaled_counts <- t(t(counts) / n.umi) * median(n.umi)
meta = data.combined@meta.data
vis <- Vision(scaled_counts,
              signatures = paste0("./",i),
              min_signature_genes = 2,
              meta = meta)

options(mc.cores = 13)
vis <- analyze(vis)
projection <- data.combined@reductions$umap@cell.embeddings
vis <- addProjection(vis, "UMAP", projection)
save(vis, file=paste0("vis_mono_",i,".RData"))

}

library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(SCENIC)

cell.info <- data.combined@meta.data
dge <- as.matrix(data.combined@assays$RNA@counts)

scenicOptions <- initializeScenic(
  org="mgi", 
  dbDir="./cisTarget_databases", # RcisTarget databases location
  dbs="hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",# file name of motif database
  datasetTitle="SCENIC on Mouse Cell Atlas", # choose a name for your analysis
  nCores=1
)


genesKept <- geneFiltering(
  exprMat=dge, 
  scenicOptions=scenicOptions,
  minCountsPerGene = 1,
  minSamples = 20
)

length(genesKept)
dge <- dge[genesKept, ]
dim(dge)

if(!dir.exists("output")) {
  dir.create("output")
}
saveRDS(cell.info, "output/s1_cell.info.rds")
source("utils/add_cellAnnotation.R")
saveLoom <- function(exprMat, output){
  cellInfo2 <- data.frame(
    row.names = rownames(cell.info),
    cellType =cell.info$fincell
  )
  loom <- build_loom(output, dgem=exprMat)
  loom <- add_cellAnnotation(loom, cellInfo2)
  close_loom(loom)
}

saveLoom(dge, "output/s1_avg20_rep1.loom")
loom <- build_loom("output/s1_exprMat.loom", dgem=dge)
loom <- add_cellAnnotation(loom, cell.info)
close_loom(loom)
library(ggplot2)
library(plyr)
library(ggsci)
emb.tsne <- read.table("output/s3_avg20_rep1.tsne.txt", sep = "\t", row.names = 1, header = T)
emb.umap <- read.table("output/s3_avg20_rep1.umap.txt", sep = "\t", row.names = 1, header = T)
colnames(emb.tsne) <- paste0("tSNE_", 1:2)
colnames(emb.umap) <- paste0("UMAP_", 1:2)
cell.info <- readRDS("output/s1_cell.info.rds")
cell.info <- cbind(cell.info, emb.tsne[rownames(cell.info), ])
cell.info <- cbind(cell.info, emb.umap[rownames(cell.info), ])
saveRDS(cell.info, "output/s4_cell.info.rds")
cell.info <- readRDS("output/s4_cell.info.rds")
get_label_pos <- function(data, emb = "tSNE", group.by="fincell") {
  new.data <- data[, c(paste(emb, 1:2, sep = "_"), group.by)]
  colnames(new.data) <- c("x","y","cluster")
  clusters <- names(table(new.data$cluster))
  new.pos <- lapply(clusters, function(i) {
    tmp.data = subset(new.data, cluster == i)
    data.frame(
      x = median(tmp.data$x),
      y = median(tmp.data$y),
      label = i)
  })
  do.call(rbind, new.pos)
}
ggplot(cell.info, aes(tSNE_1, tSNE_2, color=as.character(fincell))) + 
  geom_point(size=.1) + scale_color_manual(values = celltype_col)+
  geom_text(inherit.aes = F, data = get_label_pos(cell.info, emb = "tSNE"), aes(x,y,label=label), size=3) + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black")
  )
ggplot(cell.info, aes(UMAP_1, UMAP_2, color=fincell)) + 
  geom_point(size=.1) + scale_color_manual(values = celltype_col)+
  geom_text(inherit.aes = F, data = get_label_pos(cell.info, emb = "UMAP"), aes(x,y,label=label), size=3) + 
  theme_bw(base_size = 15) + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black")
  )
umapPlot.tissue <- function(data, tissue.use, topn) {
  clusters.use <- table(subset(data, Tissue == tissue.use)$ClusterID)
  if(topn > length(clusters.use)) {
    topn <- length(clusters.use)
  }
  clusters.use <- names(head(sort(clusters.use, decreasing = T), n=topn))
  data$new.cluster <- ifelse(data$Tissue == tissue.use & data$ClusterID %in% clusters.use, data$CellType, "Others")
  data$new.cluster <- factor(data$new.cluster, levels = c(setdiff(names(table(data$new.cluster)), "Others"), "Others"))
  data$pt.size <- ifelse(data$new.cluster == "Others", 0.05, 0.2)
  ggplot(data, aes(UMAP_1, UMAP_2, color=new.cluster)) + 
    geom_point(size=data$pt.size) + 
    scale_color_manual(values = c(pal_d3("category10")(topn), "grey")) + 
    guides(colour = guide_legend(keyheight = .7, 
                                 keywidth = .1, 
                                 override.aes = list(size=3))) + 
    theme_bw(base_size = 12) + 
    ggtitle(tissue.use) + 
    theme(legend.justification = c(0,0), 
          legend.position = c(0,0),
          legend.title = element_blank(),
          legend.key = element_rect(fill = alpha("white", 0)),
          legend.background = element_rect(fill=alpha('white', 0)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          plot.title = element_text(hjust = .5, face = "bold")
    )
}

cell.info$Tissue=cell.info$clone_group1
cell.info$ClusterID=cell.info$fincell
cell.info$CellType=cell.info$fincell
umapPlot.tissue(cell.info, "clone", 20)
umapPlot.cellType <- function(data, cell.type, topn) {
  tissues.use <- table(subset(data, CellType == cell.type)$Tissue)
  if(topn > length(tissues.use)) {
    topn <- length(tissues.use)
  }
  tissues.use <- names(head(sort(tissues.use, decreasing = T), n=topn))
  data$new.tissue <- ifelse(data$CellType == cell.type & data$Tissue %in% tissues.use, data$Tissue, "Others")
  data$new.tissue <- factor(data$new.tissue, levels = c(tissues.use, "Others"))
  data$pt.size <- ifelse(data$new.tissue == "Others", 0.05, 0.2)
  
  ggplot(data, aes(UMAP_1, UMAP_2, color=new.tissue)) + 
    geom_point(size=data$pt.size) + 
    scale_color_manual(values = c(pal_d3("category10")(topn), "grey")) + 
    
    guides(colour = guide_legend(keyheight = .7, 
                                 keywidth = .1,
                                 override.aes = list(size=3))) + 
    theme_bw(base_size = 12) + 
    ggtitle(cell.type) + 
    theme(legend.justification = c(0,0), 
          legend.position = c(0,0),
          legend.title = element_blank(),
          legend.key = element_rect(fill = alpha("white", 0)),
          legend.background = element_rect(fill=alpha('white', 0)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          plot.title = element_text(hjust = .5, face = "bold")
    )
}