# Dimensionality reduction using UMAP and concomitant FlowSOM clustering of mass cytometry data presented in Ingelfinger, Gerdes et al. "Twin study reveals non-heritable immune perturbations in multiple sclerosis"
# This code is losely based on Hartmann, F. J. et al. High-dimensional single-cell analysis reveals the immune signature of narcolepsy. J. Exp. Med. (2016) doi:10.1084/jem.20160897.


# load library
library(RColorBrewer)
library(ggplot2)
library(umap)
library(ggthemes)
library(reshape2)
library(gplots)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(FlowSOM)
library(flowCore)
library(ncdfFlow)

# read data
load(file="/Immunology_shares/Becherlab/People/Ingelfinger/DATA/Multiple Sclerosis/Twin Study/Pooled runs/Surface/0_transformation/2019_06_05_data_twins_merged.RDa")
setwd("/Immunology_shares/Becherlab/People/Ingelfinger/DATA/Multiple Sclerosis/Twin Study/Pooled runs/Surface/1_tsne_flowSOM_main")

# exclude healthy  PBMCs used for batch controls
data <- data[data[,"sample_id"]!= 67,]
data <- data[data[,"sample_id"]!= 68,]
data <- data[data[,"sample_id"]!= 135,]
data <- data[data[,"sample_id"]!= 136,]


# subsample cohort for dimensionality reduction
data_tsne <- data
n_samples <- unique(data[,"sample_id"])

data_sub <- data.frame(NULL)
n_sub_i <- 1500
for(i in n_samples){
  data_i <- data[data[,"sample_id"]==i,]
  n_i <- nrow(data_i)
  set.seed(123)
  if(n_sub_i > n_i){
    ix_i <- sample(1:n_i, n_i)
  }
  else{
    ix_i <- sample(1:n_i, n_sub_i)
  }
  data_sub_i <- data_i[ix_i,]
  data_sub <- rbind(data_sub, data_sub_i)
}

clustering_cols <- colnames(data_sub)[c(2:37)]
data_dimred <- data_sub[ , clustering_cols]


# UMAP
data_umap <-  umap(data_dimred, random_state=123, verbose =T)
umap <- as.data.frame(data_umap$layout)
colnames(umap) <- c("UMAP1", "UMAP2")

# plot UMAP black
umap_black <- ggplot(umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 0.5) +
  coord_fixed(ratio = 1)
umap_black

# plot UMAP with expression overlayed
data_dimred_exp <- cbind(data_dimred, umap)
data_dimred_melt <- melt(data_dimred_exp, id.vars = c("UMAP1", "UMAP2"))
data_dimred_comb  <- cbind(data_sub[,c("gate_source", "cell_id", "sample_id")], data_dimred_exp)

color_grad_flow2 <- c("#331820", "#4d4f55", "#55626b", "#5a767e", "#628b8a", "#709b91", "#82aa96", "#98b89a", "#b0c6a2", "#c9d3ab", "#e4e0b6", "#feedc3")

umap_overlay <- ggplot(data_dimred_melt, aes(x = UMAP1, y = UMAP2, color = value)) +
  geom_point(size = 0.05) +
  coord_fixed(ratio = 1) +
  scale_colour_gradientn(colours = color_grad_flow2, limits = c(0,1)) +
  facet_wrap(~ variable, ncol = 6) +
  theme_few() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
umap_overlay


# FlowSOM clustering, define markers to use for clustering
data_flow <- data.matrix(data)
clusteringcol_flow <- colnames(data_flow)[c(2, 3, 18, 20, 23, 24, 25, 28, 30, 31, 36, 37)] 

# generate flowframe
ff_new <- flowFrame(exprs = data_flow, desc = list(FIL = 1))

# run FlowSOM (with set.seed for reproducibility)
set.seed(123)
out_fSOM <- FlowSOM::ReadInput(ff_new, transform = F, scale = F, compensate = F)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = clusteringcol_flow)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)
labels <- out_fSOM$map$mapping[,1]

# make heatmap of 100 FlowSOM clusters
heat_mat <- matrix(NA, nrow = 100, ncol = length(clusteringcol_flow))
for(i in 1:100) {
  temp_mat <- data_flow[labels == i, clusteringcol_flow]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(heat_mat) <- 1:100
colnames(heat_mat) <- clusteringcol_flow

# plot heatmap of expression for 100 clusters
breaks <- seq(0, 1, by = 0.111111111)
white.black <- colorRampPalette(c("white", "black"))(n = 9)

heatmap.2(heat_mat,
          scale = "none",
          Colv = F, Rowv = T,
          trace = "none",
          col = white.black,
          breaks = breaks)

# set the max number of clusters for concomitant metaclustering
max  <- 20

# initialize matrix
gate_source <- as.factor(data[,"gate_source"])
meta_results <- data.frame(gate_source)

# metaclustering
for (i in 3:max) {
  set.seed(123)
  out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = i)
  meta_results <- cbind(meta_results, as.factor(out_meta[labels]))}

meta_results <- meta_results[,2:ncol(meta_results)]
colnames(meta_results) <- paste("k.", 3:max, sep = "")
meta_results <- cbind(meta_results, data[,"cell_id"])
colnames(meta_results)[19] <- "cell_id"

# prepare the metaclustering data
data_meta <- data_dimred_comb[,c("cell_id", "UMAP1", "UMAP2")]
meta.melt <- melt(meta_results, id.vars = "cell_id", variable.name = "k.value", value.name = "cluster.assigment")
meta.melt$cluster.assigment <- as.factor(meta.melt$cluster.assigment)
joined.meta <- merge(meta.melt, data_meta, by = "cell_id")

# plot umap with metaclusters overlayed
db13 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

umap_meta <- ggplot(joined.meta, aes(x = UMAP1, y = UMAP2, color = cluster.assigment)) +
  geom_point(size = 0.25) +
  coord_fixed(ratio = 1) +
  scale_colour_manual(name = NULL,values = c(db13)) +
  facet_wrap(~ k.value, ncol =6) +
  guides(colour = guide_legend(override.aes = list(size=5), title="Cluster"))
umap_meta


# consensus clustering to determine number of metaclusters
set.seed(123)
results <- ConsensusClusterPlus(t(heat_mat), maxK = 20, reps = 1000, pItem = 0.8,
                                pFeature = 1,  clusterAlg = "hc", verbose = F,
                                distance = "euclidean", seed = 123, plot = "png",
                                writeTable = T)


# perform metaclustering
choosen_k <- 8
set.seed(123)
out_meta <- metaClustering_consensus(out_fSOM$map$codes, k = choosen_k)
pop_labels <- out_meta[labels]

# generate a heatmap showing the expression for each cluster
heat_mat <- matrix(NA, nrow = choosen_k, ncol = length(clusteringcol_flow))
for(i in 1:choosen_k) {
  temp_mat <- data_flow[pop_labels == i, clusteringcol_flow]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(heat_mat) <- paste("cluster", 1:choosen_k, sep = "")
colnames(heat_mat) <- clusteringcol_flow

heatmap.2(heat_mat,
          scale = "none",
          Colv = F, Rowv = T,
          trace = "none",
          col = white.black,
          breaks = breaks,
          mar=c(6,9))

# manually merge and annotate clusters
CD4 <- c(2)
CD8 <- c(4)
non.conv.Tcells <- c(5)
NK.cells <- c(8, 3)
B.cells <- c(1)
Myeloid <- c(7, 6)

# assign a number for each metacluster
labels <- rep(0, length(pop_labels))
labels[pop_labels %in% CD4] <- 1
labels[pop_labels %in% CD8] <- 2
labels[pop_labels %in% non.conv.Tcells] <- 3
labels[pop_labels %in% NK.cells] <- 4
labels[pop_labels %in% B.cells] <- 5
labels[pop_labels %in% Myeloid] <- 6

# heatmap with meta clustering
heat_mat <- matrix(NA, nrow = 6, ncol = length(clusteringcol_flow))
for(i in 1:6) {
  temp_mat <- data[labels == i, clusteringcol_flow]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

clusternames <- c("CD4 T cells", "CD8 T cells", "unconventional T cells", "NK cells", "B cells", "Myeloid cells")
rownames(heat_mat) <- clusternames
colnames(heat_mat) <- clusteringcol_flow

heatmap.2(heat_mat,
          scale = "none",
          Colv = F, Rowv = T,
          trace = "none",
          col = white.black,
          breaks = breaks,
          mar=c(6,18))

# plot resulting clusters on umap
data_clust <- cbind(data, labels)
data_umap_clust <- merge(data_dimred_comb, data_clust[,c("labels", "cell_id")], by = "cell_id")
data_umap_clust$labels <- as.factor(data_umap_clust$labels)


color_qual_flow2 <- c("#557F7A", "#7AA591", "#CBD49C", "#FEEDC3", "#CC242A", "#FFD258")

umap_clust <- ggplot(data_tsne_compl, aes(x = tSNE1, y = tSNE2, color = labels)) +
  geom_point(size = .1) +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = color_qual_flow2, name= "labels", breaks=1:6,
                     labels= clusternames) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title=element_blank(),
        legend.background = element_rect())
umap_clust

# save clustered dataframe
colnames(data_clust)[ncol(data_clust)] <- "main_clusters"
save(data_clust, file="data_main_clust_k6.RDa")
