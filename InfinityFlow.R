# Preprocessing of InfinityFlow data, clustering and mapping onto mass cytometry data of the MS twin cohort as presented in 
# Ingelfinger, Gerdes et al. "Twin study reveals non-heritable immune perturbations in multiple sclerosis"

# library
library(flowCore)
library(RColorBrewer)
library(gplots) 
library(ggplot2)
library(gdata)
library(gridExtra)
library(grid)
library(ncdfFlow)
library(reshape2)
library(ggthemes)
library(umap)
library(pheatmap)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(FlowSOM)
library(vite)
library(panorama)
library(grappolo)
library(igraph)
library(ggraph)

# load data and concatenate different samples 
setwd("~/becherlab_shares/People/Ingelfinger/DATA/Multiple Sclerosis/Twin Study/infinityflow/Human Blood LegendScreen for Florian/fcs/")
files_names <- list.files(getwd(), pattern='.fcs$', full=FALSE)
fs <- read.FCS(files_names[6], transformation = F, truncate_max_range = F)

# combine data into a matrix and rename
data <- exprs(fs)
names(colnames(data)) <- NULL

colnames(data)[1:18] <- c("BB-CD11b (act)", "FSC-A", "FSC-H", "SSC-A", "SSC-W", "BUV737", "BB-CD123", "BB-HLA-DR",
                          "BB-CD5", "BB-CD163", "BB-CD45", "BB-CD2", "BB-CD1c", "BB-CD45RA", "BB-CD88", "BB-CD16", "BB-CD14", "BB-CD141")

cell_id <- 1:nrow(data)
data <- cbind(cell_id, data)

# # save files
save(data, file="/../../2021-02-09_infinityflow_untransf.RDa")

# convert to matrix
markers <- colnames(data)
untrans <- colnames(data)[c(1, 3, 4, 5, 6, 7)]
trans <- setdiff(markers, untrans)

data_matrix_surf1 <- data.matrix(data)

# transform data
asinh_scale <- 1000

colnames(data_matrix_surf1)

data.trans_surf1 <- asinh(data_matrix_surf1/ asinh_scale)
data.trans_surf1[,"BB-CD123"] <- asinh(1000*sinh(data.trans_surf1[,"BB-CD123"])/2000) 
data.trans_surf1[,"BB-CD5"] <- asinh(1000*sinh(data.trans_surf1[,"BB-CD5"])/2000)
data.trans_surf1[,"BB-CD163"] <- asinh(1000*sinh(data.trans_surf1[,"BB-CD163"])/3000) 
data.trans_surf1[,"BB-CD2"] <- asinh(1000*sinh(data.trans_surf1[,"BB-CD2"])/2000) 
data.trans_surf1[,"BB-CD1c"] <- asinh(1000*sinh(data.trans_surf1[,"BB-CD1c"])/2500) 
data.trans_surf1[,"BB-CD14"] <- asinh(1000*sinh(data.trans_surf1[,"BB-CD14"])/2000)
data.trans_surf1[,"BB-CD141"] <- asinh(1000*sinh(data.trans_surf1[,"BB-CD141"])/2000) 

data.trans_surf1[,untrans] <- data_matrix_surf1[,untrans]

# shift to 0 and normalize everything between 0 and 1
q.vector <- apply(data.trans_surf1, 2, function(x) quantile(x, 0.001, names = F))
data.shift <- data.trans_surf1
data.shift <- sweep(data.shift, 2, q.vector)

data.trans.new_surf1 <- data.shift
per.vector_surf1 <- apply(data.trans.new_surf1, 2, function(x) quantile(x, 0.9999, names = F))
data.trans.new_surf1 <- t(t(data.trans.new_surf1) / as.numeric(per.vector_surf1))

data.trans.new_surf1[,untrans] <- data_matrix_surf1[,untrans]

# make biaxial plots for backbone markers and all infinity markers
data_trans <- data.frame(data.trans.new_surf1)
bb_marker <- colnames(data_trans)[c(2, 8:19)]


db13 <- c("#1B9E77","gold" , "#7570B3", "#E7298A", "#D95F02" ,  "lightpink1", "#666666",
          "#57B0FF", "#66CD00", "#1F78B4", "#B2DF8A", "firebrick1", 
          "seagreen1", "brown" , "#6A3D9A", "mediumorchid", "grey" ,
          "#BC80BD", "#80B1D3")

plot_surf1_list <- list()
for(i in 1:length(bb_marker)){
  plot_surf1 <- ggplot(data_trans, aes_string(x = bb_marker[i], y = 'BB.CD14')) +
    geom_hex(bins = 100) +
    xlim(-0.1, 1.5) +
    ylim(-0.1, 1.5) +
    coord_fixed(ratio = 1) +
    scale_fill_gradientn(colours = db13, trans = "sqrt") +
    theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 3), axis.title = element_text(size = 3))
  plot_surf1_list[[i]] <-  plot_surf1
}

m1 <- arrangeGrob(grobs = plot_surf1_list, ncol = 1, top=textGrob("InfinityFlow", gp=gpar(fontsize=6)))
ggsave(file="marker_comp_total_BB.pdf", m1, width = 1, height = 15, limitsize = F)

plot_surf2_list <- list()
for(i in 1:length(colnames(data_trans))){
  plot_surf1 <- ggplot(f_surf1, aes_string(x = colnames(f_surf1)[i], y = 'BB.CD14')) +
    geom_hex(bins = 100) +
    xlim(-0.1, 1.5) +
    ylim(-0.1, 1.5) +
    coord_fixed(ratio = 1) +
    scale_fill_gradientn(colours = db13, trans = "sqrt") +
    theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 3), axis.title = element_text(size = 3))
  plot_surf2_list[[i]] <-  plot_surf1
}

m2 <- arrangeGrob(grobs = plot_surf2_list, ncol = 1, top=textGrob("InfinityFlow", gp=gpar(fontsize=6)))
ggsave(file="marker_comp_total.pdf", m2, width = 1, height = 350, limitsize = F)

# save transformed files
save(data_trans, file = "2021-02-10_infinityflow_transf.RDa")

# run umap using backbone markers only
clustering_cols <- trans
data_dimred <- data_trans[ , clustering_cols]

# UMAP
data_umap <-  umap(data_dimred, random_state=123, verbose =T)
umap <- as.data.frame(data_umap$layout)
colnames(umap) <- c("UMAP1", "UMAP2")

# plot UMAP black
umap_black <- ggplot(umap, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = 0.5) +
  coord_fixed(ratio = 1)
umap_black

# prepare the expression data
data_plot_u <- cbind(data_dimred, umap)
data_melt_u <- melt(data_plot_u, id.vars = c("UMAP1", "UMAP2"))

# plot UMAP with expression overlayed
color_grad_flow2 <- c("#331820", "#4d4f55", "#55626b", "#5a767e", "#628b8a", "#709b91", "#82aa96", "#98b89a", "#b0c6a2", "#c9d3ab", "#e4e0b6", "#feedc3")

umap_expr <- ggplot(subset(data_melt_u, variable %in% bb_marker), aes(x = UMAP1, y = UMAP2, color = value)) +
  geom_point(size = 0.05) +
  coord_fixed(ratio = 1) +
  scale_colour_gradientn(colours = color_grad_flow2, limits = c(0,1)) +
  facet_wrap(~ variable, ncol = 4) +
  theme_few() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
umap_expr

# save data
data_compl <- cbind(data_trans, umap)
save(data_compl, file ="2021-02-10_infinityflow_transf_umap.RDa")

# FlowSOM using backbone markers
colnames(data_compl)[352:353] <- c("sUMAP1", "sUMAP2")
data_flow <- data_compl
data_flow <- data.matrix(data_flow)
clusteringcol_flow <- bb_marker[c(1:5, 7, 8, 10, 11, 12, 13)]
ff_new <- flowFrame(exprs = data_flow, desc = list(FIL = 1))

# run FlowSOM (with set.seed for reproducibility)
set.seed(123)
out_fSOM <- FlowSOM::ReadInput(ff_new, transform = F, scale = T, compensate = F)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = clusteringcol_flow)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)
labels <- out_fSOM$map$mapping[,1]

# make heatmap
heat_mat <- matrix(NA, nrow = 100, ncol = length(clusteringcol_flow))
for(i in 1:100) {
  temp_mat <- data_flow[labels == i, clusteringcol_flow]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(heat_mat) <- 1:100
colnames(heat_mat) <- clusteringcol_flow

breaks <- seq(0, 1, by = 0.111111111)
white.black <- colorRampPalette(c("white", "black"))(n = 9)

heatmap.2(na.omit(heat_mat),
          scale = "none",
          Colv = F, Rowv = T,
          trace = "none",
          col = white.black,
          breaks = breaks)

# set the maximum number of metaclusters you want
max  <- 20

# meta clustering
data_compl[,"gate_source"] <- 1
gate_source <- as.factor(data_compl[,"gate_source"])
meta_results <- data.frame(gate_source)

for (i in 3:max) {
  set.seed(123)
  out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = i)
  meta_results <- cbind(meta_results, as.factor(out_meta[labels]))}

meta_results <- meta_results[,2:ncol(meta_results)]
colnames(meta_results) <- paste("k.", 3:max, sep = "")
meta_results <- cbind(meta_results, data_compl[,"cell_id"])
colnames(meta_results)[19] <- "cell_id"
data_meta <- data_compl[,c("cell_id", "UMAP1", "UMAP2")]

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

# consensus clustering to determine how many clusters are there
set.seed(123)
results <- ConsensusClusterPlus(t(na.omit(heat_mat)), maxK = 20, reps = 1000, pItem = 0.8,
                                pFeature = 1,  clusterAlg = "hc", verbose = F,
                                distance = "euclidean", seed = 123, plot = "png",
                                writeTable = T)


# perform actual metaclustering
choosen_k <- 15
set.seed(123)
out_meta <- metaClustering_consensus(out_fSOM$map$codes, k = choosen_k)
pop_labels <- out_meta[labels]

# heatmap of metaclusters
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

# manually merge and annotate metaclusters
ncM <- c(10)
iM <- c(7)
cM <- c(1)
pDCs <- c(14, 13, 15)
cDC1 <- c(12)
cDC2 <- c(4, 3, 6, 9, 8, 5, 2)
unassigned <- c(11)

# give them the right metacluster number
labels <- rep(0, length(pop_labels))
labels[pop_labels %in% ncM] <- 1
labels[pop_labels %in% iM] <- 2
labels[pop_labels %in% cM] <- 3
labels[pop_labels %in% pDCs] <- 4
labels[pop_labels %in% cDC1] <- 5
labels[pop_labels %in% cDC2] <- 6
labels[pop_labels %in% unassigned] <- 7

# heatmap with manually annnotated clusters
heat_mat <- matrix(NA, nrow = 7, ncol = length(clusteringcol_flow))
for(i in 1:7) {
  temp_mat <- data_flow[labels == i, clusteringcol_flow]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

clusternames <- c("non-classical monocytes", "intermediate monocytes", "classical monocytes", "pDCs", "cDC1s", "cDC2s", "unassigned")
rownames(heat_mat) <- clusternames
colnames(heat_mat) <- clusteringcol_flow

pheatmap(mat = heat_mat, 
         scale= "none",
         color = white.black,
         cluster_rows=F,
         cluster_cols=F,
         cellwidth = 10, cellheight = 10,
         border_color = NA,
         filename= "main_clusters_heatmap.pdf"
)


# heatmap for all imputed markers
all_marker <- colnames(data_flow)[c(2, 8:351)]
data_flow <- cbind(data_flow, labels)
data_flow_red <- subset(data_flow, labels != 7)

heat_mat <- matrix(NA, nrow = 6, ncol = length(all_marker))
for(i in 1:6) {
  temp_mat <- data_flow_red[data_flow_red[,"labels"] == i, all_marker]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(heat_mat) <- clusternames[1:6]
colnames(heat_mat) <- all_marker

# make a cluster heatmap
pheatmap(mat = heat_mat, 
         scale= "none",
         color = white.black,
         # breaks = breaks,
         cluster_rows=F,
         # cluster_cols=F,
         # cellwidth = 10, cellheight = 10,
         border_color = NA,
         filename= "main_clusters_heatmap_all_markers.pdf"
)

# top 100 variable markers
var.markers <- data.frame(markers = colnames(heat_mat), variance = colVars(heat_mat))
var.markers <- var.markers[order(var.markers$variance, decreasing = T),]
rownames(var.markers) <- NULL

pheatmap(mat = heat_mat[,var.markers[1:100,1]], 
         scale= "none",
         color = white.black,
         # breaks = breaks,
         cluster_rows=F,
         # cluster_cols=F,
         cellwidth = 10, cellheight = 10,
         border_color = NA,
         filename= "main_clusters_heatmap_top_markers.pdf"
)

# plot results on umap
data_clust <- cbind(data_flow, labels)

data_clust.df <- data.frame(data_clust)
data_clust.df <- subset(data_clust.df, labels != 7)
data_clust.df$labels <- factor(data_clust.df$labels, levels = 1:6, labels = clusternames[1:6])

color_qual_flow2 <- c("classical monocytes" = "#557F7A", "intermediate monocytes" = "#7AA591", "non-classical monocytes" = "#CBD49C", "cDC2s" = "#FFD258",  "pDCs" = "#FEEDC3", "cDC1s" = "grey")

umap_clust <- ggplot(data_clust.df, aes(x = UMAP1, y = UMAP2, color = labels)) +
  geom_point(size = .0001) +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = color_qual_flow2, name= "labels", breaks=1:7,
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
ggsave(umap_clust, filename="umap_main_clusters_k7.png")

# save clustered data in new FCS file
data_clustered <- data_clust.df
colnames(data_clustered)[ncol(data_clustered)] <- "main_clusters"
save(data_clustered, file="2021-02-10_data_main_clust_k7.RDa")


# mapping onto mass cytometry data of the MS twin cohort
# transfer cluster labels to untransformed subset
load(file = "2021-02-09_infinityflow_untransf.RDa")
data_inf_untrans_clust <- merge(data, data_clustered[,c("cell_id", "main_clusters")], by= "cell_id")

# export landmark nodes
n_landmark <- as.numeric(unique(data_inf_untrans_clust[,"main_clusters"]))
clusternames <- c("non-classical monocytes", "intermediate monocytes", "classical monocytes", "pDCs", "cDC1s", "cDC2s", "unassigned")

colnames(data_inf_untrans_clust)[c(2, 8:19, 34, 204, 70, 196, 230, 191, 147, 194, 30, 195, 335, 86)] <- c("CD11b.b", "CD123.b", "HLA-DR.b", "CD5.b", "CD163.b", "CD45.b", "CD2.b", "CD1c.b", "CD45RA.b", "CD88.b", "CD16.b", "CD14.b", "CD141.b", "CD11b.sim",
                                                                                                                 "CD207", "CD45R0", "CCR7", "CD268", "CXCR3", "CD127", "CCR5", "CD8", "CCR6", "TCRgd", 
                                                                                                                 "CD56")
data_inf_untrans_clust.m <- data.matrix(data_inf_untrans_clust)

for (i in n_landmark) {
  sample_i <- data_inf_untrans_clust.m[data_inf_untrans_clust.m[,"main_clusters"]==i,]
  ff_i <- flowFrame(sample_i)
  print(paste("sample_infinity_", clusternames[i], ".fcs", sep=""))
  write.FCS(ff_i, filename =  paste("sample_surface_", clusternames[i], ".fcs", sep=""))
}

# generate fcs file for infinity flow data
colnames(data_clustered)[c(2, 8:19, 34, 204, 70, 196, 230, 191, 147, 194, 30, 195, 335, 86)] <- c("CD11b.b", "CD123.b", "HLA-DR.b", "CD5.b", "CD163.b", "CD45.b", "CD2.b", "CD1c.b", "CD45RA.b", "CD88.b", "CD16.b", "CD14.b", "CD141.b", "CD11b.sim",
                                                                                                          "CD207", "CD45R0", "CCR7", "CD268", "CXCR3", "CD127", "CCR5", "CD8", "CCR6", "TCRgd", "CD56")
                      
data_clustered_filt <- subset(data_clustered, main_clusters != 7)
                                                                                                                                                                                                                                                                         
ff_i <- flowFrame(data.matrix(data_clustered_filt))
write.FCS(ff_i, filename =  ""../"blood_infinity.fcs")

# generate MS twin file
# load diffcyt clusters of myeloid cells and bring into format
load(file="/Immunology_shares/Becherlab/People/Ingelfinger/DATA/Multiple Sclerosis/Twin Study/Pooled runs/Surface/11_diffcyt/2019-12-16_DS_clustered_full.RDa")
data_MS_clust <- subset(data_incl_clust, main_clusters == 6 | cluster_id == 77)

n_samples <- unique(data_MS_clust[,"sample_id"])

data_sub <- data.frame(NULL)
n_sub_i <- 5000
for(i in n_samples){
  data_i <- data_MS_clust[data_MS_clust[,"sample_id"]==i,]
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

ff_i <- flowFrame(data.matrix(data_sub))
write.FCS(ff_i, filename =  "MS_twins.fcs")

# Transfer of diffcyt clusters
colnames(data_sub)[c(3, 13, 19, 20, 26, 38)] <- c("CD11b.b", "CD45RA.b", "CD14.b", "CD5.b", "CD16.b", "HLA-DR.b")
col.names <- colnames(data_sub)[3:38]

markers_heat <- col.names
heat_mat <- matrix(NA, nrow = length(unique(data_sub$cluster_id)), ncol = length(markers_heat))
for(i in 1:length(unique(data_sub$cluster_id))) {
  temp_mat <- data_sub[data_sub[,"cluster_id"] == unique(data_sub$cluster_id)[i], markers_heat]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

colnames(heat_mat) <- col.names
rownames(heat_mat) <- unique(data_sub$cluster_id)

pop_size <- data.frame(table(data_sub$cluster_id))
colnames(pop_size) <- c("cellType", "popsize")

data_diff_exp <- cbind(unique(data_sub$cluster_id), data.frame(heat_mat), rep("data_diffcyt_pooled.fcs", nrow(heat_mat)))
colnames(data_diff_exp) <- c("cellType", col.names, "sample")
data_diff_exp <- merge(data_diff_exp, pop_size)


# remove clusters with less than 100 cells
data_diff_exp_filt <- subset(data_diff_exp, popsize > 100)

# export clustered file
write.table(data_diff_exp_filt, file="data_diffcyt_pooled.fcs.clustered.txt", sep="\t", row.names=F, col.names = T, quote = F)

# load landmarks
landmarks.data <- load_landmarks_from_dir("gated/", asinh.cofactor = 1500, transform.data = T)

col.names.com <- colnames(data_diff_exp)[colnames(data_diff_exp) %in% colnames(landmarks.data$landmarks.data)][2:33]

# perform clara clustering on infinity flow data
clust.col <- colnames(data_clustered)[c(2, 8:19)]
files.list <- list.files(path = getwd(), pattern = ".fcs$")
cluster_fcs_files(files.list[1], num.cores = 4, col.names = clust.col, num.clusters = 100,
                  asinh.cofactor = NULL)

# compute scaffolds
input.files <- list.files(pattern = ".txt")
run_scaffold_analysis(input.files, ref.file = input.files[1], 
                      landmarks.data = landmarks.data, col.names = col.names.com, process.clusters.data = F)

# generate networks
g_MS_diff <- read_graph(file="scaffold_result/data_diffcyt_pooled.fcs.graphml", format= "graphml")
g_inf <- read_graph(file="scaffold_result/blood_infinity.fcs.graphml", format= "graphml")

l_g_MS_diff <- create_layout(g_MS_diff, layout="fr")
l_g_MS_diff$x <- V(g_MS_diff)$x
l_g_MS_diff$y <- V(g_MS_diff)$y

l_g_inf <- create_layout(g_inf, layout="fr")
l_g_inf$x <- V(g_inf)$x
l_g_inf$y <- V(g_inf)$y

ann_col <- c("landmark" = "#557F7A", "diffcyt"= "#FEEDC3", "selected" = "#CC242A")
sig_clust <- paste("c", c(15, 6, 37, 23, 13, 34, 63, 77), sep="")

# plot networks
graph_diffcyt_overview <- ggraph(l_g_MS_diff) +
  geom_edge_link(alpha=.15) + 
  geom_node_point(data = subset(l_g_MS_diff, type== "cluster"), aes(size= popsize, fill="diffcyt"), shape=21)+
  geom_node_point(data = subset(l_g_MS_diff, type== "cluster" & Label %in% sig_clust), aes(size= popsize, fill = "selected"), shape=21)+
  geom_node_point(data= subset(l_g_MS_diff, type== "landmark"), aes(fill = "landmark"), shape=21, size=4) +
  scale_fill_manual(values=ann_col)+
  scale_size_continuous(range = c(4, 14)) +
  geom_node_text(data= subset(l_g_MS_diff, type == "landmark"), aes(label = name), size = 6, repel=T)+
  theme_graph(base_family = 'Helvetica') + theme(legend.position="none")
graph_diffcyt_overview
ggsave(filename="Scaffold_MS_twins_diff.pdf", graph_diffcyt_overview, width=12, height= 8, units="in", useDingbats=FALSE)

graph_inf <- ggraph(l_g_inf) +
  geom_edge_link(alpha=.15) + 
  geom_node_point(data = subset(l_g_inf, type== "cluster"), aes(size= popsize, fill="diffcyt"), shape=21)+
  geom_node_point(data= subset(l_g_inf, type== "landmark"), aes(fill = "landmark"), shape=21, size=4) +
  scale_fill_manual(values=ann_col)+
  scale_size_continuous(range = c(4, 14)) +
  geom_node_text(data= subset(l_g_inf, type == "landmark"), aes(label = name), size = 6, repel=T)+
  theme_graph(base_family = 'Helvetica') + theme(legend.position="none")
graph_inf
ggsave(filename="Scaffold_MS_twins_infinity.pdf", graph_inf, width=12, height= 8, units="in", useDingbats=FALSE)

# overlay expression of individual signature markers
graph_CD88 <- ggraph(l_g_inf) +
  geom_edge_link(alpha=.15) + 
  geom_node_point(data = subset(l_g_inf, type== "cluster"), aes(size= popsize, fill= CD88), shape=21)+
  geom_node_point(data= subset(l_g_inf, type== "landmark"), shape=21, size=4, fill = "black") +
  scale_fill_gradientn(colors = color_grad_flow2, limits = c(0,1)) +
  scale_size_continuous(range = c(4, 14)) +
  geom_node_text(data= subset(l_g_inf, type == "landmark"), aes(label = name), size = 6, repel=T)+
  theme_graph(base_family = 'Helvetica') + theme(legend.position="none")
graph_CD88
ggsave(filename="Scaffold_MS_twins_infinity_CD88.pdf", graph_CD88, width=12, height= 8, units="in", useDingbats=FALSE)

graph_HLADQ <- ggraph(l_g_inf) +
  geom_edge_link(alpha=.15) + 
  geom_node_point(data = subset(l_g_inf, type== "cluster"), aes(size= popsize, fill= `HLA.DQ`), shape=21)+
  geom_node_point(data= subset(l_g_inf, type== "landmark"), shape=21, size=4, fill = "black") +
  scale_fill_gradientn(colors = color_grad_flow2, limits = c(0,1)) +
  scale_size_continuous(range = c(4, 14)) +
  geom_node_text(data= subset(l_g_inf, type == "landmark"), aes(label = name), size = 6, repel=T)+
  theme_graph(base_family = 'Helvetica') + theme(legend.position="none")
graph_HLADQ
ggsave(filename="Scaffold_MS_twins_infinity_HLADQ.pdf", graph_HLADQ, width=12, height= 8, units="in", useDingbats=FALSE)

graph_Fcer1a <- ggraph(l_g_inf) +
  geom_edge_link(alpha=.15) + 
  geom_node_point(data = subset(l_g_inf, type== "cluster"), aes(size= popsize, fill= `FceRIa`), shape=21)+
  geom_node_point(data= subset(l_g_inf, type== "landmark"), shape=21, size=4, fill = "black") +
  scale_fill_gradientn(colors = color_grad_flow2, limits = c(0,1)) +
  scale_size_continuous(range = c(4, 14)) +
  geom_node_text(data= subset(l_g_inf, type == "landmark"), aes(label = name), size = 6, repel=T)+
  theme_graph(base_family = 'Helvetica') + theme(legend.position="none")
graph_Fcer1a
ggsave(filename="Scaffold_MS_twins_infinity_Fcer1a.pdf", graph_Fcer1a, width=12, height= 8, units="in", useDingbats=FALSE)

graph_CD89 <- ggraph(l_g_inf) +
  geom_edge_link(alpha=.15) + 
  geom_node_point(data = subset(l_g_inf, type== "cluster"), aes(size= popsize, fill= `CD89`), shape=21)+
  geom_node_point(data= subset(l_g_inf, type== "landmark"), shape=21, size=4, fill = "black") +
  scale_fill_gradientn(colors = color_grad_flow2, limits = c(0,1)) +
  scale_size_continuous(range = c(4, 14)) +
  geom_node_text(data= subset(l_g_inf, type == "landmark"), aes(label = name), size = 6, repel=T)+
  theme_graph(base_family = 'Helvetica') + theme(legend.position="none")
graph_CD89
ggsave(filename="Scaffold_MS_twins_infinity_CD89.pdf", graph_CD89, width=12, height= 8, units="in", useDingbats=FALSE)