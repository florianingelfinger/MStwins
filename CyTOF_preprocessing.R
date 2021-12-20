# Preprocessing of mass cytometry data presented in Ingelfinger, Gerdes et al. "Twin study reveals non-heritable immune perturbations in multiple sclerosis"
# This code is losely based on Hartmann, F. J. et al. High-dimensional single-cell analysis reveals the immune signature of narcolepsy. J. Exp. Med. (2016) doi:10.1084/jem.20160897.

# load library
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


# set working directory
setwd("/Immunology_shares/Becherlab/People/Ingelfinger/DATA/Multiple Sclerosis/Twin Study/Pooled runs/Surface/0_transformation")


# read fcs files of the two batches
setwd("/Immunology_shares/Becherlab/People/Ingelfinger/DATA/Multiple Sclerosis/Twin Study/2019-02-26_twins_surf/normed/gated/debarcoded/new")
files_names_surf1 <- list.files(getwd(), pattern='.fcs$', full=FALSE)

fs_surf1  <- read.ncdfFlowSet(files = list.files(path = getwd(), pattern = ".fcs$"),
                             transformation = F,
                             phenoData = ,
                             truncate_max_range = F)

setwd("/Immunology_shares/Becherlab/People/Ingelfinger/DATA/Multiple Sclerosis/Twin Study/2019-04-29_twins_surf/normed/gated/debarcoded/new")

files_names_surf2 <- list.files(getwd(), pattern='.fcs$', full=FALSE)

fs_surf2  <- read.ncdfFlowSet(files = list.files(path = getwd(), pattern = ".fcs$" ),
                             transformation = F,
                             phenoData = ,
                             truncate_max_range = F)


# combine data into a matrix and rename
data_surf1 <- fsApply(fs_surf1, exprs)
data_surf2 <- fsApply(fs_surf2, exprs)


cn_surf1 <- colnames(data_surf1)
cn_surf1 <- c("CD11b", "CD20", "FoxP3", "CD207", "CD45R0", "CD86", "CCR7", "Ki67", "IgD", "CD38", "CD45RA", 
              "CD268", "CXCR3", "CCR4", "CD69", "CD127", "CD14", "CD5", "CD19", "CCR2", "CCR5", "CD4", "CD8", 
              "CD16", "CD27", "CCR6", "CD11c", "CD25", "TCRgd", "CD3", "CD116", "CD33", "CD138", "IgM", "CD56", "HLA-DR")
length(colnames(data_surf1)) == length(cn_surf1)

data_surf2 <- data_surf2[,c(1:18, 20:37)]
cn_surf2 <- colnames(data_surf2)
cn_surf2 <- cn_surf1
colnames(data_surf1)[c(1:36)] <- cn_surf1

colnames(data_surf2) <- cn_surf2
head(data_surf2)
dim(data_surf2)


# add an ID for each cell and each patient
dim.vec_surf1 <- fsApply(fs_surf1, dim)
dim.vec_surf1 <- as.numeric(dim.vec_surf1[,1])
gate.source_surf1 <- as.vector(x = NULL)
for(i in 1:length(dim.vec_surf1)) {temp.source <- rep(i, dim.vec_surf1[i])
gate.source_surf1 <- c(gate.source_surf1, temp.source)}

data_surf1 <- cbind(data_surf1, gate.source_surf1)
cell_id_surf1 <- 1:dim(data_surf1)[1]
data_surf1 <- cbind(data_surf1, cell_id_surf1)
colnames(data_surf1)[37:38] <- c("gate_source", "cell_id")

dim.vec_surf2 <- fsApply(fs_surf2, dim)
dim.vec_surf2 <- as.numeric(dim.vec_surf2[,1])
gate.source_surf2 <- as.vector(x = NULL)
for(i in 1:length(dim.vec_surf2)) {temp.source <- rep(i, dim.vec_surf2[i])
gate.source_surf2 <- c(gate.source_surf2, temp.source)}

data_surf2 <- cbind(data_surf2, gate.source_surf2)
cell_id_surf2 <- 1:dim(data_surf2)[1]
data_surf2 <- cbind(data_surf2, cell_id_surf2)
colnames(data_surf2)[37:38] <- c("gate_source", "cell_id")

# continue with gate.source and cell_id for run 2 to avoid duplicates when analyzing together
data_surf2[, "gate_source"] <-  data_surf2[, "gate_source"]+max(data_surf1[, "gate_source"])
data_surf2[, "cell_id"] <-  data_surf2[, "cell_id"]+max(data_surf1[, "cell_id"])


### extract the respective sample_id from .fcs filename
sample_id_surf1 <- as.numeric(gsub("[^0-9]", "", files_names_surf1))
sample_gate_surf1 <- cbind(unique(data_surf1[,"gate_source"]), sample_id_surf1)
colnames(sample_gate_surf1) <- c("gate_source", "sample_id")
data_surf1 <- merge(data_surf1, sample_gate_surf1, by= "gate_source")

sample_id_surf2 <- as.numeric(gsub("[^0-9]", "", files_names_surf2))
sample_gate_surf2 <- cbind(unique(data_surf2[,"gate_source"]), sample_id_surf2)
colnames(sample_gate_surf2) <- c("gate_source", "sample_id")
sample_gate_surf2[,2] <- as.numeric(sample_gate_surf2[,2]) + 68
data_surf2 <- merge(data_surf2, sample_gate_surf2, by= "gate_source")

# convert to matrix
data_matrix_surf1 <- data.matrix(data_surf1)
data_matrix_surf2 <- data.matrix(data_surf2)


# transform data with global arcsin cofactor of 5 and adjustments of individual channels
# note: cofactors have been chosen on na indiividual basis depending on the separation of positive and negative fraction
# as observed in the output pdf file.

asinh_scale <- 5

data.trans_surf1 <- asinh(data_matrix_surf1/ asinh_scale)
data.trans_surf1[,"CD3"] <- asinh(5*sinh(data.trans_surf1[,"CD3"])/20)
data.trans_surf1[,"CD20"] <- asinh(5*sinh(data.trans_surf1[,"CD20"])/20)
data.trans_surf1[,"FoxP3"] <- asinh(5*sinh(data.trans_surf1[,"FoxP3"])/50)
data.trans_surf1[,"CD207"] <- asinh(5*sinh(data.trans_surf1[,"CD207"])/10)
data.trans_surf1[,"CD45R0"] <- asinh(5*sinh(data.trans_surf1[,"CD45R0"])/20)
data.trans_surf1[,"CCR7"] <- asinh(5*sinh(data.trans_surf1[,"CCR7"])/12)
data.trans_surf1[,"CD45RA"] <- asinh(5*sinh(data.trans_surf1[,"CD45RA"])/20)
data.trans_surf1[,"CD268"] <- asinh(5*sinh(data.trans_surf1[,"CD268"])/3)
data.trans_surf1[,"CXCR3"] <- asinh(5*sinh(data.trans_surf1[,"CXCR3"])/8)
data.trans_surf1[,"CCR4"] <- asinh(5*sinh(data.trans_surf1[,"CCR4"])/13)
data.trans_surf1[,"CD5"] <- asinh(5*sinh(data.trans_surf1[,"CD5"])/12)
data.trans_surf1[,"CCR2"] <- asinh(5*sinh(data.trans_surf1[,"CCR2"])/12)
data.trans_surf1[,"CCR5"] <- asinh(5*sinh(data.trans_surf1[,"CCR5"])/8)
data.trans_surf1[,"CD4"] <- asinh(5*sinh(data.trans_surf1[,"CD4"])/20)
data.trans_surf1[,"CD8"] <- asinh(5*sinh(data.trans_surf1[,"CD8"])/30)
data.trans_surf1[,"CD16"] <- asinh(5*sinh(data.trans_surf1[,"CD16"])/10)
data.trans_surf1[,"CD27"] <- asinh(5*sinh(data.trans_surf1[,"CD27"])/7)
data.trans_surf1[,"TCRgd"] <- asinh(5*sinh(data.trans_surf1[,"TCRgd"])/8)
data.trans_surf1[,"CD56"] <- asinh(5*sinh(data.trans_surf1[,"CD56"])/3)


data.trans_surf2 <- asinh(data_matrix_surf2 / asinh_scale)
data.trans_surf2[,"CD3"] <- asinh(5*sinh(data.trans_surf2[,"CD3"])/25)
data.trans_surf2[,"CD20"] <- asinh(5*sinh(data.trans_surf2[,"CD20"])/10)
data.trans_surf2[,"FoxP3"] <- asinh(5*sinh(data.trans_surf2[,"FoxP3"])/45)
data.trans_surf2[,"CD45R0"] <- asinh(5*sinh(data.trans_surf2[,"CD45R0"])/20)
data.trans_surf2[,"CD45RA"] <- asinh(5*sinh(data.trans_surf2[,"CD45RA"])/15)
data.trans_surf2[,"CCR4"] <- asinh(5*sinh(data.trans_surf2[,"CCR4"])/10)
data.trans_surf2[,"CD5"] <- asinh(5*sinh(data.trans_surf2[,"CD5"])/10)
data.trans_surf2[,"CCR2"] <- asinh(5*sinh(data.trans_surf2[,"CCR2"])/10)
data.trans_surf2[,"CD4"] <- asinh(5*sinh(data.trans_surf2[,"CD4"])/20)
data.trans_surf2[,"CD8"] <- asinh(5*sinh(data.trans_surf2[,"CD8"])/27)
data.trans_surf2[,"CD16"] <- asinh(5*sinh(data.trans_surf2[,"CD16"])/10)
data.trans_surf2[,"CD27"] <- asinh(5*sinh(data.trans_surf2[,"CD27"])/3)
data.trans_surf2[,"HLA-DR"] <- asinh(5*sinh(data.trans_surf2[,"HLA-DR"])/8)


data.trans_surf1[,"gate_source"] <- data_surf1[,"gate_source"]
data.trans_surf1[,"cell_id"] <- data_surf1[,"cell_id"]
data.trans_surf1[,"sample_id"] <- data_surf1[,"sample_id"]
data.trans_surf2[,"gate_source"] <- data_surf2[,"gate_source"]
data.trans_surf2[,"cell_id"] <- data_surf2[,"cell_id"]
data.trans_surf2[,"sample_id"] <- data_surf2[,"sample_id"]


# normalize based on the 99.99th percentile to have everything between 0 and 1
data.trans.new_surf1 <- data.trans_surf1
per.vector_surf1 <- apply(data.trans.new_surf1, 2, function(x) quantile(x, 0.9999, names = F))
data.trans.new_surf1 <- t(t(data.trans.new_surf1) / as.numeric(per.vector_surf1))
data.trans.new_surf1[,"gate_source"] <- data_surf1[,"gate_source"]
data.trans.new_surf1[,"cell_id"] <- data_surf1[,"cell_id"]
data.trans.new_surf1[,"sample_id"] <- data_surf1[,"sample_id"]


data.trans.new_surf2 <- data.trans_surf2
per.vector_surf2 <- apply(data.trans.new_surf2, 2, function(x) quantile(x, 0.9999, names = F))
data.trans.new_surf2 <- t(t(data.trans.new_surf2) / as.numeric(per.vector_surf2))
data.trans.new_surf2[,"gate_source"] <- data_surf2[,"gate_source"]
data.trans.new_surf2[,"cell_id"] <- data_surf2[,"cell_id"]
data.trans.new_surf2[,"sample_id"] <- data_surf2[,"sample_id"]


# adjustments of single channels
data.trans.new_surf2[,"CD27"] <- data.trans_surf2[,"CD27"]/quantile(data.trans_surf2[,"CD27"], 0.99999, names = F)


# select batch control samples (identical pbmc samples outside the twinn cohort to assess residual batch effect between individual runs)
data.trans.norm1 <- data.frame(data.trans.new_surf1[data.trans.new_surf1[,"gate_source"]==67,])
data.trans.norm2 <- data.frame(data.trans.new_surf2[data.trans.new_surf2[,"gate_source"]==135,])

# check means per channel between batches
means1 <- apply(data.trans.norm1, 2, mean)
means2 <- apply(data.trans.norm1, 2, mean)

# plot CD3 vs all remaining markers of normalization controls
db13 <- c("#1B9E77","gold" , "#7570B3", "#E7298A", "#D95F02" ,  "lightpink1", "#666666",
          "#57B0FF", "#66CD00", "#1F78B4", "#B2DF8A", "firebrick1", 
          "seagreen1", "brown" , "#6A3D9A", "mediumorchid", "grey" ,
          "#BC80BD", "#80B1D3")

norm_surf1_list <- list()
norm_surf2_list <- list()
for(i in 1:length(colnames(data.trans.norm1))){
  norm_surf1 <- ggplot(data.trans.norm1, aes_string(x = colnames(data.trans.norm1)[i], y = "CD3")) +
    geom_hex(bins = 100) +
    xlim(-0.1, 1.5) +
    ylim(-0.1, 1.5) +
    coord_fixed(ratio = 1) +
    scale_fill_gradientn(colours = db13, trans = "sqrt") +
    theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 3), axis.title = element_text(size = 3))
  norm_surf1_list[[i]] <-  norm_surf1
  
  norm_surf2 <- ggplot(data.trans.norm2, aes_string(x = colnames(data.trans.norm2)[i], y = "CD3")) +
    geom_hex(bins = 100) +
    xlim(-0.1, 1.5) +
    ylim(-0.1, 1.5) +
    coord_fixed(ratio = 1) +
    scale_fill_gradientn(colours = db13, trans = "sqrt") +
    theme_bw()+ theme(legend.position="none", axis.text = element_text(size = 3), axis.title = element_text(size = 3))
  norm_surf2_list[[i]] <-  norm_surf2
}

setwd("/Immunology_shares/Becherlab/People/Ingelfinger/DATA/Multiple Sclerosis/Twin Study/Pooled runs/Surface/0_transformation")

g_m1 <- arrangeGrob(grobs = norm_surf1_list, ncol = 1, top=textGrob("Experiment 1", gp=gpar(fontsize=6)))
g_m2 <- arrangeGrob(grobs = norm_surf2_list, ncol = 1, top=textGrob("Experiment 2", gp=gpar(fontsize=6)))
m <- arrangeGrob(g_m1, g_m2, ncol= 2)
ggsave(file="marker_nrom2.pdf", m, width = 2, height = 34)

# merge both runs into a single dataframe
data.trans.new_surf1 <- cbind(data.trans.new_surf1, rep(1, nrow(data.trans.new_surf1)))
data.trans.new_surf2 <- cbind(data.trans.new_surf2, rep(2, nrow(data.trans.new_surf2)))
colnames(data.trans.new_surf1)[40] <- "n_run"
colnames(data.trans.new_surf2)[40] <- "n_run"

data <-  rbind(data.trans.new_surf1, data.trans.new_surf2)
save(data, file="2019_06_05_data_twins_merged.RDa")

