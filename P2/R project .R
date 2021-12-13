#PART 1

#PART 2

#install useful packages (not useful if already installed)
install.packages("data.table")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("viridis")

library(data.table)
library(pheatmap)
library(viridis)
library("ggplot2")

#replace by your path to the project
setwd("/Users/Djay/Desktop/BA5/G&g/R-project/P2")

#PART 3

# 1
#symetric matrix, approx the same size 
#value => contact matrix = interaction ? 

#10kb
Rv1_in <- fread("data/HiC/22Rv1_chr12_10kb_hic_matrix.txt", sep = "\t", header=T, data.table = F, stringsAsFactors = F)
rownames(Rv1_in) = Rv1_in$V1
Rv1_in = data.matrix(Rv1_in[,2:ncol(Rv1_in)])
Rv1 = Rv1_in[1:500,1:500]
C42B_in <- fread("data/HiC/C42B_chr12_10kb_hic_matrix.txt", sep = "\t", header=T, data.table = F, stringsAsFactors = F)
rownames(C42B_in) = C42B_in$V1
C42B_in = data.matrix(C42B_in[,2:ncol(C42B_in)])
C42B = C42B_in[1:500,1:500]
RWPE1_in <- fread("data/HiC/RWPE1_chr12_10kb_hic_matrix.txt", sep = "\t", header=T, data.table = F, stringsAsFactors = F)
rownames(RWPE1_in) = RWPE1_in$V1
RWPE1_in = data.matrix(RWPE1_in[,2:ncol(RWPE1_in)])
RWPE1 = RWPE1_in[1:500,1:500]

#40kb
Rv140_in <- fread("data/HiC/22Rv1_chr12_40kb_hic_matrix.txt", sep = "\t", header=T, data.table = F, stringsAsFactors = F)
rownames(Rv140_in) = Rv140_in$V1
Rv140_in = data.matrix(Rv140_in[,2:ncol(Rv140_in)])
Rv140 = Rv140_in[1:500,1:500]
C42B40_in <- fread("data/HiC/C42B_chr12_40kb_hic_matrix.txt", sep = "\t", header=T, data.table = F, stringsAsFactors = F)
rownames(C42B40_in) = C42B40_in$V1
C42B40_in = data.matrix(C42B40_in[,2:ncol(C42B40_in)])
C42B40 = C42B40_in[1:500,1:500]
RWPE140_in <- fread("data/HiC/RWPE1_chr12_40kb_hic_matrix.txt", sep = "\t", header=T, data.table = F, stringsAsFactors = F)
rownames(RWPE140_in) = RWPE140_in$V1
RWPE140_in = data.matrix(RWPE140_in[,2:ncol(RWPE140_in)])
RWPE140 = RWPE140_in[1:500,1:500]

#2

Rv1map = log2(Rv1+1)
C42Bmap = log2(C42B+1)
RWPE1map = log2(RWPE1+1)

breaksList = seq(0, 15, by = 1)

png(filename = "pheatmap RV1.png")
pheatmap(Rv1map, labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap C42B.png")
pheatmap(C42Bmap, labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks =breaksList)
dev.off()
png(filename = "pheatmap RWPE1.png")
pheatmap(RWPE1map, labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()


#PART 4

#1

Vanilla_coverage = function(matrix) {
  R = diag(1/rowSums(matrix))
  #we have C = R cause it's symetric
  M = R %*% matrix %*% R
  M = M * sum(matrix)/sum(M)
  return(M)
}

Rv1_norm = Vanilla_coverage(Rv1)
C42B_norm = Vanilla_coverage(C42B)
RWPE1_norm = Vanilla_coverage(RWPE1)
Rv140_norm = Vanilla_coverage(Rv140)
C42B40_norm = Vanilla_coverage(C42B40)
RWPE140_norm = Vanilla_coverage(RWPE140)

#2

RWPE140map = log2(RWPE140+1)
RWPE1_norm_map = log2(RWPE1_norm+1)
RWPE140_norm_map = log2(RWPE140_norm+1)


png(filename = "pheatmap RWP1 norm.png")
pheatmap(RWPE1_norm_map, labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap RWP1 40.png")
pheatmap(RWPE140map, labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap RWP1 40 norm.png")
pheatmap(RWPE140_norm_map, labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()

#3 BONUS

#PART 5

#1

RWPE1_N_cut = data.matrix(RWPE1_in[12380:12779,12380:12779])
C42B_N_cut = data.matrix(C42B_in[12380:12779,12380:12779])
Rv1_N_cut = data.matrix(Rv1_in[12380:12779,12380:12779])

RWPE1_N_cut_map = log2(RWPE1_N_cut+1)
C42B_N_cut_map = log2(C42B_N_cut+1)
Rv1_N_cut_map = log2(Rv1_N_cut+1)

png(filename = "pheatmap RWP1 norm cut.png")
pheatmap(RWPE1_N_cut_map, labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap C42B norm cut.png")
pheatmap(C42B_N_cut_map, labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap Rv1 norm cut .png")
pheatmap(Rv1_N_cut_map, labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()

#2

#3
RW_RV = RWPE1_N_cut_map-Rv1_N_cut_map
RW_CB = RWPE1_N_cut_map-C42B_N_cut_map
RV_CB = Rv1_N_cut_map-C42B_N_cut_map

png(filename = "pheatmap RW_RV.png")
pheatmap(RW_RV, labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap RW_CB.png")
pheatmap(RW_CB, labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap RV_CB .png")
pheatmap(RV_CB, labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()

#PART 6 

#1
#2
#3

#PART 7 BONUS