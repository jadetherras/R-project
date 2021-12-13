#PART 1

#PART 2

#install useful packages (not useful if already installed)
install.packages("data.table")
install.packages("ggplot2")
install.packages("pheatmap")

library(data.table)
library(pheatmap)
library("ggplot2")

#replace by your path to the project
setwd("/Users/Djay/Desktop/BA5/G&g/R-project/P2")

#PART 3

# 1
#symetric matrix, approx the same size 
#value => contact matrix = interaction ? 

#10kb
Rv1 <- fread("data/HiC/22Rv1_chr12_10kb_hic_matrix.txt", sep = "\t", header=T, data.table = F, stringsAsFactors = F)
rownames(Rv1) = Rv1$V1
Rv1 = Rv1[,2:ncol(Rv1)]
C42B <- fread("data/HiC/C42B_chr12_10kb_hic_matrix.txt", sep = "\t", header=T, data.table = F, stringsAsFactors = F)
rownames(C42B) = C42B$V1
C42B = C42B[,2:ncol(C42B)]
RWPE1 <- fread("data/HiC/RWPE1_chr12_10kb_hic_matrix.txt", sep = "\t", header=T, data.table = F, stringsAsFactors = F)
rownames(RWPE1) = RWPE1$V1
RWPE1 = RWPE1[,2:ncol(RWPE1)]

#40kb
Rv140 <- fread("data/HiC/22Rv1_chr12_40kb_hic_matrix.txt", sep = "\t", header=T, data.table = F, stringsAsFactors = F)
rownames(Rv140) = Rv140$V1
Rv140 = Rv140[,2:ncol(Rv140)]
C42B40 <- fread("data/HiC/C42B_chr12_40kb_hic_matrix.txt", sep = "\t", header=T, data.table = F, stringsAsFactors = F)
rownames(C42B40) = C42B40$V1
C42B40 = C42B40[,2:ncol(C42B40)]
RWPE140 <- fread("data/HiC/RWPE1_chr12_40kb_hic_matrix.txt", sep = "\t", header=T, data.table = F, stringsAsFactors = F)
rownames(RWPE140) = RWPE140$V1
RWPE140 = RWPE140[,2:ncol(RWPE140)]

#2
Rv1map = Rv1[1:500,1:500]
C42Bmap = C42B[1:500,1:500]
RWPE1map = RWPE1[1:500,1:500]
for (i in 1:500) {
  for (j in 1:500) {
    if(Rv1map[i,j] != 0){
      Rv1map[i,j] = log2(Rv1map[i,j])
    }
    if(C42Bmap[i,j] != 0){
      C42Bmap[i,j] = log2(C42Bmap[i,j])
    } 
    if(RWPE1map[i,j] != 0){
      RWPE1map[i,j] = log2(RWPE1map[i,j])
    } 
  }
}

#pheatmap(Rv1map, labels_row = '', labels_col = '', breaks = seq(0,10,by = 0.5), color = red)

png(filename = "pheatmap.png")
par(mfrow = c(1,3))
pheatmap(Rv1map, labels_row = '', labels_col = '')
pheatmap(C42Bmap, labels_row = '', labels_col = '')
pheatmap(RWPE1map, labels_row = '', labels_col = '')
dev.off()

#PART 4

#1

Vanilla_coverage = function(matrix) {
  R = diag(1/rowSums(matrix))
  #we have C = R cause it's symetric
  print(R)
  M = R %*% matrix %*% R
  
  M = sum(matrix)/sum(M)
  
  return(M)
}

Rv1_norm = Vanilla_coverage(Rv1)
C42B_norm = Vanilla_coverage(C42B)
RWPE1_norm = Vanilla_coverage(RWPE1)
Rv140_norm = Vanilla_coverage(RV140)
C42B40_norm = Vanilla_coverage(C42B40)
RWPE140_norm = Vanilla_coverage(RWPE140)

v = c(1,1,1,1)
m = matrix(v,2,2)

Vanilla_coverage(matrix)

#2


RWPE1map = RWPE1[1:500,1:500]
RWPE140map = RWPE140[1:500,1:500]
RWPE1_norm_map = RWPE1_norm[1:500,1:500]
RWPE140_norm_map = RWPE140_norm[1:500,1:500]
for (i in 1:500) {
  for (j in 1:500) {
    if(RWPE1map[i,j] != 0){
      RWPE1map[i,j] = log2(RWPE1map[i,j])
    }
    if(RWPE140map[i,j] != 0){
      RWPE140map[i,j] = log2(RWPE140map[i,j])
    } 
    if(RWPE1_norm_map[i,j] != 0){
      RWPE1_norm_map[i,j] = log2(RWPE1_norm_map[i,j])
    }
    if(RWPE140_norm_map[i,j] != 0){
      RWPE140_norm_map[i,j] = log2(RWPE140_norm_map[i,j])
    }
  }
}
pheatmap(RWPE1map, labels_row = '', labels_col = '')
pheatmap(Rv1map_norm, labels_row = '', labels_col = '')
pheatmap(RWPE140map, labels_row = '', labels_col = '')
pheatmap(Rv1map40_norm_map, labels_row = '', labels_col = '')



#3 BONUS

#PART 5
#PART 6 
#PART 7