#PART 1

#PART 2

#install useful packages (not useful if already installed)
install.packages("data.table")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("viridis")
install.packages("gridExtra")

library("gridExtra")
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

max_list = 15
breaksList = seq(0, max_list, by = 1)

png(filename = "pheatmap RV1.png")
pheatmap(Rv1map, color=colorRampPalette(c("white", "red"))(max_list), main = 'pheatmap 22Rv1 10kb', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap C42B.png")
pheatmap(C42Bmap,  color=colorRampPalette(c("white", "red"))(max_list),main = 'pheatmap C42B,10kb', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks =breaksList)
dev.off()
png(filename = "pheatmap RWPE1.png")
pheatmap(RWPE1map, color=colorRampPalette(c("white", "red"))(max_list),main = 'pheatmap RWPE1 10kb', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
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
pheatmap(RWPE1_norm_map, color=colorRampPalette(c("white", "red"))(max_list),main = 'pheatmap RWPE1 10kb normalised', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap RWP1 40.png")
pheatmap(RWPE140map, color=colorRampPalette(c("white", "red"))(max_list),main = 'pheatmap RWPE1 40kb', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap RWP1 40 norm.png")
pheatmap(RWPE140_norm_map, color=colorRampPalette(c("white", "red"))(max_list),main = 'pheatmap RWPE1 40kb normalised', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()

#3 BONUS

#PART 5

#1

start_num <- which(colnames(RWPE1_in) == "chr12:127000000:127010000")
stop_num <- which(rownames(RWPE1_in) == "chr12:130990000:131000000")

RWPE1_N_cut = Vanilla_coverage(data.matrix(RWPE1_in[start_num:stop_num, start_num:stop_num]))
C42B_N_cut = Vanilla_coverage(data.matrix(C42B_in[start_num:stop_num, start_num:stop_num]))
Rv1_N_cut = Vanilla_coverage(data.matrix(Rv1_in[start_num:stop_num, start_num:stop_num]))


RWPE1_N_cut_map = log2(RWPE1_N_cut+1)
C42B_N_cut_map = log2(C42B_N_cut+1)
Rv1_N_cut_map = log2(Rv1_N_cut+1)

png(filename = "pheatmap C42B norm cut.png")
pheatmap(C42B_N_cut_map, color=colorRampPalette(c("white", "red"))(max_list),main = 'pheatmap C42B 10kb normalised ', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap Rv1 norm cut.png")
pheatmap(Rv1_N_cut_map, color=colorRampPalette(c("white", "red"))(max_list),main = 'pheatmap 22Rv1 10kb normalised ', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap RWP1 norm cut.png")
pheatmap(RWPE1_N_cut_map, color=colorRampPalette(c("white", "red"))(max_list),main = 'pheatmap RWPE1 10kb normalised ', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()

#2

#3

#signe ? 
RW_RV = sign(RWPE1_N_cut-Rv1_N_cut)*log2(abs(RWPE1_N_cut-Rv1_N_cut) +1)
RW_CB = sign(RWPE1_N_cut-C42B_N_cut)*log2(abs(RWPE1_N_cut-C42B_N_cut)+1)
RV_CB = sign(Rv1_N_cut-C42B_N_cut)*log2(abs(Rv1_N_cut-C42B_N_cut)+1)
#RW_RV = log2(abs(RWPE1_N_cut-Rv1_N_cut) +1)
#RW_CB = log2(abs(RWPE1_N_cut-C42B_N_cut)+1)
#RV_CB = log2(abs(Rv1_N_cut-C42B_N_cut)+1)

png(filename = "pheatmap RW_RV.png")
pheatmap(RW_RV, color=magma(max_list),main = 'RWPE1 & 22Rv1', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap RW_CB.png")
pheatmap(RW_CB, color=magma(max_list),main = 'RWPE1 & C42B', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap RV_CB.png")
pheatmap(RV_CB, color=magma(max_list),main = '22Rv1 & C42B', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()

#PART 6 

#1

r = 40*10^3
w = 2*10^6
nb_bin = w/r

Directionality_Index = function(matrix) {
  DI = c()
  I = ncol(matrix)
  for (i in 1:I) {
    A = sum(matrix[i,max(1,(i-nb_bin)):max(1,(i-1))])
    B = sum(matrix[i,min(I,(i+1)):min(I,(i+nb_bin))])
    E = (A + B)/2
    DI = c(DI, (B-A)*((A-E)^2/E+(B-E)^2/E)/abs(B-A))
  }
  return(data.frame("Bins" = c(1:I), "DI" = DI))
}

#2 and 3

start_num <- which(colnames(Rv140_in) == "chr12:127000000:127040000")
stop_num <- which(rownames(Rv140_in) == "chr12:130960000:131000000")

# with all data

for_DI_Rv140 = Vanilla_coverage(Rv140_in)
for_DI_C42B40 = Vanilla_coverage(C42B40_in)
for_DI_RWPE140 = Vanilla_coverage(RWPE140_in)

for_DI_Rv140_map = log2(for_DI_Rv140+1)
for_DI_C42B40_map = log2(for_DI_C42B40+1)
for_DI_RWPE140_map = log2(for_DI_RWPE140+1)

max_list = 15
breaksList = seq(0, max_list, by = 1)

png(filename = "pheatmap for_DI_Rv140.png")
pheatmap(for_DI_Rv140_map,color=colorRampPalette(c("white", "red"))(max_list),main = 'for_DI_Rv140', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap for_DI_C42B40.png")
pheatmap(for_DI_C42B40_map,color=colorRampPalette(c("white", "red"))(max_list),main = 'for_DI_C42B40', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()
png(filename = "pheatmap for_DI_RWPE140.png")
pheatmap(for_DI_RWPE140_map,color=colorRampPalette(c("white", "red"))(max_list),main = 'for_DI_RWPE140', labels_row = '', labels_col = '', cluster_cols = F, cluster_rows  = F,breaks = breaksList)
dev.off()

DI_Rv1_all  = Directionality_Index(for_DI_Rv140)[start_num:stop_num,]
DI_C42B_all  = Directionality_Index(for_DI_C42B40)[start_num:stop_num,]
DI_RWPE1_all  = Directionality_Index(for_DI_RWPE140)[start_num:stop_num,]


l = 1

DI_plot_Rv1_all = data.frame(xm = (DI_Rv1_all$Bins-start_num)*l, xM = (DI_Rv1_all$Bins-start_num+1)*l, ym = 0, yM = DI_Rv1_all$DI)
DI_plot_C42B_all = data.frame(xm = (DI_C42B_all$Bins-start_num)*l, xM = (DI_C42B_all$Bins-start_num+1)*l, ym = 0, yM = DI_C42B_all$DI)
DI_plot_RWPE1_all = data.frame(xm = (DI_RWPE1_all$Bins-start_num)*l, xM = (DI_RWPE1_all$Bins-start_num+1)*l, ym = 0, yM = DI_RWPE1_all$DI)

png(filename = "DI_Rv1_all.png")
p <- ggplot() + geom_rect(data=DI_plot_Rv1_all, mapping=aes(xmin = xm, xmax=xM, ymin=ym, ymax=yM, fill = as.factor(sign(yM)))) + ggtitle("DI for 22Rv1 with all the matrix")
p
dev.off()

png(filename = "DI_C42B_all.png")
p <- ggplot() + geom_rect(data=DI_plot_C42B_all, mapping=aes(xmin = xm, xmax=xM, ymin=ym, ymax=yM,fill =  as.factor(sign(yM))))+ ggtitle("DI for C42B with all the matrix")
p
dev.off()

png(filename = "DI_RWPE1_all.png")
p <- ggplot() + geom_rect(data=DI_plot_RWPE1_all, mapping=aes(xmin = xm, xmax=xM, ymin=ym, ymax=yM,fill =  as.factor(sign(yM))))+ ggtitle("DI for RWPE1 with all the matrix")
p
dev.off()

#only with the cut data

DI_Rv1 = Directionality_Index(Vanilla_coverage(Rv140_in[start_num:stop_num, start_num:stop_num]))
DI_C42B = Directionality_Index(Vanilla_coverage(C42B40_in[start_num:stop_num, start_num:stop_num]))
DI_RWPE1 = Directionality_Index(Vanilla_coverage(RWPE140_in[start_num:stop_num, start_num:stop_num]))

DI_plot_Rv1 = data.frame(xm = (DI_Rv1$Bins-start_num)*l, xM = (DI_Rv1$Bins-start_num+1)*l, ym = 0, yM = DI_Rv1$DI)
DI_plot_C42B = data.frame(xm = (DI_C42B$Bins-1)*l, xM = (DI_C42B$Bins)*l, ym = 0, yM = DI_C42B$DI)
DI_plot_RWPE1 = data.frame(xm = (DI_RWPE1$Bins-1)*l, xM = (DI_RWPE1$Bins)*l, ym = 0, yM = DI_RWPE1$DI)


png(filename = "DI_Rv1.png")
p <- ggplot() + geom_rect(data=DI_plot_Rv1, mapping=aes(xmin = xm, xmax=xM, ymin=ym, ymax=yM, fill = as.factor(sign(yM)))) + ggtitle("DI for 22Rv1") + coord_cartesian(ylim = c(-2000, 2000))
p
dev.off()

png(filename = "DI_C42B.png")
p <- ggplot() + geom_rect(data=DI_plot_C42B, mapping=aes(xmin = xm, xmax=xM, ymin=ym, ymax=yM,fill = as.factor(sign(yM))))+ ggtitle("DI for C42B") + coord_cartesian(ylim = c(-2000, 2000))
p
dev.off()

png(filename = "DI_RWPE1.png")
p <- ggplot() + geom_rect(data=DI_plot_RWPE1, mapping=aes(xmin = xm, xmax=xM, ymin=ym, ymax=yM,fill = as.factor(sign(yM))))+ ggtitle("DI for RWPE1") + coord_cartesian(ylim = c(-2000, 2000))
p
dev.off()

#4

triangle = function(M) {
  xlast = 0
  x = c()
  y = c()
  g = c()
  nb = 1
  for (i in 2:nrow(M)) {
    if (sign(M$DI[i-1]) != sign(M$DI[i])) {
      x = c(x,xlast,(xlast+i)/2,i)
      xlast = i
      y = c(y, 0,1,0)
      g = c(g,nb,nb,nb)
      nb = nb +1
    }
  }
  x = c(x,xlast,(xlast+nrow(M))/2,nrow(M))
  y = c(y, 0,1,0)
  g = c(g,nb,nb,nb)
  return(data.frame('x'=x,'y'=y, 'g'=g))
}

# with only selected data

png(filename = "triangle_and_tract_22RV1.png")
p1 <- ggplot(DI_Rv1, aes(x = Bins, y = sign(DI))) + geom_point() + ggtitle("tract for 22Rv1")
p2 <- ggplot() + geom_polygon(data = triangle(DI_Rv1), mapping=aes(x=x, y=y, group=g)) + ggtitle("triangle for 22Rv1")
grid.arrange(p1,p2, nrow=2)
dev.off()

png(filename = "triangle_and_tract_C42B.png")
p1 <- ggplot(DI_C42B, aes(x = Bins, y = sign(DI))) + geom_point() + ggtitle("tract for C42B")
p2 <- ggplot() + geom_polygon(data = triangle(DI_C42B), mapping=aes(x=x, y=y, group=g)) + ggtitle("triangle for C42B")
grid.arrange(p1,p2, nrow=2)
dev.off()

png(filename = "triangle_and_tract_RWPE1.png")
p1 <- ggplot(DI_RWPE1, aes(x = Bins, y = sign(DI))) + geom_point() + ggtitle("tract for RWPE1")
p2 <- ggplot() + geom_polygon(data = triangle(DI_RWPE1), mapping=aes(x=x, y=y, group=g)) + ggtitle("triangle for RWPE1")
grid.arrange(p1,p2, nrow=2)
dev.off()

#with all data

png(filename = "triangle_and_tract_22RV1_all.png")
p1 <- ggplot(DI_Rv1_all, aes(x = Bins, y = sign(DI))) + geom_point() + ggtitle("tract for 22Rv1")
p2 <- ggplot() + geom_polygon(data = triangle(DI_Rv1_all), mapping=aes(x=x, y=y, group=g)) + ggtitle("triangle for 22Rv1")
grid.arrange(p1,p2, nrow=2)
dev.off()

png(filename = "triangle_and_tract_C42B_all.png")
p1 <- ggplot(DI_C42B_all, aes(x = Bins, y = sign(DI))) + geom_point() + ggtitle("tract for C42B")
p2 <- ggplot() + geom_polygon(data = triangle(DI_C42B_all), mapping=aes(x=x, y=y, group=g)) + ggtitle("triangle for C42B")
grid.arrange(p1,p2, nrow=2)
dev.off()

png(filename = "triangle_and_tract_RWPE1_all.png")
p1 <- ggplot(DI_RWPE1_all, aes(x = Bins, y = sign(DI))) + geom_point() + ggtitle("tract for RWPE1")
p2 <- ggplot() + geom_polygon(data = triangle(DI_RWPE1_all), mapping=aes(x=x, y=y, group=g)) + ggtitle("triangle for RWPE1")
grid.arrange(p1,p2, nrow=2)
dev.off()


#PART 7 BONUS

#2

start_stop = function(u) {
  start_num <- which(u$V2 >= 127000000)[1]
  stop_num <- which(u$V3 > 131000000)[1] -1
  return(u[start_num:stop_num,])
}

Rv1_H3K9me3_chr12 <- start_stop(fread("data/bonus_chipseq/22Rv1_H3K9me3_chr12.bedGraph", sep = "\t", data.table = F, stringsAsFactors = F))
Rv1_H3K27me3_chr12 <- start_stop(fread("data/bonus_chipseq/22Rv1_H3K27me3_chr12.bedGraph", sep = "\t", data.table = F, stringsAsFactors = F))
Rv1_H3K36me3_chr12 <- start_stop(fread("data/bonus_chipseq/22Rv1_H3K36me3_chr12.bedGraph", sep = "\t", data.table = F, stringsAsFactors = F))

C42B_H3K9me3_chr12 <- start_stop(fread("data/bonus_chipseq/C42B_H3K9me3_chr12.bedGraph", sep = "\t", data.table = F, stringsAsFactors = F))
C42B_H3K27me3_chr12 <- start_stop(fread("data/bonus_chipseq/C42B_H3K27me3_chr12.bedGraph", sep = "\t", data.table = F, stringsAsFactors = F))
C42B_H3K36me3_chr12 <- start_stop(fread("data/bonus_chipseq/C42B_H3K36me3_chr12.bedGraph", sep = "\t", data.table = F, stringsAsFactors = F))

RWPE1_H3K9me3_chr12 <- start_stop(fread("data/bonus_chipseq/RWPE1_H3K9me3_chr12.bedGraph", sep = "\t", data.table = F, stringsAsFactors = F))
RWPE1_H3K27me3_chr12 <- start_stop(fread("data/bonus_chipseq/RWPE1_H3K27me3_chr12.bedGraph", sep = "\t", data.table = F, stringsAsFactors = F))
RWPE1_H3K36me3_chr12 <- start_stop(fread("data/bonus_chipseq/RWPE1_H3K36me3_chr12.bedGraph", sep = "\t", data.table = F, stringsAsFactors = F))


#3

data_plot = function(frame,title){
  plotty = data.frame(xm = frame$V2, xM = frame$V3, ym = 0, yM = log2(frame$V4+1))
  p <- ggplot() + geom_rect(data=plotty , mapping=aes(xmin = xm, xmax=xM, ymin=ym, ymax=yM))+ xlim(127000000, 131000000) + ggtitle(title)
  return(p)
}


png(filename = "Chip-sequencing_Rv1.png")

p1 = data_plot(Rv1_H3K9me3_chr12,"Rv1_H3K9me3_chr12")
p2 = data_plot(Rv1_H3K27me3_chr12,"Rv1_H3K27me3_chr12")
p3 = data_plot(Rv1_H3K36me3_chr12, "Rv1_H3K36me3_chr12")

grid.arrange(p1,p2,p3, nrow = 3)

dev.off()


png(filename= "Chip-sequencing_C42B.png")
p4 = data_plot(C42B_H3K9me3_chr12,"C42B_H3K9me3_chr12")
p5 = data_plot(C42B_H3K27me3_chr12,"C42B_H3K27me3_chr12")
p6 = data_plot(C42B_H3K36me3_chr12,"C42B_H3K36me3_chr12")

grid.arrange(p4,p5,p6, nrow=3)
dev.off()


png(filename= "Chip-sequencing_RWPE1.png")
p7 = data_plot(RWPE1_H3K9me3_chr12,"RWPE1_H3K9me3_chr12")
p8 = data_plot(RWPE1_H3K27me3_chr12,"RWPE1_H3K27me3_chr12")
p9 = data_plot(RWPE1_H3K36me3_chr12,"RWPE1_H3K36me3_chr12")

grid.arrange(p7,p8,p9, nrow = 3)

dev.off()

