#PART 1

#install useful packages (not useful if already installed)
install.packages("data.table")
install.packages("ggplot2")

library(data.table)
library("ggplot2")

#replace by your path to the project
setwd("/Users/Djay/Desktop/BA5/G&g/R-project")

#PART 2

#1

#we decide to download the resources folder and at it at the project level
#replace by your path to the resources
covariates <- fread("resources/covariates.txt", sep=",", header=T, data.table = F, stringsAsFactors = F)
phenotypes = fread("resources/phenotypes.txt", sep=",", header=T, data.table = F, stringsAsFactors = F)
genotypes= fread("resources/genotypes.vcf", sep="\t", header=T,  data.table = F, stringsAsFactors = F, na.strings=getOption("datatable.na.strings","."))

#2

#compute callrate for all SNP
callrate = c()
for(i in  1:nrow(genotypes)){
  callrate = c(callrate, 1 - (sum(is.na(genotypes[i,10:ncol(genotypes)])))/(ncol(genotypes)-9))
}

#add in a dataframe
df <- data.frame(callrate)

#create the histogram
png(filename = "callrate_histogram.png")
histograme <- ggplot(df, aes(x=callrate)) + geom_histogram(bins = 30) + ggtitle("callrate histogram")
histograme
dev.off()

#filtred the SNP in function of callrate
filt_genotype = genotypes[rowSums(is.na(genotypes[,10:ncol(genotypes)])) == 0,]

callratefiltred = c()
for(i in  1:nrow(filt_genotype)){
  callratefiltred = c(callratefiltred, 1 - (sum(is.na(filt_genotype[i,10:ncol(filt_genotype)])))/(ncol(filt_genotype)-9))
}
callratefiltred

supprimed = nrow(genotypes) - nrow(filt_genotype) # nombre de SNP retirés
supprimed
#3
altered = function(x) {
  vect = c()
  for(i in 1:ncol(x)) {
    if(x[i] == "1/1") {
      vect = c(vect,2)
    }
    else if(x[i] == "0/1" || x[i] == "1/0") {
      vect = c(vect,1)
    } else {
      vect = c(vect,0)
    }
  }
  return(vect)
}

MAF = c()
for(i in  1:nrow(filt_genotype)){
  AF = sum(altered(filt_genotype[i,10:ncol(filt_genotype)]))/((ncol(filt_genotype)-9)*2)
  if(AF <= 0.5){
    MAF = c(MAF,  AF)
  } else {
    MAF = c(MAF,  1-AF)
  }
}

df <- data.frame(MAF)

png(filename = "MAF histogram.png")
histograme_maf <- ggplot(df, aes(x=MAF)) + geom_histogram(bins = 30) + ggtitle("MAF histogram")
histograme_maf
dev.off()


final_genotype = filt_genotype[MAF >= 0.01,]

supprimed2 = nrow(filt_genotype) - nrow(final_genotype) # nombre de SNP retirés
supprimed2

#PART 3

#1

model = lm(phenotypes[,2] ~ covariates[,2])
#ou utiliser phenotypes$Cholesterol et same pour covariates
summary(model)
R = summary(model)$r.squared
R

#2

png(filename = "boxplot_cholesterol__gender.png")
gen_chol = merge(phenotypes,covariates)
gen_chol$gender = as.factor(gen_chol$gender)
head(gen_chol)

gcboxplot = ggplot(gen_chol, aes(x=gender, y=Cholesterol)) + 
  geom_boxplot() + ggtitle("boxplot of cholesterol level for different genders")

gcboxplot
dev.off()

png(filename = "density_cholesterol.png")

density_chol = gen_chol[,2:3]
density_chol$gender = as.factor(density_chol$gender)
density_plot = ggplot(density_chol, aes(x = Cholesterol, color = gender)) + geom_density() + ggtitle("distribution of cholesterol level for different genders")
density_plot
dev.off()

#3
data_for_pca = t(final_genotype[, 10:ncol(final_genotype)])
colnames(data_for_pca) = final_genotype$ID

# Change the string values in 3 numeric values 
# (Homozygous ref = 0, Heterozygous = 1, Homozygous alt = 2)
data_for_pca[data_for_pca == "0/0"] = 0
data_for_pca[data_for_pca == "1/0"] = 1
data_for_pca[data_for_pca == "0/1"] = 1
data_for_pca[data_for_pca == "1/1"] = 2

#The class is still string
str(data_for_pca)
storage.mode(data_for_pca) <- "numeric" # This change the storage mode
#cannot use as.numeric, that would transform the matrix in a vector
str(data_for_pca)

#Running PCA
# The data is ready, run the PCA (can take few seconds)
data_pca = prcomp(data_for_pca, center = T)

data_plot = data.frame(data_pca$x)

gender = as.factor(covariates$gender)

png(filename = "PCA.png")
p <- ggplot(data_plot, aes(x = PC1, y = PC2)) +
  geom_point(aes(color=gender)) + ggtitle("PCA performed on genotype data")
p
dev.off()

#4

p_val = c()
beta = c()
for (i in  1:ncol(data_for_pca)) {
  model_vp = lm(phenotypes[,2] ~ data_for_pca[,i])
  beta = c(beta,summary(model_vp)$coefficients[2,1])
  p_val = c(p_val,summary(model_vp)$coefficients[2,4])
}

#5 

position_manhattan = c()
Xaxis = c()
CHR = c()
current_chr = 0
relative_pos = 0
pos_before = 0
for (i in 1:nrow(final_genotype)) {
  if (current_chr != final_genotype[i,1]) {
    if(current_chr != 0) {
      Xaxis = c(Xaxis, (relative_pos)/2 + pos_before)
    }
    pos_before = relative_pos + pos_before
    relative_pos = final_genotype[i,2]
    CHR = c(CHR,final_genotype[i,1])
    current_chr = final_genotype[i,1]
  } else {
    relative_pos = final_genotype[i,2]-final_genotype[i-1,2] + relative_pos
  }
  position_manhattan = c(position_manhattan,final_genotype[i,2] +pos_before)
  }
Xaxis = c(Xaxis, (relative_pos)/2 + pos_before)


for_manhattan = data.frame("P" = p_val, "ID" = final_genotype$ID, "Chr" = final_genotype$`#CHROM`,"Position" = final_genotype$POS, "Chromosome_position" = position_manhattan)

sig_data = for_manhattan[for_manhattan$P < 0.05/length(p_val),]

Manhattan <- ggplot(for_manhattan, aes(x=Chromosome_position, y=-log10(P)))+ ggtitle("Manhattan plot without covariates") +
# Show all points
geom_point( aes(color=as.factor(Chr)), alpha=0.8, size=1.3) +
scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) + 
geom_hline(yintercept = -log10(0.05/length(p_val))) +
  
  # custom X axis:
  scale_x_continuous( label = CHR, breaks= Xaxis) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  
#Custom the theme:
theme_bw() +
theme( 
  legend.position="none",
  panel.border = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank()
) + geom_point(data =sig_data, aes(x=Chromosome_position, y=-log10(P)), color="orange", size=2) # Add highlighted points


png(filename = "Manhattan.png")
Manhattan 
dev.off()

#6 

p_val2 = c()
beta2 = c()
for (i in  1:ncol(data_for_pca)) {
  model_vp2 = lm(phenotypes[,2] ~ data_for_pca[,i]+ data_plot$PC1 + data_plot$PC2 + data_plot$PC3 + data_plot$PC4 + data_plot$PC5 + data_plot$PC6 + data_plot$PC7 + data_plot$PC8+ data_plot$PC9+ data_plot$PC10)
  beta2 = c(beta2,summary(model_vp2)$coefficients[2,1])
  p_val2 = c(p_val2,summary(model_vp2)$coefficients[2,4])
}


for_manhattan2 = data.frame("P" = p_val2,"ID" = final_genotype$ID, "Chr" = final_genotype$`#CHROM`,"Position" = final_genotype$POS, "Chromosome_position" = position_manhattan)

sig_data2 = for_manhattan2[for_manhattan2$P < 0.05/length(p_val2),]

Manhattan2 <- ggplot(for_manhattan2, aes(x=Chromosome_position, y=-log10(P)))+ ggtitle("Manhattan plot with covariates")  +
  
  # Show all points
  geom_point( aes(color=as.factor(Chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) + 
  geom_hline(yintercept = -log10(0.05/length(p_val2))) +
  
  # custom X axis:
  scale_x_continuous( label = CHR, breaks= Xaxis) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  
  #Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + geom_point(data =sig_data2, aes(x=Chromosome_position, y=-log10(P)), color="orange", size=2) # Add highlighted points


png(filename = "Manhattan with covariates.png")
Manhattan2 
dev.off()

#7 in the report 

#8

qqpoints_exp = ppoints(length(p_val2))
qqpoints_exp = sort(-log10(qqpoints_exp))

qqpoints_obs = sort(-log10(p_val2))

qqpoints_all = data.frame("Exp" = qqpoints_exp, "Obs" = qqpoints_obs)


png(filename = "Q-Q plot of GWAS p-values.png")
qqplot <- ggplot(qqpoints_all, aes(x=Exp,y=Obs)) + geom_point() + geom_abline(colour = "red") + ggtitle("Q-Q plot of GWAS p-values (-log10(p))")

qqplot

dev.off()

