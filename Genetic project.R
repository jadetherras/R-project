#PART 1

install.packages("data.table")
install.packages("ggplot2")

library(data.table)
library("ggplot2")

setwd("/Users/Djay/Desktop/BA5/G&g/R-project")

#PART 2

#1
covariates <- fread("resources/covariates.txt", sep=",", header=T, data.table = F, stringsAsFactors = F)
phenotypes = fread("resources/phenotypes.txt", sep=",", header=T, data.table = F, stringsAsFactors = F)
genotypes = fread("resources/genotypes.vcf", sep="\t", header=T,  data.table = F, stringsAsFactors = F, na.strings=getOption("datatable.na.strings","."))

#2
callrate = c()
for(i in  1:nrow(genotypes)){
  callrate = c(callrate, 1 - (sum(is.na(genotypes[i,10:ncol(genotypes)])))/(ncol(genotypes)-9))
}

df <- data.frame(callrate)

png(filename = "callrate_histogram.png")
histograme <- ggplot(df, aes(x=callrate)) + geom_histogram(bins = 30) + ggtitle("callrate histogram")
histograme
dev.off()

filt_genotype = genotypes[rowSums(is.na(genotypes[,10:ncol(genotypes)])) == 0,]

callratefiltred = c()
for(i in  1:nrow(filt_genotype)){
  callratefiltred = c(callratefiltred, 1 - (sum(is.na(filt_genotype[i,10:ncol(filt_genotype)])))/(ncol(filt_genotype)-9))
}

supprimed = nrow(genotypes) - nrow(filt_genotype) # nombre de SNP retirés

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
#demander encore au prof s'il faut tester pour trouver le maf (ou si on peut juste partir du principe que c'est bon)

#PART 3

#1

model = lm(phenotypes[,2] ~ covariates[,2])
#ou utiliser phenotypes$Cholesterol et same pour covariates

R = summary(model)$r.squared

#R^2 = -0.0002128, which indicates a high probability that the gender does impact the cholesterol level

#2

png(filename = "boxplot_cholesterol__gender.png")
gen_chol = merge(phenotypes,covariates)
gen_chol$gender = as.factor(gen_chol$gender)
head(gen_chol)

gcboxplot = ggplot(gen_chol, aes(x=gender, y=Cholesterol)) + 
  geom_boxplot()

gcboxplot
dev.off()

png(filename = "density_cholesterol.png")

density_chol = gen_chol[,2:3]
density_chol$gender = as.factor(density_chol$gender)
density_plot = ggplot(density_chol, aes(x = Cholesterol, color = gender)) + geom_density()
density_plot
dev.off()

#As there is a shift in the densities between the two genders, we can conclude that the gender is a covariate of cholesterol


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

png(filename = "PCA.png")
p <- ggplot(data_plot, aes(x = PC1, y = PC2)) +
  geom_point()
p
dev.off()



# We see 3 clusters. Each cluster represent one population cluster that have preferencialy specific SNP.
#no cause 3 ethnicies, 3 clusters

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
pos = 0
for (i in 1:nrow(final_genotype)) {
  pos = pos + final_genotype[i,2]
  position_manhattan = c(position_manhattan,pos)
}

for_manhattan = data.frame("P" = p_val, "Chr" = final_genotype$`#CHROM`,"Position" = final_genotype$POS, "Chromosome_position" = position_manhattan)

Manhattan <- ggplot(for_manhattan, aes(x=Chromosome_position, y=-log10(P))) +
  
# Show all points
geom_point( aes(color=as.factor(Chr)), alpha=0.8, size=1.3) +
scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) + 
geom_hline(yintercept = -log10(0.05/length(p_val))) +
  
#Custom the theme:
theme_bw() +
theme( 
  legend.position="none",
  panel.border = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank()
)

png(filename = "Manhattan.png")
Manhattan #BONUS plot colors
dev.off()

#6 

p_val2 = c()
beta2 = c()
for (i in  1:ncol(data_for_pca)) {
  model_vp2 = lm(phenotypes[,2] ~ data_for_pca[,i]+ data_plot$PC1 + data_plot$PC2 + data_plot$PC3 + data_plot$PC4 + data_plot$PC5 + data_plot$PC6 + data_plot$PC7 + data_plot$PC8+ data_plot$PC9+ data_plot$PC10)
  beta2 = c(beta2,summary(model_vp2)$coefficients[2,1])
  p_val2 = c(p_val2,summary(model_vp2)$coefficients[2,4])
}

position_manhattan2 = c()
pos2 = 0
for (i in 1:nrow(final_genotype)) {
  pos2 = pos2 + final_genotype[i,2]
  position_manhattan2 = c(position_manhattan2,pos2)
}

for_manhattan2 = data.frame("P" = p_val2, "Chr" = final_genotype$`#CHROM`,"Position" = final_genotype$POS, "Chromosome_position" = position_manhattan2)

Manhattan2 <- ggplot(for_manhattan2, aes(x=Chromosome_position, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(Chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) + 
  geom_hline(yintercept = -log10(0.05/length(p_val2))) +
  
  #Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

png(filename = "Manhattan with covariates.png")
Manhattan2 #BONUS plot colors
dev.off()



