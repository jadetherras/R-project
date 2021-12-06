#PART 1

install.packages("data.table")
install.packages("ggplot2")

library(data.table)
library("ggplot2")

setwd("/Users/Djay/Desktop/BA5/G&g/R-project")

#PART 2

#1
covariates = fread("resources/covariates.txt", sep=",", header=T, data.table = F, stringsAsFactors = F)
phenotypes = fread("resources/phenotypes.txt", sep=",", header=T, data.table = F, stringsAsFactors = F)
genotypes = fread("resources/genotypes.vcf", sep="\t", header=T,  data.table = F, stringsAsFactors = F, na.strings=getOption("datatable.na.strings","."))

#2
callrate = c()
for(i in  1:nrow(genotypes)){
  callrate = c(callrate, 1 - (sum(is.na(genotypes[i,10:ncol(genotypes)])))/(ncol(genotypes)-9))
}

df <- data.frame(callrate)

png(filename = "callrate histogram")
histograme <- ggplot(df, aes(x=callrate)) + geom_histogram(bins = 30) + ggtitle("callrate histogram")
histograme
dev.off()

filt_genotype = genotypes[rowSums(is.na(genotypes[,10:ncol(genotypes)])) == 0,]

callrate_filtred = c()
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
  MAF = c(MAF,  (sum(altered(filt_genotype[i,10:ncol(filt_genotype)]))/((ncol(filt_genotype)-9)*2)))
}
New_MAF = 1-MAF

df <- data.frame(New_MAF)

png(filename = "MAF histogram")
histograme_maf <- ggplot(df, aes(x=New_MAF)) + geom_histogram(bins = 30) + ggtitle("MAF histogram")
histograme_maf
dev.off()

final_genotype = filt_genotype[New_MAF[] > 0.1,]

supprimed2 = nrow(filt_genotype) - nrow(final_genotype) # nombre de SNP retirés

#PART 3

#1

model = lm(phenotypes.Cholesterol, covariates.gender)
summary(model)

#2
Cholesterol_F = 
Cholesterol_M =
geom_boxplot()
geom_density()


#3
#4
#5
#6
#7
#8

