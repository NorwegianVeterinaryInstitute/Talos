#!/usr/bin/env Rscript

# a script to process the output from hulk smash
# And vizualize them as a heatmap and a ordination plot

#functions
# check.packages function: install and load multiple R packages. 
# (downloaded from: https://gist.github.com/smithdanielle/9913897)
# Check to see if packages are installed. Install them if they are not, then load them into the R session.
# NOTE CRAN MIRROR SELECT IS norway MIRROR !!!
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos="https://cran.uib.no/")
  sapply(pkg, require, character.only = TRUE)
}

#Libraries to load and install if needed.
packages<-c("heatmap3")
check.packages(packages)

# loading and processing the data
all_samples <- read.csv("all_samples.Weighted_Jaccard.hulk-matrix.csv", head=TRUE)
colnames(all_samples) <- sub(".clean.fq.gz.json","",colnames(all_samples))
rownames(all_samples) <- colnames(all_samples)

samples <- ncol(all_samples)

#transforming the data to a matrix
all_samples_matrix <- as.matrix(all_samples)

#creating a heatmap using the Jaccard distances, and saving it.
colors <- colorRampPalette(c("steelblue", "white","orange2"))(100)

pdf(file= "Heatmap.all_samples.WeightedJaccard.pdf",
    width = 12, height = 9, bg="white") 

heatmap3(all_samples_matrix, scale="none", col=colors,cexRow=(0.8), cexCol=0.8, margins=c(10,6))
dev.off()