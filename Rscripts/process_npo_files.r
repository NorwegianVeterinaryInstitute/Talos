#!/usr/bin/env Rscript

# a script to process the data from nonpareil and generate a meaningful figure

#functions
# check.packages function: install and load multiple R packages. 
# (downloaded from: https://gist.github.com/smithdanielle/9913897)
# Check to see if packages are installed. Install them if they are not, then load them into the R session.
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos="https://cran.uib.no/")
  sapply(pkg, require, character.only = TRUE)
}

#Libraries to load and install if needed.
packages<-c("Nonpareil")
check.packages(packages)

# create list of npo files to process with this script, and count element of list
npo_files <- list.files(".", pattern = "*.npo",full.names = TRUE)

number <- length(npo_files)

# Create a list of colors to use.
all_colors <- colorRampPalette(c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                                 "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
                                 "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"))

# plotting each sample as an individual graph with an dispersion distribution around it.
#system("mkdir single_plots")

for (i in 1:number){
  myfile <- file.path("single_plots", paste0(npo_files[i], ".np_plot.png"))
  png(filename= myfile,
      width = 12, height = 9, units="in", res=300, bg="white") 
  p <- Nonpareil.curve(npo_files[i], plot.dispersion="sd", col="black")
  Nonpareil.legend(p)
  dev.off()
}
## printing all samples in one graph
png(filename= "non_pareil_plot.png",
    width = 12, height = 9, units="in", res=300, bg="white") 
p <- Nonpareil.curve.batch(npo_files, col=all_colors(number))
dev.off() 
