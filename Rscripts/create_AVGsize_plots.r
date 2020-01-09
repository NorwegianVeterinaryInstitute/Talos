#!/usr/bin/env Rscript

# a script to process the output from microbecensus
# And vizualize the output

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
packages<-c("reshape", "ggplot2")
check.packages(packages)

# create list of npo files to process with this script, and count element of list
avg_files <- list.files(".", pattern = "*.txt",full.names = TRUE)

number <- length(avg_files)

# Create a list of colors to use.
all_colors <- colorRampPalette(c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                                 "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
                                 "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"))

# making an big datafram with results of all data

datalist = list()

for (i in 1:number) {
      # ... make some data
      dat <- read.delim(avg_files[i], skip=10, sep=":", head=FALSE)
      dat2 <- data.frame(dat$V2)
      rownames(dat2) <- dat$V1
      name <- sub(".avgs_estimate.txt", "",avg_files[i]) #getting the name of the file and add to dataframe
      colnames(dat2) <- sub("./", "", name)# maybe you want to keep track of which iteration produced it?
      datalist[[i]] <- dat2 # add it to your list
    }
    
avg_data = do.call(cbind, datalist)
avg_data <- data.frame(t(avg_data))  # to get the datafram in a proper orientation.
avg_data$samples <- rownames(avg_data)

#melting dataframe
avg_data_melted <- melt(avg_data, id=c("samples"))


#plotting average genome sizes
pdf(file= "Average_genome_size.all_samples.pdf",
    width = 12, height = 9, bg="white") 

avg_size <- subset(avg_data_melted, variable == "average_genome_size")
size <- round_any((max(avg_size$value)), 1000000, ceiling)  ## determine the max value for settings the y-scale
ggplot(avg_size, aes(x=samples, y=(avg_size$value/1000000)))+
  geom_point(aes(color=samples)) +
  xlab("Samples") +
  ylab("Average genome size (Mbp)")+
  ggtitle("Average genome sizes of metagenomes")+
  ylim(0,(size/1000000))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # angled to fit more labels

dev.off()

#plotting genome equivalents
pdf(file= "Genome equivalents.all_samples.pdf",
    width = 12, height = 9, bg="white") 

genome_eq <- subset(avg_data_melted, variable == "genome_equivalents")
size <- round_any((max(genome_eq$value)), 10)  ## determine the max value for settings the y-scale
ggplot(genome_eq, aes(x=samples, y=value)) +
  geom_point(aes(color=samples)) +
  xlab("Samples") +
  ylab("Genome Equivalents")+
  ggtitle("Genome equivalents of metagenomes")+
  ylim(0,size)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # angled to fit more labels

dev.off()
