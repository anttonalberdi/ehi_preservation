
setwd("/Users/anttonalberdi/bamse_faeces_animals")

#Load libraries
library(hilldiv)

counts <- read.csv("ASV_counts.csv",row.names=1)
metadata <- read.csv("metadata.csv")

#Remove samples with <1000 reads
counts <- depth_filt(counts, 1000)

#Buffer test
divtest <- div_test(counts,qvalue=q,hierarchy=metadata[,c(1,2)])
div_test_plot(divtest)

#Species test
divtest <- div_test(counts,qvalue=1,hierarchy=metadata[,c(1,3)])
div_test_plot(divtest)

#Time test
divtest <- div_test(counts,qvalue=1,hierarchy=metadata[,c(1,4)])
div_test_plot(divtest)

#Protocol test
divtest <- div_test(counts,qvalue=1,hierarchy=metadata[,c(1,5)])
div_test_plot(divtest)
