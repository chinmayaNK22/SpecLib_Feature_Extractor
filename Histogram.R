library(ggplot2)

setwd("D:\\OSCC_DIA/DIA_Proteomics_April_2021/OSCC_SpecLib_Features/")


data <- read.delim("OSCC_Crude_Chewers_Smokers_Multi-Consesus_Spec_LibSpecLib_peptide_fragment_frequency.txt", sep = '\t')

colnames(data)

ggplot(data, aes(No..of.Fragments)) +
  geom_histogram(binwidth = 1)

# For histograms with tick marks between each bin, use `geom_bar()` with
# `scale_x_binned()`.
ggplot(data, aes(Length)) + geom_bar() + scale_x_binned()
