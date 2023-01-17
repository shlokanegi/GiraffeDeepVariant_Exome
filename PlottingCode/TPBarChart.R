'''
Author: Shloka Negi, shnegi@ucsc.edu
Usage: Run line by line in RStudio
Purpose: 
  1. Bar Charts for comparing True Positive Variants between BWA-MEM and Giraffe-PE, separated for SNPs and INDELs
  2. Bar Charts for showing performance evaluation report by Google.
Input files required: extended.csv, summary.csv (for Agilent and Ref-Seq both)

NOTE : GiraffePE files are named as vg.rerun
       GiraffenoPE files are named as vg.andrew
       BWA-MEM files are named as bwa.andrew

'''

# Installing required packages and Importing required libraries
devtools::install_github("Illumina/happyR")
library(happyR)
library(ggplot2)
theme_set(theme_minimal())
library(viridis)
library(hrbrthemes)

########## 1. Bar Charts for comparing True Positive Variants between BWA-MEM and Giraffe-PE, separated for SNPs and INDELs #################
# Initialize and set working directory
Rplots <- "~/Desktop/patenlab_rotation/scripts/Rplots"
setwd(Rplots)

# Prepare datasets for plotting boxplots
giraffe_PE <- happyR::read_happy("agilent/vg.rerun_truth")
giraffe_PE <- giraffe_PE$summary
giraffe_PE <- giraffe_PE[giraffe_PE$Filter=="PASS",]
giraffe_PE <- cbind("Mode"=rep("giraffe-PE", nrow(giraffe_PE)), giraffe_PE)
bwa_mem <- happyR::read_happy("agilent/bwa.andrew_truth")
bwa_mem <- bwa_mem$summary
bwa_mem <- bwa_mem[bwa_mem$Filter=="PASS",]
bwa_mem <- cbind("Mode"=rep("BWA-MEM", nrow(bwa_mem)), bwa_mem)
both_pr <- rbind(giraffe_PE, bwa_mem)

# Plotting Bar-chart
ggplot(both_pr, aes(fill=Mode, y=TRUTH.TP, x=Mode, label=TRUTH.TP, color=Mode)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(fontface="bold"), size = 5, vjust=-0.25) +
  theme_minimal() + 
  facet_wrap(~Type, scales = "free", ncol = 2, strip.position = "top") +
  scale_fill_manual(values=c("blue", "red")) +
  scale_color_manual(values=c("blue", "red")) +
  labs(x='Mode', y='No. of TP variants', title='Comparison of True Positive Variants', subtitle = '(Agilent capture-specific region)') +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size=18, face='bold', family = "Times")) +
  theme(plot.subtitle = element_text(size = 16, family = 'Times')) +
  theme(axis.title = element_text(size=14, face="bold")) +
  theme(axis.text = element_text(size = 10)) +
  theme(strip.background = element_rect(fill = "grey22")) +
  theme(panel.border =  element_rect(color="black", fill = NA)) +
  theme(strip.text = element_text(face = "bold", size = 12, colour = "white"))

  
########## 2. Prepare Barcharts for showing Google's performance evaluation report ###############
# Making a dataframe from scratch
df <- data.frame("Mode"=c(rep("Giraffe-NoPE",2), rep("BWA-MEM", 2)), 
                 "Type"=c("INDEL", "SNP", "INDEL", "SNP"),
                 "TRUTH.TP"=c(3870, 49804, 9742, 89434),
                 "METRIC.Recall"=c(0.2719, 0.4225, 0.4075, 0.5142),
                 "METRIC.Precision"=c(0.8084, 0.9313, 0.9007, 0.9545))

# Plotting
ggplot(df, aes(fill=Mode, y=TRUTH.TP, x=Mode, label=TRUTH.TP, color=Mode)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(fontface="bold"), size = 6, vjust=-0.25) +
  theme_minimal() + 
  facet_wrap(~Type, scales = "free", ncol = 2, strip.position = "top") +
  scale_fill_manual(values=c("blue", "green4")) +
  scale_color_manual(values=c("blue", "green4")) +
  labs(x='Mode', y='No. of TP variants', title='Comparison of True Positive Variants', subtitle = '(Agilent capture-specific region)') +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size=18, face='bold', family = "Times")) +
  theme(plot.subtitle = element_text(size = 16, family = 'Times')) +
  theme(axis.title = element_text(size=14, face="bold")) +
  theme(axis.text = element_text(size = 10)) +
  theme(strip.background = element_rect(fill = "grey22")) +
  theme(panel.border =  element_rect(color="black", fill = NA)) +
  theme(strip.text = element_text(face = "bold", size = 12, colour = "white"))
