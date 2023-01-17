'''
Author: Shloka Negi, shnegi@ucsc.edu
Usage: Run line by line in RStudio
Purpose: Plotting QQ value v/s F1 score curve combined for 2 modes (GiraffePE and BWA-MEM), separated for SNPs and INDELs
Input files required: extended.csv, roc.all.csv.gz, summary.csv (for Agilent and Ref-Seq both)

NOTE : GiraffePE files are named as vg.rerun
       GiraffenoPE files are named as vg.andrew
       BWA-MEM files are named as bwa.andrew
'''

# Installing required packages and Importing required libraries
devtools::install_github("Illumina/happyR")
library(happyR)
library(ggplot2)
library(magrittr)
theme_set(theme_minimal())

# Initialize and set working directory
Rplots <- "~/Desktop/patenlab_rotation/scripts/Rplots"
setwd(Rplots)

############## QQ v/s F1_score Plots ##############
######### 1. AGILENT-capture specific region ###########
giraffe_PE <- happyR::read_happy("agilent/vg.rerun_truth")
giraffe_PE <- pr_data(giraffe_PE, filter = "PASS") # 204
sum(giraffe_PE$Type == "INDEL") # 113
sum(giraffe_PE$Type == "SNP") # 149
giraffe_PE <- cbind("Mode"=rep("giraffe-PE", nrow(giraffe_PE)), giraffe_PE)

bwa_mem <- happyR::read_happy("agilent/bwa.andrew_truth")
bwa_mem <- pr_data(bwa_mem, filter = "PASS") # 204
sum(bwa_mem$Type == "INDEL") # 87
sum(bwa_mem$Type == "SNP") # 117
bwa_mem <- cbind("Mode"=rep("BWA-MEM", nrow(bwa_mem)), bwa_mem)

# Combine giraffe_PE and bwa_mem dfs
both_pr <- rbind(giraffe_PE, bwa_mem)

# Plotting
ggplot(both_pr, aes(x = QQ, y = METRIC.F1_Score, col= Mode)) +
  geom_line(size = 1) + theme_minimal() +
  scale_color_manual(values=c("blue", "red")) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "QQ_value v/s F1_Score Curve", subtitle = "(For Agilent capture-specific region)") +
  theme(plot.title = element_text(face="bold", size = 20, family = "Times New Roman")) + 
  theme(plot.subtitle = element_text(face = "bold", size = 16, family = "Times New Roman")) +
  theme(axis.title = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(face = "bold", size = 18, family = "Times New Roman")) +
  theme(legend.text = element_text(size = 12, face = "bold", family = "Times New Roman")) +
  theme(legend.background = element_rect(size = 0.5, linetype = "solid", colour = "white")) +
  guides(color = guide_legend(override.aes = list(linetype = 1))) +
  facet_wrap(~Type, scales = "free", ncol = 2, strip.position = "top") +
  theme(strip.background = element_rect(fill = "grey22")) +
  theme(panel.background = element_rect(fill = "seashell")) +
  theme(panel.border =  element_rect(color="black", fill = NA)) +
  theme(strip.text = element_text(face = "bold", size = 12, family = "Times", colour = "white"))

######### 2. REF-SEQ region ###########
giraffe_PE <- happyR::read_happy("ref-seq/vg.rerun_truth")
giraffe_PE <- pr_data(giraffe_PE, filter = "PASS") # 271
sum(giraffe_PE$Type == "INDEL") # 122
sum(giraffe_PE$Type == "SNP") # 149
giraffe_PE <- cbind("Mode"=rep("giraffe-PE", nrow(giraffe_PE)), giraffe_PE)

bwa_mem <- happyR::read_happy("ref-seq/bwa.andrew_truth")
bwa_mem <- pr_data(bwa_mem, filter = "PASS") # 223
sum(bwa_mem$Type == "INDEL") # 106
sum(bwa_mem$Type == "SNP") # 117
bwa_mem <- cbind("Mode"=rep("BWA-MEM", nrow(bwa_mem)), bwa_mem)

# Combine giraffe_PE and bwa_mem dfs
both_pr <- rbind(giraffe_PE, bwa_mem)

# Plotting
ggplot(both_pr, aes(x = QQ, y = METRIC.F1_Score, col= Mode)) +
  geom_line(size = 1) + theme_minimal() +
  scale_color_manual(values=c("blue", "red")) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "QQ_value v/s F1_Score Curve", subtitle = "(For general Ref-Seq region)") +
  theme(plot.title = element_text(face="bold", size = 20, family = "Times New Roman")) + 
  theme(plot.subtitle = element_text(face = "bold", size = 16, family = "Times New Roman")) +
  theme(axis.title = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(face = "bold", size = 18, family = "Times New Roman")) +
  theme(legend.text = element_text(size = 12, face = "bold", family = "Times New Roman")) +
  theme(legend.background = element_rect(size = 0.5, linetype = "solid", colour = "white")) +
  guides(color = guide_legend(override.aes = list(linetype = 1))) +
  facet_wrap(~Type, scales = "free", ncol = 2, strip.position = "top") +
  theme(strip.background = element_rect(fill = "grey22")) +
  theme(panel.background = element_rect(fill = "seashell")) +
  theme(panel.border =  element_rect(color="black", fill = NA)) +
  theme(strip.text = element_text(face = "bold", size = 12, family = "Times", colour = "white"))
