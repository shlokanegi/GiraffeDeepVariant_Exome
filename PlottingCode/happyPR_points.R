'''
Author: Shloka Negi, shnegi@ucsc.edu
Usage: Run line by line in RStudio
Purpose: 
  1. Precision v/s Recall plot for all 3 modes (GiraffePE, GiraffenoPE and BWA-MEM)
  2. Precision v/s Recall plot to display performance evaluation report by Google.
Input files required: extended.csv, roc.all.csv.gz, summary.csv (for Agilent and Ref-Seq both)

NOTE : GiraffePE files are named as vg.rerun
       GiraffenoPE files are named as vg.andrew
       BWA-MEM files are named as bwa.andrew

Tool Documentation : https://github.com/Illumina/happyR 
'''

# Installing required packages and Importing required libraries
devtools::install_github("Illumina/happyR")
library(happyR)
library(tidyverse, quietly = TRUE)
require("ggrepel")


###########  1. Precision v/s Recall plot for all 3 modes (GiraffePE, GiraffenoPE and BWA-MEM) ##############
# Initialize and set working directory
Rplots <- "~/Desktop/patenlab_rotation/scripts/Rplots"
setwd(Rplots)

########## For ref-seq target region ###########
# Prepare dataset
samplesheet <- tibble::tribble(
  ~group_id,  ~replicate_id, ~happy_prefix,
  "giraffe-PE", "vg.rerun_truth", paste(Rplots, "ref-seq", "vg.rerun_truth", sep = "/"),
  "giraffe-noPE", "vg.andrew_truth", paste(Rplots, "ref-seq", "vg.andrew_truth", sep = "/"),
  "bwa-mem", "bwa.andrew_truth",  paste(Rplots, "ref-seq", "bwa.andrew_truth", sep = "/")
)
hap_samplesheet <- read_samplesheet_(samplesheet)
summary <- extract_results(hap_samplesheet$results, table = "summary") %>% 
  inner_join(samplesheet, by = "happy_prefix") %>% 
  filter(Filter == "PASS")

# Coordinates corresponding to Recall and Precision
R <- c(0.28, 0.44, 0.27, 0.42, 0.27, 0.42)
P <- c(0.75, 0.89, 0.80, 0.93, 0.85, 0.93)
coords = paste(R,P,sep=", ")
summary_coords <- cbind(coords, summary)

# Plotting
ggplot(data = summary_coords, aes(x = METRIC.Recall, y = METRIC.Precision, color = group_id, shape = Type)) +
  geom_point(size = 4) + theme_minimal() +
  geom_text_repel(aes(label = coords), size=4, box.padding = 1, fontface = "bold") +
  xlim(NA, 1) + ylim(NA, 1) +
  scale_color_manual(values=c("blue", "forestgreen", "red")) +
  labs(x = "Recall", y = "Precision", color = "mapping-algorithm") +
  ggtitle("Precision vs. Recall of mapping algorithms with DeepVariant workflow",
          "General ref-seq target region") +
  theme(plot.title = element_text(face="bold", size = 20, family = "Times New Roman")) + 
  theme(plot.subtitle = element_text(face = "bold", size = 16, family = "Times New Roman")) +
  theme(axis.title = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(face = "bold", size = 16, family = "Times New Roman")) +
  theme(legend.text = element_text(size = 12, family = "Times New Roman")) +
  theme(legend.background = element_rect(size = 0.5, linetype = "solid", colour = "white")) +
  guides(color = guide_legend(override.aes = list(linetype = 1)))


######## For agilent capture-specific target region #########
Rplots <- "~/Desktop/patenlab_rotation/scripts/Rplots"
setwd(Rplots)

# Prepare dataset
samplesheet <- tibble::tribble(
  ~group_id,  ~replicate_id, ~happy_prefix,
  "giraffe-PE", "vg.rerun_truth", paste(Rplots, "agilent", "vg.rerun_truth", sep = "/"),
  "giraffe-noPE", "vg.andrew_truth", paste(Rplots, "agilent", "vg.andrew_truth", sep = "/"),
  "bwa-mem", "bwa.andrew_truth",  paste(Rplots, "agilent", "bwa.andrew_truth", sep = "/")
)
hap_samplesheet <- read_samplesheet_(samplesheet)
summary <- extract_results(hap_samplesheet$results, table = "summary") %>% 
  inner_join(samplesheet, by = "happy_prefix") %>% 
  filter(Filter == "PASS")

# Coordinates corresponding to Recall and Precision
R <- c(0.93, 0.99, 0.80, 0.95, 0.81, 0.95)
P <- c(0.95, 0.991, 0.96, 0.997, 0.993, 0.997)
coords = paste(R,P,sep=", ")
summary_coords <- cbind(coords, summary)

# Plotting
ggplot(data = summary_coords, aes(x = METRIC.Recall, y = METRIC.Precision, color = group_id, shape = Type)) +
  geom_point(size = 4) + theme_minimal() +
  geom_text_repel(aes(label = coords), size=4, box.padding = 1, fontface = "bold") +
  xlim(NA, 1) + ylim(NA, 1) +
  scale_color_manual(values=c("blue", "forestgreen", "red")) +
  labs(x = "Recall", y = "Precision", color = "mapping-algorithm") +
  ggtitle("Precision vs. Recall of mapping algorithms with DeepVariant workflow",
          "Agilent capture-specific target region") +
  theme(plot.title = element_text(face="bold", size = 20, family = "Times New Roman")) + 
  theme(plot.subtitle = element_text(face = "bold", size = 16, family = "Times New Roman")) +
  theme(axis.title = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(face = "bold", size = 16, family = "Times New Roman")) +
  theme(legend.text = element_text(size = 12, family = "Times New Roman")) +
  theme(legend.background = element_rect(size = 0.5, linetype = "solid", colour = "white")) +
  guides(color = guide_legend(override.aes = list(linetype = 1)))


###########  2. Precision v/s Recall plot to display performance evaluation report by Google. ##############
# Making a dataframe from scratch
df <- data.frame("Mode"=c(rep("Giraffe-NoPE",2), rep("BWA-MEM", 2)), 
                 "Type"=c("INDEL", "SNP", "INDEL", "SNP"),
                 "TRUTH.TP"=c(3870, 49804, 9742, 89434),
                 "METRIC.Recall"=c(0.2719, 0.4225, 0.4075, 0.5142),
                 "METRIC.Precision"=c(0.8084, 0.9313, 0.9007, 0.9545))

# Defining coordinates
R <- c(0.27, 0.42, 0.41, 0.51)
P <- c(0.81, 0.93, 0.90, 0.95)
coords = paste(R,P,sep=", ")
df_coord = cbind(coords, df)

# Plotting
ggplot(data = df_coord, aes(x = METRIC.Recall, y = METRIC.Precision, color = Mode, shape = Type)) +
  geom_point(size = 6) + theme_minimal() +
  geom_text_repel(aes(label = coords), size=4, box.padding = 1, fontface = "bold") +
  xlim(NA, 1) + ylim(NA, 1) +
  scale_color_manual(values=c("blue", "forestgreen")) +
  labs(x = "Recall", y = "Precision", color = "mapping-algorithm") +
  ggtitle("Precision vs. Recall of mapping algorithms with DeepVariant workflow",
          "Agilent capture-specific target region") +
  theme(plot.title = element_text(face="bold", size = 20, family = "Times New Roman")) + 
  theme(plot.subtitle = element_text(face = "bold", size = 16, family = "Times New Roman")) +
  theme(axis.title = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(face = "bold", size = 16, family = "Times New Roman")) +
  theme(legend.text = element_text(size = 12, family = "Times New Roman")) +
  guides(color = guide_legend(override.aes = list(linetype = 1)))

