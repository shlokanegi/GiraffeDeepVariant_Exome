'''
Author: Shloka Negi, shnegi@ucsc.edu
Usage: Execute whatever is required
Purpose: Miscellaneous code snippets for mentioned tasks.
Conda environments: /public/home/shnegi/.conda/envs/giraffe
'''

#1. Make chromosome identifiers consistent for Giraffe and Truth
# Giraffe says GRCh38.chr1 and Truth says chr1. So, making Giraffe also chr1
sed -i 's/GRCh38.//' HG003.agilent.abra.vg.vg_exome_model.vcf

#2. Sort a BED file by chromosome
sort -V -k1,1 -k2,2 input.bed > output.bed

#3. Sort a VCF file by chromosome and position
bcftools sort input.vcf > output.vcf

#4. One way to subset vcf file based on given BED regions
vcftools --vcf vg_truth_TP.vcf --bed vg_truth_TP_only.bed --out vg_truth_TP_only --recode --keep-INFO-all
vcftools --vcf bwa_truth_TP.vcf --bed bwa_truth_TP_only.bed --out bwa_truth_TP_only --recode --keep-INFO-all
# We get "bwa_truth_TP_only.recode.vcf" & "vg_truth_TP_only.recode.vcf" files
# But this is not working on UCSC-GB as it needs VCF with headers

#5. How to read lines without the header
zcat vg_truth_TP_only_header.vcf.gz | grep -v "^#" | wc -l



