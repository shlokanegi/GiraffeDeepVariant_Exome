'''
Author: Shloka Negi, shnegi@ucsc.edu
Usage: Execute line by line, or blocks of code as required
Purpose: Analysis of Benchmarking results, such as finding intersecting Variants between Giraffe and BWA, etc.
Input file requirements: output VCFs from hap.py comparison
Conda environment: /public/home/shnegi/.conda/envs/giraffe
'''

### Working dir
SHNEGI="/public/groups/cgl/graph-genomes/shnegi"
### Intitialize directories
HAPPY="/public/home/anovak/build/hap.py/install"
INPUTS="/public/groups/cgl/graph-genomes/shnegi/inputs"
ref_genomes="/public/groups/cgl/graph-genomes/shnegi/ref-genomes"
SDK="/public/groups/cgl/graph-genomes/shnegi/google-cloud-sdk/bin"
SAMTOOLS="/public/groups/cgl/graph-genomes/shnegi/samtools-1.16.1"

cd ${SHNEGI}/happyout/vg_rerun/

#####1. Extracting the TP variants
awk '($10 ~ /TP/ && $11 ~ /TP/)' giraffePE_truth.vcf > giraffePE_truth_TP_TP.vcf
awk '($10 ~ /TP/ && $11 ~ /TP/)' giraffenoPE_truth.vcf > giraffenoPE_truth_TP_TP.vcf
awk '($10 ~ /TP/ && $11 ~ /TP/)' bwa_truth.vcf > bwa_truth_TP_TP.vcf

#####2. Add headers to the TP_TP VCFs (bedtools and UCSC-GB requires VCF files with header)
grep "^#" giraffePE_truth.vcf > header.txt; cat giraffePE_truth_TP_TP.vcf >> header.txt; mv header.txt giraffePE_truth_TP_TP.vcf
grep "^#" giraffenoPE_truth.vcf > header.txt; cat giraffenoPE_truth_TP_TP.vcf >> header.txt; mv header.txt giraffenoPE_truth_TP_TP.vcf
grep "^#" bwa_truth.vcf > header.txt; cat bwa_truth_TP_TP.vcf >> header.txt; mv header.txt bwa_truth_TP_TP.vcf

#####3. ANALYSIS - Comparing OLD giraffe TP_TP VCF with NEW Giraffe TP_TP VCF
# NEW GIRAFFE variants only
bedtools subtract -a giraffePE_truth_TP_TP.vcf -b giraffenoPE_truth_TP_TP.vcf > giraffePE_truth_TP_TP_only.vcf
grep "^#" giraffePE_truth.vcf > header.txt; cat giraffePE_truth_TP_TP_only.vcf >> header.txt; mv header.txt giraffePE_truth_TP_TP_only.vcf
# OLD GIRAFFE variants only
bedtools subtract -a giraffenoPE_truth_TP_TP.vcf -b giraffePE_truth_TP_TP.vcf > giraffenoPE_truth_TP_TP_only.vcf
grep "^#" giraffenoPE_truth.vcf > header.txt; cat giraffenoPE_truth_TP_TP_only.vcf >> header.txt; mv header.txt giraffenoPE_truth_TP_TP_only.vcf

#####4. Find BWA TP variants which were neither detected by OLD giraffe nor NEW giraffe
bedtools subtract -a bwa_truth_TP_TP.vcf -b giraffenoPE_truth_TP_TP.vcf > bwa_only_w_OLDgiraffe.vcf
grep "^#" bwa_truth_TP_TP.vcf > header.txt; cat bwa_only_w_OLDgiraffe.vcf >> header.txt; mv header.txt bwa_only_w_OLDgiraffe.vcf
bedtools subtract -a bwa_only_w_OLDgiraffe.vcf -b giraffePE_truth_TP_TP.vcf > bwa_only_w_OLDgiraffe_NEWgiraffe.vcf
grep "^#" bwa_only_w_OLDgiraffe.vcf > header.txt; cat bwa_only_w_OLDgiraffe_NEWgiraffe.vcf >> header.txt; mv header.txt bwa_only_w_OLDgiraffe_NEWgiraffe.vcf
# 768 BWA TP variants were still un-detected (VCF)

#####5. Convert to BED (required later for subsetting VCF based on BED regions)
vcf2bed --deletions <bwa_only_w_OLDgiraffe_NEWgiraffe.vcf> bwa_only_w_OLDgiraffe_NEWgiraffe_deletions.bed
vcf2bed --insertions <bwa_only_w_OLDgiraffe_NEWgiraffe.vcf> bwa_only_w_OLDgiraffe_NEWgiraffe_insertions.bed
vcf2bed --snvs <bwa_only_w_OLDgiraffe_NEWgiraffe.vcf> bwa_only_w_OLDgiraffe_NEWgiraffe_snvs.bed
bedops --everything bwa_only_w_OLDgiraffe_NEWgiraffe_{deletions,snvs,insertions}.bed > bwa_only_w_OLDgiraffe_NEWgiraffe_allvariants.bed ; rm bwa_only_w_OLDgiraffe_NEWgiraffe_deletions.bed bwa_only_w_OLDgiraffe_NEWgiraffe_snvs.bed bwa_only_w_OLDgiraffe_NEWgiraffe_insertions.bed
# 774 BWA TP variants according to the BED file

#####6. Subsetting the NEW giraffe "original" VCF at these BED regions to see why these variants were missed.
bcftools filter vg.HG003.sorted.vcf.gz -R bwa_only_w_OLDgiraffe_NEWgiraffe_allvariants.bed > vg.HG003.sorted.undetected.vcf
# 647 in VCF, 672 in BED

#####7. Now get those BWA-TP only variants which we couldn't find in the NEW giraffe "original" VCF
bedtools subtract -a bwa_only_w_OLDgiraffe_NEWgiraffe_allvariants.bed -b vg.HG003.sorted.undetected.bed > out.vcf

