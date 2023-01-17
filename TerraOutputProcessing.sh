'''
Author: Shloka Negi, shnegi@ucsc.edu
Usage: Execute line by line, or blocks of code as required
Purpose: Post-processing the Terra output BAM and VCF files, before using for IGV visualizations and hap.py comparisons
TERRA specifications:
    Workspace: https://app.terra.bio/#workspaces/firecloud-cgl/hprc-evaluation
    Google Bucket: https://console.cloud.google.com/storage/browser/fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3?authuser=shnegi@ucsc.edu
    Bucket Name: fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3
    Workflow Used: GiraffeDeepVariant, based on WDL - https://github.com/vgteam/vg_wdl/blob/f2af46c56ea4aef50e6b437e93586ed0de59620e/workflows/giraffe_and_deepvariant.wdl
    Inputs JSON: provided on GitHub repo

'''

### Working dir
SHNEGI="/public/groups/cgl/graph-genomes/shnegi"
### Intitialize directories
INPUTS="/public/groups/cgl/graph-genomes/shnegi/inputs"
ref_genomes="/public/groups/cgl/graph-genomes/shnegi/ref-genomes"
SDK="/public/groups/cgl/graph-genomes/shnegi/google-cloud-sdk/bin"
SAMTOOLS="/public/groups/cgl/graph-genomes/shnegi/samtools-1.16.1"


############# NEW GIRAFFE VCF AND BAMs generated (outputs of Terra workflow) after providing the pairing information ################
#1. Get all inputs into the ${INPUTS}/HG003.agilent.PE directory
cd ${INPUTS}/HG003.agilent.PE
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-concatClippedVCFChunks/vg.HG003.vcf.gz .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-concatClippedVCFChunks/vg.HG003.vcf.gz.tbi .

$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-0/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr1.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-1/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr10.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-2/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr11.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-3/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr12.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-4/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr13.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-5/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr14.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-6/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr15.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-7/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr16.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-8/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr17.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-9/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr18.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-10/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr19.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-11/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr2.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-12/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr20.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-13/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr21.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-14/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr22.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-15/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr3.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-16/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr4.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-17/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr5.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-18/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr6.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-19/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr7.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-20/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr8.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-21/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chr9.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-22/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chrX.indel_realigned.bam .
$SDK/gsutil cp gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/submissions/e3b7bbbb-8e70-49a7-a182-a02519d1e754/GiraffeDeepVariant/002f5721-0a20-429f-9747-5374be1b9a57/call-runAbraRealigner/shard-23/cacheCopy/glob-f0d62d34bedc83144872e9cd1657dafc/vg.HG003.GRCh38.chrY.indel_realigned.bam .

#2. Merge all BAMs
samtools merge vg.HG003.GRCh38.merged.bam *.bam
#3. Sort the merged BAM
samtools sort -o vg.HG003.GRCh38.merged.sorted.bam vg.HG003.GRCh38.merged.bam

#4. Sort the VCF and change the chromosome names to a consistent one
grep "^#" vg.HG003.vcf > vg.HG003.sorted.vcf
grep -v "^#" vg.HG003.vcf| sort -k1,1V -k2,2g >> vg.HG003.sorted.vcf
sed -i 's/GRCh38.//' vg.HG003.sorted.vcf