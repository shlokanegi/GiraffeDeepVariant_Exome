'''
Author: Shloka Negi, shnegi@ucsc.edu
Usage: Execute line by line, or blocks of code as required
Purpose: Convert BAM file (missing pairing information) to paired-end fastQ files (suitable for upload on Terra)
Input file requirements: Unpaired BAM file
Conda environment: /public/home/shnegi/.conda/envs/giraffe

Agilent vg data (generated by google) downloaded from: 
    https://storage.googleapis.com/brain-genomics/awcarroll/share/ucsc/vg_exomes/vcf/HG003.agilent.abra.vg.vg_exome_model.vcf.gz
Google Bucket :
    fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/exomes/
'''

### Working dir
SHNEGI="/public/groups/cgl/graph-genomes/shnegi"
### Intitialize directories
INPUTS="/public/groups/cgl/graph-genomes/shnegi/inputs"
SDK="/public/groups/cgl/graph-genomes/shnegi/google-cloud-sdk/bin"
SAMTOOLS="/public/groups/cgl/graph-genomes/shnegi/samtools-1.16.1"

#1. Convert VG BAM to interleaved FASTQs
$SAMTOOLS/samtools view -O SAM ${INPUTS}/agilent/HG003.agilent.abra.vg.bam | cut -f1 | sort | uniq -c | grep -v '^ *1' | rev | cut -f1 | rev > ${SHNEGI}/fastq-ilv/paired-reads_vg.txt
$SAMTOOLS/samtools view -O BAM -N ${SHNEGI}/fastq-ilv/paired-reads_vg.txt ${INPUTS}/agilent/HG003.agilent.abra.vg.bam | $SAMTOOLS/samtools collate -O - > ${SHNEGI}/fastq-ilv/temp_vg.bam
~anovak/.local/bin/bamToFastq -i ${SHNEGI}/fastq-ilv/temp_vg.bam -fq ${SHNEGI}/fastq-ilv/vg.reads.interleaved.fq

#2. Convert interleaved fastQs to paired-end fastQs
paste - - - - - - - - < vg.reads.interleaved.fq \
    | tee >(cut -f 1-4 | tr "\t" "\n" > vg.read1.fq) \
    | cut -f 5-8 | tr "\t" "\n" > vg.read2.fq

#3. Convert them to zipped format (Giraffe-DV workflow only takes fastq.gz files)
gzip vg.read1.fq.gz
gzip vg.read2.fq.gz

#4. Tranfer the fastq.gz files to google bucket
$SDK/gsutil cp vg.read1.fq.gz gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/exomes/vg.read1.fq.gz
$SDK/gsutil cp vg.read2.fq.gz gs://fc-0e9b62da-2dd8-42e5-bb50-8d08e7cdebc3/exomes/vg.read2.fq.gz
