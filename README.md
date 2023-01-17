# GoogleDeepVariant_Exome
Performance Evaluation of  Giraffe-DeepVariant Workflow on WES data

## Project Background
### Variation Graph and Pangenome
Variation graphs provide a succinct encoding of the sequences of many genomes. <br>
A Pangenome is a collection of genomes from individuals in a population and is represented as a variation graph. 

### Giraffe (Short-Read Mapping Algorithm)
[Giraffe](https://github.com/vgteam/vg) is a short-read mapping algorithm that finds the walk through the pangenomic graph that best matches the read.

<img width="790" alt="image" src="https://user-images.githubusercontent.com/66521525/212864627-6a62166c-8013-4c33-9273-c9b8b79200ca.png">

### Google’s DeepVariant (Variant Calling Algorithm)
Read alignments from Giraffe can be passed to [DeepVariant](https://google.github.io/deepvariant/), which visualizes these into pileup images (can be 6 or 7 channels deep) , and then uses a CNN to classify the images into variant calls.

<img width="827" alt="image" src="https://user-images.githubusercontent.com/66521525/212865383-301fe8a0-0ddd-4027-bcc4-c8a2a9a02b79.png">

# Giraffe-DV v/s BWA-DV Performance Statistics by Google

<img width="870" alt="image" src="https://user-images.githubusercontent.com/66521525/212865991-4ef25bbc-9f7c-4996-9da5-34ba17d33408.png">

# Methodology
1. Conversion of BAM to paired-end fastQ files. 
2. Attempted the [GiraffeDeepVariant](https://github.com/vgteam/vg_wdl/blob/f2af46c56ea4aef50e6b437e93586ed0de59620e/tasks/deepvariant.wdl) workflow on Terra.
3. Fixed the input parameter file to succeed.
4. Fetched Giraffe-PE VCF. Used it for [hap.py](https://github.com/Illumina/hap.py) benchmarking comparisons. 
5. Investigations of exclusive BWA-only variants on IGV.

# Results
## hap.py Benchmarking Comparison
<img width="824" alt="image" src="https://user-images.githubusercontent.com/66521525/212867347-3353c9a6-99a8-4455-88a2-eb7814d0d6b8.png">

## hap.py Benchmarking Results
<img width="803" alt="image" src="https://user-images.githubusercontent.com/66521525/212867480-a3802615-50a8-4417-b5e9-fd7cb4adabcd.png">

<img width="830" alt="image" src="https://user-images.githubusercontent.com/66521525/212867905-916d909f-0f31-4a74-a665-7972927d7120.png">

<img width="808" alt="Screenshot 2023-01-17 at 2 00 54 AM" src="https://user-images.githubusercontent.com/66521525/212868210-7d497191-6b6b-4cb7-a7fe-ceea574b7011.png">

<img width="784" alt="image" src="https://user-images.githubusercontent.com/66521525/212868398-a624a0b9-a689-453c-8fc0-5ca49824e8e3.png">

## IGV Visualizations
<img width="867" alt="image" src="https://user-images.githubusercontent.com/66521525/212868522-17591683-5182-4330-9534-ac49fbd669cb.png">

<img width="865" alt="image" src="https://user-images.githubusercontent.com/66521525/212868591-070aa35d-e074-49a5-96ab-2a7c3516892b.png">

# Conclusion and Future Work
* DeepVariant Model Retraining on correct paired-end exome data.
* Repeat analysis
* Evaluate performance of Giraffe on other kits, (Nextera and IDT) as well as the general ref-seq region.

