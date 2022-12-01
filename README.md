# SNPcalling_2022.md
## Quality control

I started by creating folders called myanalyses and sourcefiles to store my ouput files.

```
mkdir myanalyses
mkdir sourcefiles
cd sourcefiles
```

I then subsetted the large files (.gz) to 250,000 bp reads (.fq) as this will be quicker to test the initial structure and adapter content before working on the whole data. As this is paired-end data I repeated this on read 1 and read 2 independently.

```
#!/bin/sh  
zcat  GK_GM_S1_R1_001.fastq.gz | head -n 1000000 > GK_GM_S1_R1_sample.fq 

#!/bin/sh  
zcat  GK_GM_S1_R2_001.fastq.gz | head -n 1000000 > GK_GM_S1_R2_sample.fq 
```

I then wanted to check the quality of the data using fastqc. 

```
#!/bin/sh 
module load FastQC 
fastqc *fq 
```

The fastqc output for read 1 and read 2 shows quite a lot of adapter contamination by an Illumina Universal Adapter. 

Looking at the Illumina support website, this adapter sequence is “AGATCGGAAGAG”. 

I then decided to retain sequences with a minimum (-m) of 30 base pairs in length from read 1 and read 2, using cutadapt on a subsample of the data. This acts to remove the last bases in each lane which contained a lot of adapter contamination. 

```
#!/bin/sh 
module load cutadapt 
cutadapt -j 4  -a AGATCGGAAGAG -A AGATCGGAAGAG -m 30:30 -o trimmed_GK_GM_S1_R1_sample.fq -p trimmed_GK_GM_S1_R2_sample.fq GK_GM_S1_R1_sample.fq GK_GM_S1_R2_sample.fq 
```

I then tested whether this was successful in removing adapter contamination by running fastqc on the output of cutadapt.

```
fastqc trimmed_GK_GM_S1_R1_sample.fq trimmed_GK_GM_S1_R2_sample.fq 
```

Output by fastqc for read 1 and read 2 show adapter contamination was successfully removed.

I then ran cutadapt on the whole data and checked this was successful in removing adapter contamination using fastqc.

```
cutadapt -j 4  -a AGATCGGAAGAG -A  AGATCGGAAGAG -m 30:30 -o trimmed_GK_GM_S1_R1_.fq -p trimmed_GK_GM_S1_R2_.fq GK_GM_S1_R1_001.fastq.gz GK_GM_S1_R2_001.fastq.gz
fastqc trimmed_GK_GM_S1_R1_.fq trimmed_GK_GM_S1_R2_.fq 
```

```
Read 1 
17,204,395,192 bp in 
13,788,154,797 bp out  
 
Read 2  
17,204,395,192 bp in 
13,836,853,326 bp out  
```

## Demultiplexing

I could then demultiplex my samples based on their barcode. I started by making two new folders; "raw" for inputting my files to be demultiplexed, and "samples" for the output of demultiplexing. I then created links to put my trimmed whole data for read 1 and read 2 into "raw".

```
#!/bin/sh 
mkdir raw samples 
cd raw 
ln -s ../trimmed_GK_GM_S1_R1_.fq trimmed_GK_GM_S1_R1_001.fastq 
ln -s ../trimmed_GK_GM_S1_R2_.fq trimmed_GK_GM_S1_R2_001.fastq 
cd ..  
```

Then I extracted the .key file containing the barcodes used for each of read 1 and read 2 of samples. I loaded this file into a text editor for formatting the data, according to section 1.4.2 of the Stacks manual. I named this file "barcodes.txt" and loaded it into the sourcefiles folder.

I then loaded Stacks to use the programme process_radtags for demultiplexing.

```
#!/bin/sh 
module load Stacks
```

For demultiplexing, I specified the data as paired end (-P), the barcodes.txt file for barcodes, the restriction enzyme used as PstI (-e)  and for the program to rescue barcodes and cut sites (-r). I also wanted to have higher quality data so I specified (-c) to remove any reads with uncalled bases and (-q) to discard reads with low quality scores.

I also specified the barcode type according to section 1.4.2 of the Stacks manual, as paired end with inline barcodes on the single and paired-ends (--inline_inline)

```
process_radtags –p raw/ -P -b barcodes.txt -o ./samples/ -e PstI –r –c –q --inline_inline 
```

This retains 90% of reads.

```
227773798 total sequences 
20262420 barcode not found drops (8.9%) 
    99519 low quality read drops (0.0%) 
  2344424 RAD cutsite not found drops (1.0%) 
205067435 retained reads (90.0%) 
```

## Alignment and variant calling

Sample files were then clean and able to be aligned to the reference genome.

The first step was to index the genome. To do this I used Burrows-Wheeler Aligner (BWA).

```
#!/bin/sh 
module load BWA 
bwa index stoneflygenomeassemblyv1.fasta 
```

I then created a new folder "alignment" and moved the indexed genome to this folder to conduct alignment.

```
mkdir alignment 
$ cd alignment 
cp ../^C 
cp ../stoneflygenomeassemblyv1.fasta* ./ 
```

The BWA-MEM algorithm and SAMtools were then used to align each of read 1 and read 2 for every sample to the reference stonefly genome using a loop command. I specified the number of CPU's the programme was running as 4 (-t). All commands can be found in align.sh.

```
module load SAMtools 
```

```
#!/bin/sh 
src=../samples/ 
bwa_db=stoneflygenomeassemblyv1.fasta 
for sample in $files 
do  
    echo $sample 
    bwa mem -t 4 $bwa_db $src/${sample}.1.fq.gz $src/${sample}.2.fq.gz |   samtools view -b | samtools sort --threads 4 > ${sample}.bam 
done 
```

I then wanted to run the Stacks programme refmap to filter out any low-quality individuals. To do this I created a folder "aligned_samples" and moved the BAM output files from reference alignment into this folder, to use as input for refmap. I also created the folder "output_refmap" for the output of refmap.

```
#!/bin/sh 
module load Stacks 
mkdir output_refmap 
mkdir aligned_samples 
```

I then used a text editor to make a file called popmap.txt. This file excluded 8 samples that were not part of my data but were run on the same lane. I imported popmap.txt into the alignment folder. I was then able to run refmap, specifying the number of CPU's the programme was running as 4 (-t).

```
ref_map.pl --samples ./aligned_samples/ --popmap ./popmap.txt -T 4 -o ./output_refmap 
```

An error at sample bj_23 prevented the programme refmap from finishing. Using the command ```ls -lh``` revealed samples that had been run by refmap and their respective read sizes. I identified samples that had low read sizes, and used Fastqc to inspect the total sequences. If samples had less than 200K reads and less than 200bp total sequence on Fastqc, they were removed from the popmap.txt file. Four samples were removed: bj_23, bj_24, gor_07, jim_19. This allowed the programme refmap to finish running.

```
ls -lh

fastqc bj_23.bam
fastqc bj_24.bam
fastqc gor_07.bam
fastqc jim_19.bam
```

I then ran populations using Stacks to obtain a variant call format (VCF) of my samples for future analyses. I used some parameters for filtering, for example I specified a maximum observed heterozygosity rate of 0.65 (--max-obs-het), and tolerated a maximum of 25% missing data at any given site (-r). I also used a programme that inferred patterns of migration and splitting events in population history (--treemix). 

```
populations -P output_refmap/ -M popmap.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -r 0.75  -O output_refmap 
```

```
77678 variant sites remained 
```

I then wanted to find out more about each individual, in case I needed to filter out any more low-quality individuals. To do this I ran VCFtools.

```
module load VCFtools 
cd output_refmap 
vcftools --vcf populations.snps.vcf --missing-indv
```

I then inspected the ouput folder created by VCFtools called "out.imiss" to look for individuals that were missing a lot of data.

Some low-quality individuals with around 20% missing data which I retained as they still had plenty of data: 

```
coal_03   15931     0.203189 
coal_07   15877     0.2025 
coal_01   16436     0.209629 
elb_09	  17196     0.219323 
coal_18   20491     0.261348 
tim_01	   15886    0.202615 
gor_11	   19168    0.244474 
gor_18	   18475    0.235635  
```

Some low-quality individuals with more than 50% missing data which I removed:

```
tim_18	    42644     0.543894 
wil_20	    54882     0.699981 
elb_15      60707     0.774275 
gor_17      70752     0.902391 
gor_19	     71412    0.910809 
gor_09	     73397    0.936127 
bj_18	     77118    0.983585
tim_19	     77409    0.987297 
elb_12	     77486    0.988279 
```

I then removed the nine low-quality individuals from “popmap.txt” and saved this as “popmap_nooutliersmissing.txt” which I used for further analyses.

```
cat popmap.txt | grep -v tim_18 | grep -v wil_20 | grep -v elb_15 | grep -v gor_17 | grep -v gor_19 | grep -v gor_09 | grep -v bj_18 | grep -v tim_19 | grep -v elb_12 > popmap_nooutliersmissing.txt 
```

I re-ran populations using Stacks, tolerating a maximum of 20% missing data at any given site (-r). 

```
populations -P output_refmap/ -M popmap_nooutliersmissing.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -r 0.80  -O output_refmap 
```

68011 variant sites remained for 171 individuals  

I then downloaded the following files:

populations.plink.ped, 
populations.plink.map, 
populations.structure, 
populations.snps.vcf.gz




