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

Then I extracted the .key file containing the barcodes used for each of read 1 and read 2 of samples. I loaded this file into a text editor for formatting the data, according to section 1.4.2 of the Stacks manual. I named this file "barcodes.txt" and loaded it into the myanalyses folder.

I then loaded Stacks to use the programme process_radtags for demultiplexing.

```
#!/bin/sh 
module load Stacks/2.58-gimkl-2020a 
```

For demultiplexing, I specified the data as paired end (-P), the barcodes.txt file for barcodes, the restriction enzyme used as PstI (-e)  and for the program to rescue barcodes and cut sites (-r). I also wanted to have higher quality data so I specified (-c) to remove any reads with uncalled bases and (-q) to discard reads with low quality scores.

I also specified the barcode type according to section 1.4.2 of the Stacks manual, as paired end with inline barcodes on the single and paired-ends (--inline_inline)

```
process_radtags –p raw/ -P -b barcodes.txt -o ./samples/ -e PstI –r –c –q --inline_inline 
```

This retains about 90% of reads.

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

BWA and samtools were then used to align each of read 1 and read 2 for every sample to the reference stonefly genome using a loop command. All commands can be found in align.sh.

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

I then wanted to run the Stacks programme refmap to identify any low-quality individuals. To do this I firstly created a folder "aligned_samples" and moved the BAM output files from reference alignment into this folder, to be used as input for refmap. I also created the folder "output_refmap" for the output of refmap.

```
#!/bin/sh 
module load Stacks 
mkdir output_refmap 

mkdir aligned_samples 
```

I then made a text file in a text editor called popmap.txt, not including 8 samples that were not part of my data but were run on the same lane, and imported this into the alignment folder. I was then able to run refmap, specifying the number of CPU's the programme was running as 4 (-T).

```
ref_map.pl --samples ./aligned_samples/ --popmap ./popmap.txt -T 4 -o ./output_refmap 
```

Using the command ls -lh reveals samples that had been run by refmap and their respective read sizes. I identified samples that had low read sizes, and used Fastqc to inspect the total sequences. If samples had less than 200K reads and less than 200bp sequence on Fastqc, they were removed. These sequences were removed (bj_23, bj_24, gor_07, jim_19) before continuing. 

```
ls -lh

fastqc bj_23.bam
fastqc bj_24.bam
fastqc gor_07.bam
fastqc jim_19.bam
```

I then ran populations using stacks to obtain a variant call format (VCF) of my samples for future analyses. I specified a maximum observed hetrozygosiy rate of 0.65 (--max-obs-het), tolerating a maximum of 25% missing data (-r), and using a programme that infers patterns of migration and splitting events in population history (--treemix). 

```
populations -P output_refmap/ -M popmap.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -r 0.75  -O output_refmap 
```

```
77678 variant sites remained 
```

I then wanted to find out more about each individual, in case I needed to filter out any low quality individuals. To do this I ran vcftools.

```
module load VCFtools 
cd output_refmap 
vcftools --vcf populations.snps.vcf --missing-indv
```

I then inspected the ouput folder created by vcftools called "out.imiss" to look for individuals that were missing a lot of data.

Some low-quality individuals with around 20% missing data which I retained as they still had plenty of data: 

```
coal_03   15729    0.20249 
coal_07   15695    0.202052 
coal_01   16204    0.208605 
elb_09	   16926    0.2179 
coal_18   20333    0.26176 
tim_01	   15798    0.203378 
gor_11	   19060    0.245372 
gor_18	   18369    0.236476 
```

Some low-quality individuals with more than 50% missing data which I removed:

```
tim_18	    42185    0.543075 
wil_20	    54362    0.699838 
elb_15     60330    0.776668 
gor_17     70113    0.902611 
gor_19	    70788    0.9113 
gor_09	    72940    0.939005 
bj_18	     76409    0.983663 
tim_19	    76685    0.987216 
elb_12	    76747    0.988015 
```

I then removed the nine low-quality individuals from “popmap.txt” and saved this as “popmap_nooutliersmissing.txt” which I used for further analyses.

```
cat popmap.txt | grep -v tim_18 | grep -v wil_20 | grep -v elb_15 | grep -v gor_17 | grep -v gor_19 | grep -v gor_09 | grep -v bj_18 | grep -v tim_19 | grep -v elb_12 > popmap_nooutliersmissing.txt 
```









