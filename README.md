# SNPcalling_2022.md
## Quality control

Firstly I subsetted the large files (.gz) to 250,000 bp reads as this will be quicker to test the initial structure and adapter content before working on the whole data. As this is paired-end data I repeated this on read 1 and read 2 independently.

```
#!/bin/sh 
# in sourcefiles 
zcat  GK_GM_S1_R1_001.fastq.gz | head -n 1000000 > GK_GM_S1_R1_sample.fq 

#!/bin/sh 
# in sourcefiles 
zcat  GK_GM_S1_R2_001.fastq.gz | head -n 1000000 > GK_GM_S1_R2_sample.fq 
```

```
#!/bin/sh 
module load FastQC 
fastqc *fq 
```
I then wanted to check the quality of the data using fastqc. 

The fastqc output for read 1
https://jupyter.nesi.org.nz/user/krogr057/lab/workspaces/auto-E/tree/uoo03677/myanalyses/sourcefiles/GK_GM_S1_R1_sample_fastqc.html 
and read 2 https://jupyter.nesi.org.nz/user/krogr057/lab/workspaces/auto-E/tree/uoo03677/myanalyses/sourcefiles/GK_GM_S1_R2_sample_fastqc.html 

shows quite a lot of adapter contamination by an Illumina Universal Adapter. 

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

Output by fastqc for read 1 (https://jupyter.nesi.org.nz/user/krogr057/lab/tree/uoo03677/myanalyses/sourcefiles/trimmed_GK_GM_S1_R1__fastqc.html) and read 2 (https://jupyter.nesi.org.nz/user/krogr057/lab/tree/uoo03677/myanalyses/sourcefiles/trimmed_GK_GM_S1_R2__fastqc.html) show adapter contamination was successfully removed.

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

I started by making two new folders; "raw" for inputting my files for demultiplexing, and "samples" for the output of demultiplexing. I then created links to put my trimmed whole data for read 1 and read 2 into "raw".

```
#!/bin/sh 
mkdir raw samples 
cd raw 
ln -s ../trimmed_GK_GM_S1_R1_.fq trimmed_GK_GM_S1_R1_001.fastq 
ln -s ../trimmed_GK_GM_S1_R2_.fq trimmed_GK_GM_S1_R2_001.fastq 
cd .. sourcefiles 
```

```
#!/bin/sh 
module load Stacks/2.58-gimkl-2020a 
 
```

Then I extracted the .key file containing the barcodes used for each of read 1 and read 2 of my samples. I loaded this file into a text editor "Sublime Text" for formatting the data, according to section 1.4.2 of the Stacks manual. For specifying combinatorial barcodes with sample names, there is one barcode per column, with two columns for each of the read 1 and read 2 barcodes, and sample names in a separate column with each column separated by a tab.

I then named this file "barcodes.txt" and loaded it onto Stacks.

I used process_radtags for demultiplexing, specifying the data is paired end (-P), the barcodes.txt file for barcodes, the restriction enzyme used as PstI (-e)  and rescue barcodes and cut sites (-r). I also wanted to have higher quality data so I specified (-c) to remove any reads with uncalled bases and (-q) to discard reads with low quality scores.

Next, I specified the barcode type according to section 1.4.2 of the Stacks manual, as paired end with inline barcodes on the single and paired-ends (--inline_inline)

```
process_radtags –p raw/ -P -b barcodes.txt -o ./samples/ -e PstI –r –c –q –inline_inline 
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

Sample files are now clean and able to be aligned to the reference genome.

The first step is to index the genome. To do this I need to use Burrows-Wheeler Aligner (BWA).

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

```
module load SAMtools 
```

BWA was then used to align each of read 1 and ead 2 for every sample to the reference stonefly genome using a loop. One example command below; all commands can be found in align.sh.

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










