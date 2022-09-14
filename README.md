# SNPcalling_2022.md
## Quality control

Firstly I subset the large files (.gz) recieved from X to 250,000 bp reads. As this is paired-end data I repeated this on read 1 and read 2 independently.

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

The fastqc output for read 1 https://jupyter.nesi.org.nz/user/krogr057/lab/workspaces/auto-E/tree/uoo03677/myanalyses/sourcefiles/GK_GM_S1_R1_sample_fastqc.html and read 2 https://jupyter.nesi.org.nz/user/krogr057/lab/workspaces/auto-E/tree/uoo03677/myanalyses/sourcefiles/GK_GM_S1_R2_sample_fastqc.html shows quite a lot of adapter contamination by an Illumina Universal Adapter. 

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

I started by extracting the .key file containing the barcodes used for each of read 1 and read 2 for the samples in my data. I then loaded this file into a text editor "Sublime Text" for formatting this data, according to section 1.4.2 of the Stacks manual. For specifying combinatorial barcodes with sample names; there is one barcode per column, with two columns for each of the read 1 and read 2 barcodes, and sample names in a separate column included in the barcode file, with each column separated by a tab. 

I then specified the barcode type according to section 1.4.2 of the Stacks manual, as paired end with inline barcodes on the single and paired-ends (--inline_inline)

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


