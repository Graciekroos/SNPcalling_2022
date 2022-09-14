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
I then wanted to check the quality of the data using fastqc. The output () shows quite a lot of adapter contamination by an Illumina Universal Adapter. This adapter sequence is “AGATCGGAAGAG”. I then wanted to remove this adapter contamination. My approach was to retain sequences with a minimum (-m) of 30 base pairs in length. 


```
#!/bin/sh 
module load cutadapt 
cutadapt -j 4  -a AGATCGGAAGAG -A AGATCGGAAGAG -m 30:30 -o trimmed_GK_GM_S1_R1_sample.fq -p trimmed_GK_GM_S1_R2_sample.fq GK_GM_S1_R1_sample.fq GK_GM_S1_R2_sample.fq 
fastqc trimmed_GK_GM_S1_R1_sample.fq trimmed_GK_GM_S1_R2_sample.fq 
```

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


