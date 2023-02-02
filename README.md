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

I then needed to clean the output of process_radtags.

#!/usr/bin/env python

Clean3pAdapteronShortDemuxReads.py
Ludovic Dutoit, 2022, email dutoit.ludovic@gmail.com for questions
##To clean the output of process_radtags (Stacks, Catchen et al.) This script is a quick solution to remove the second barcode on the first read or the first barcode on the reverse read.
It is not barcode specific. It is simply looking for a few random bases trailing the second restriction cutsite before the end of the sequence.
The files have to be cleaned from adapters before using this script.
In this case, THE ENZYME ON BOTH SIDE IS THE SAME (i.e. GBS protocol, Elshire et al. 2011), it would have to be adapted for different enzyymes on 3' of forward or reverse.
INPUT files, typically output of process radtags. Files have to be gzipped.
OUTPUT files are not gzipped, it is a new folder with all the sequence files in the input folder but cleaned of trailing barcodes
Parameters lines 16-20, to adapt


import gzip,re,os, argparse


##PARAMETERS TO ADAPT




#parser
parser = argparse.ArgumentParser() # add the parser
parser.add_argument("input_folder",help="Folder with all the gzipped output files of process_radtags") # add the parser
parser.add_argument("cleaned_folder", help=" output folder",type=str)
parser.add_argument("-e","-enzyme_cutsite", help="Enzyme cutsite (3' of both F and R reabd (GBS protocol)", type=str,default="TGCA")
parser.add_argument("-m","-max_barcode_length",help="Maximum length of barcode in the dataset",  type=int,default=8) 
args = parser.parse_args()

#input_folder="samples_all" # demultiplexed folder with output of process_radtags
#cleaned_folder="samples_all_cleaned" # clean folder for output
#enzyme_cutsite="TGCA" #pst1 # in my case, same enzyme for F and R reads.
#max_barcode_length=8 # Maximum length of barcode to look for.#



#Initialise a bunch of counts
n_matches=0
n_reads=0
total_length=0 # in bp
total_length_matches=0 #in bp
n_files=0

#Find all gzip files in input folder
files = [args.input_folder+"/"+file for file in os.listdir(args.input_folder) if file.endswith("gz")]

#If output folder does not exist, make it

if not os.path.exists(args.cleaned_folder):
		os.mkdir(args.cleaned_folder)

# For each gzip file
for file in files:
	n_files+=1
	print(file,n_files,"out of",len(files))
	original=gzip.open(file,"rt")
	output=open(args.cleaned_folder+"/"+os.path.basename(file)[:-3],"w") # not gzip output
	# the lines below loop record by record
	while True:
		line1 = original.readline() # @ ID line
		line2 = original.readline() #  sequence line ATGGG
		line3 = original.readline() # +
		line4 = original.readline() # Phred quality line
		if not line1 or not line2 or not line3 or not line4: break
		#print(line1)
		n_reads+=1
		total_length+=len(line2.strip())
		 # look for cutsite + at least one trailing bases or a max of args.max_barcode_length before the end of the read
		match=re.findall(args.enzyme_cutsite+"[ATGC]{1,"+str(args.max_barcode_length)+"}$",line2)
		if match: # There is a match
			n_matches+=1
			length_to_remove = len(match[0]) # figure out how much to remove
			total_length_matches+=length_to_remove 
			#remove it but not the enzyme cutsie
			#line2=args.enzyme_cutsite="TGCA"
			line2=line2[0:len(line2)-(length_to_remove)+len(args.enzyme_cutsite)-1]+"\n" # remove the bases from the seq
			line4=line2[0:len(line4)-(length_to_remove)+len(args.enzyme_cutsite)-1]+"\n" # remove the bases from the qual line
		output.write(line1+line2+line3+line4) # write the modified or unmodified line to output
	output.close()

#summary
print("removed patterns",n_matches, "times out of", n_reads," for a total of", (total_length_matches/total_length)*100,"percent of bases")





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




