### Exon Level Quantification from RSEM BAM output 

We are aiming to have a consistent method of comparison between gene, transcript and exon concatenation per TF for stage 2 training. As a result, we chose the RSEM (as recommended by ENCODE) to generate the expression levels per gene and transcript. However, the RSEM output does not include explicit exon level quantification for the RNA experiment. RSEM outputs a genome based bam file and a transcript based bam file. -- The idea is to use transcript based bam file to compute the weighted counts per exon. 

We use the transcript BAM file output by RSEM (rsem.transcript.bam) to generate exon level quantifications. 

### Install

- Bioinformatics:
	- SAMtools (http://www.htslib.org/download/)
	- bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html)
  - BEDOPS (https://bedops.readthedocs.io/en/latest/)

### Prepare exon interval file "exonFile"
Using the gtf file and the transcript rsem bam file generated by RSEM for each RNA experiment, we will generate the weighted exon counts -- "expected counts per exon"

### Prepare the RNAseq file from the rsem.transcript.bam
1. Look at reads on forward strand
      
      samtools view -h -F 20 rsem.transcript.bam | samtools view -S -b > rsem.transcript.forward.bam
 
 2. Convert bam to bed file we used bedops in the script bam2bed.sh
      
      bam2bed  < rsem.transcript.forward.bam > rsem.transcript.forward.bed

 3. bedtools intersect -a $exonFile -b $bedFile -wb -wa -sorted
 4. Get weighted count quantifyExons.py
