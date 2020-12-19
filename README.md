# Project-1-semester-IB
Authors: 

- Danko Katerina 
- Anton Sidorin
- Sogomonyan Karina
- Ilyutkin Stanislav

### Introduction 
Habitat adaptation is one of the key factors in evolutionary success and species radiation. Voles inhabit different biotopes worldwide. Close relatives with a little differences in genome live in forests and in rocky mountains. Moreover cryptic species live in the same еcological niches, but differ genetically. Vole adaptation to different habitats occurs due to changes in gene expression. However, the differentially expressed genes were still not fully understood for voles from different habitats and cryptic species from the same habitat. 
We performed differential expression analysis for close relative species, inhabiting forests (*C. glareolus*), forest mountains (*C. nivalis*), rocky mountains (*A. lemminus*) and two cryptic species *L. raddei* & *L. gregalis* from steppe.

**Goal**
To carry out a comparative analysis of the differential expression of the vole representatives living in different niches

**Tasks**
1. Estimate the quality of raw reads
2. Align reads on the reference genome of Microtus ochrogaster
3. Find differentially expressed genes
4. Carry out the gene ontology analysis

### Data
For all species, except for *C.glareolus*, tissues of the following organs were sequenced: muscle, heart, lungs, testes, and brain. *C. glareolus* tissues were not sequenced by us, but taken from the NCBI database. Paired reads of the heart, liver and spleen were found for this species.

![](Personal_comparison.jpg)

### Quality control
First of all, we checked the quality of our reads using FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). We used the MultiQC tool (https://multiqc.info/) for a more visual presentation of the quality control data (see attached materials) . According to the preliminary assessment of the reads, it can be argued that the data is of good quality.  However, at the beginning of the reads, an uneven distribution of ~ 15 bases. This may indicate the presence of adapters, as well as simply the poor quality of the reads at the beginning of the sequences. We also got rid of very short reads (less than 30 in length), and those with several nucleotides in a row are of low quality.  Samples were trimmed by Trimmomatic-0.39 (http://bitinfo.colorado.edu/biofrontiers-core-facility-workshops/previous-biofrontiers-workshops/short-read-2015/day-4/trimmomatic-manual) For adapter removing we additionally used AdapterRemoval (v. 2.3.1)

Example of command for trimming:
```
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar  PE -phred33 5441_S37_R1_001.fastq.gz  5441_S37_R2_001.fastq.gz C_nivalis_pairs_1-2.fastq.gz  C_nivalis_single_1-2.fastq.gz C_nivalis_pairs_2-2.fastq.gz  C_nivalis_single_2-2.fastq.gz HEADCROP:10  SLIDINGWINDOW:4:15 MINLEN:30
```
New trimmed files were checked with FastQC again  (see attached materials) . As a result, we got rid of low quality nucleotides at the beginning of the reads, adapters, short reads and low quality reads.

### Alignment
We decided to align reads to the genome of well characterised species - Microtus ochrogaster, which is the closest relative to our experimental species and has the most fully assembled genome and its annotation.
Reference genome - M.ochragaster
- GCF_000317375.1_MicOch1.0_genomic.fna.gz (genome, fasta format)
- GCF_000317375.1_MicOch1.0_genomic.gff.gz (genome annotation, GFF format)

We mapped our reads via HISAT2 programm. In order to start the alignment, indexes of the reference genome must be built:
```
hisat2-build /path/to/reference_genome.fna genome_index
```
For species closer to the reference (C.nivalis), we used the command:
```
hisat2 -p 2 -x genome_index -1 <forward reads> -2  <reverse_reads>| samtools view -b > aligned.bam
```
For more distant (C. glareolus, L. gregalis A ,  L. raddei, A.lemminus ) we used more soft parameters , as the overall alignment with standard parameters rate was low:
```
hisat2 --mp 2,1 --sp 1,1 -p 2  -x genome_index -1 <forward reads> -2  <reverse_reads>| samtools view -b > aligned.bam
```
Results of alignment:

