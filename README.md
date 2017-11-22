# ALLHiC
A package to scaffolding polyploidy or heterozygosis genome using Hi-C data 

## Introduction
    

### Algorithm demo  

### Citation
    

### Installation

### Running the pipeline

- **Data Preparation**  
    - Closely related 'diploid' species (chromosomal level) with gene annotation
    - Hi-C reads (with some coverage)
    - Contig level assembly of target genome with decent N50 and gene annotation  

- **Identification of alleles based on BLAST results**  
> Blast CDS in target genome to CDS file in close related reference  
> Note: Please modify cds name before running BLAST. The cds name should be same with gene name present in GFF3   

```
$ blastn -query rice.cds -db Bd.cds -out rice_vs_Sb.blast.out -evalue 0.001 -outfmt 6 -num_threads 4 -num_alignments 1
```
> Remove blast hits with identity < 60% and coverage < 80%  
```
./blastn_parse.pl -i rice_vs_Bd.blast.out -o Erice_vs_Bd.blast.out -q ../4_rice_alleles/T2/riceT2.cds.fasta -b 1 -c 0.6 -d 0.8 
```
> Classify alleles based on BLAST results
```
./classify -i Eblast.out -p 2 -r Bdistachyon_314_v3.1.gene.gff3 -g riceT2.gff3   
```
> After running the scripts above, two tables will be genrated. Allele.gene.table lists the allelic genes in the order of diplod refernece genome and Allele.ctg.table lists corresponding contig names in the same order.   

- **Map Hi-C reads to draft assembly** 
> use bwa index and samtools faidx to index your draft genme assembly  
```
bwa index -a bwtsw draft.asm.fasta  
samtools faidx draft.asm.fasta  
```
> Aligning Hi-C reads to the draft assembly  
```
bwa aln -t 24 draft.asm.fasta reads_R1.fastq.gz > sample_R1.sai  
bwa aln -t 24 draft.asm.fasta reads_R2.fastq.gz > sample_R2.sai  
bwa sampe draft.asm.fasta sample_R1.sai sample_R2.sai reads_R1.fastq.gz reads_R2.fastq.gz > sample.bwa_aln.sam  
```
> Filtering SAM file 
```
PreprocessSAMs.pl sample.bwa_aln.sam draft.asm.fasta MBOI
perl ~/software/script/filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam  
samtools view -bt draft.asm.fasta.fai sample.clean.sam > sample.clean.bam  
```

- **Prune BAM files**  
> Next, we will used Allele.ctg.table to prune 1) signals that link alleles and 2) weak signals. 
```  
./prunning.pl -i Allele.ctg.table -b bam.list -r draft.asm.fasta   
```
- **Partion contigs and Hi-C reads into allele groups**
>And then use Allele.gene.table to partition contigs and prunning bam into allele groups.
```
./partition.pl -g Allele.gene.table -r draft.asm.fasta
```

- **Optimize order and orientation**
- **Chromosome linking**
- **Postprocessing**
