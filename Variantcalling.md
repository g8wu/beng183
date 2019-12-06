# Variant Calling
BENG 183 Group 1
Authors: Haoyin Xu, Hongru Yu, Ginny Wu

## Introduction
The genomes of individuals and overall populations are all incredibly similar; humans share 99.9% of our DNA while the remaining 0.1% of variations instruct our diversity.[<sup>[1]</sup>](https://www.genome.gov/17516714/2006-release-about-whole-genome-association-studies) These variations arise from random mutations as well as from gene recombinations in the germ line. When we compare genomes, the variations can be searched for and used to study anything from diseases to body development. For example, variations can be used to map the development of cell lineages in an embryo or the growth of a cancerous tumor.

<p float="left">
  <img src='/pictures/snp.png' width='300'/>
  <img src='/pictures/indel.png'  width='500'/>
</p>

*Figure 1a: SNP Example*

*Figure 1b: Indel Example*

  
The specific ways variant DNA sequences appear can be categorized into three groups: SNPs, indels, and structural variations. Single nucleotide polymorphisms (SNPs) represent differences of a single nucleotide. Indels are insertions and deletions of segments of DNA and not as common as SNPs. Structural variations are much larger and typically characterized as more than 1 kb in length. These segments can be inverted, translocated, or copied redundantly within the genome. Variant calling is the process by which these variations are identified from sequence data.

## Strategies

Different variant calling methods rely on several kinds of general strategies, including probabilistic strategy, heuristic strategy, and machine learning. Each of these approaches has its own advantages and disadvantages, and researchers may choose depending on the actual data and sample type.
#### Bayes' Method
The probabilistic approach takes a Bayesian perspective on the data. Researchers use the data to generate prior estimates for genotype probabilities (**P(G)**), create error models for data observations (**P(D|G)**), and combine these steps to calculate the probabilities of variants at certain loci. During these calculations, researchers have to consider the effects of linkage disequilibrium, which makes genotypes at adjacent loci not independent.

<img src='/pictures/bayes.svg'>

*Figure 2: Bayes’ Theorem*
<br>

#### Heuristic Method
Heuristic based algorithms serve as an alternative method. Instead of calculating genotype possibilities, researchers would use a list of heuristic factors to set the bounds for variant calling. Those factors might include minimum allele counts, read quality cut-offs, and depth levels of read coverage. Though a relatively unpopular approach, the method could robustly outly data that violate the assumptions of probabilistic models.

#### Machine Learning Method
Machine learning represents researchers’ recent attempts to optimize the current variant calling methods. Relying on convolutional neural network (CNN), the method is able to magically output genotype likelihoods. Researchers currently have few practical ways to understand the nature of these neural networks, but try their best to make sure that the accuracy of input data meet their expectations.

## Procedure (to be filled)
After the whole genome or exome is sequenced, the raw reads in FASTQ files are quality checked and aligned to the reference genome, resulting in BAM or CRAM files. The alignment shows how the sequences are different from the reference genome, and these variants can be further analyzed to confirm their significance.

Because this is a comparative analysis, the algorithms can differ depending on sample type. There are many packages and pipelines that have been developed to accommodate for diploidy, somatic cells, and germline cells.

## Significance (to be filled)
More discussion on applications of variant calling. Examples include [Feliciano, 2018].(https://doi.org/10.1016/j.tube.2018.04.003)

## Demo
We demonstrate an analysis pipeline for identifying tuberculosis related SNPs starting from analysis ready reads to final VCF file visualizations. The tools used:
*  **FastQC**: quality check and trimming of reads[<sup>[3]</sup>](http://www.ncbi.nlm.nih.gov/pubmed/19451168)
* **bwa mem**: Burrows-Wheeler Alignment alignment[<sup>[3]</sup>](http://www.ncbi.nlm.nih.gov/pubmed/19451168)
* **samtools**: file conversion sam to bam[<sup>[3]</sup>](http://samtools.sourceforge.net)
* **VarScan**: variant calling[<sup>[3]</sup>](http://dkoboldt.github.io/varscan/)

Demo files:----------------

[Full pipeline script](https://github.com/g8wu/beng183/blob/master/run_variance.txt)

### FastQC before and after trimming
<p float="left">
  <img src='/pictures/fastqc1.png' width='400'/>
  <img src='/pictures/fastqc1-trim.png'  width='400'/>
  <img src='/pictures/fastqc2.png' width='400'/>
  <img src='/pictures/fastqc2-trim.png'  width='400'/>
</p>
<br>

(prefixes: ERR2432987 ERR2432988 ERR2432989 ERR2433004 ERR2433005 ERR2433006)
Setup output file and reference genome file:
```
file=log.txt
end=.fastq.gz
ref=tuberculosis.fasta
bwa index $ref
```

Run Fastqc to quality check reads
```
fastqc -o . \${prefix}_1${end} \${prefix}_2${end}
```

Using Sickle, trim ends with QC score threshold 30
```
sickle pe -q 30 -f \${prefix}_1${end} -r \${prefix}_2${end} -t sanger \
-o ${prefix}_t1.fastq -p ${prefix}_t2.fastq -s singletons.fastq \
```

Quality check reads with trimmed ends
```
fastqc -o . ${prefix}_t1.fastq ${prefix}_t2.fastq
```
	  
Using bwa align tuberculosis sequences to the reference genome
```
bwa mem $ref ${prefix}_t1.fastq ${prefix}_t2.fastq > ${prefix}.sam
```

Using samtools, quality check alignment sequences and convert sam to bam file format
```
samtools flagstat ${prefix}.sam
samtools view -S -b ${prefix}.sam > ${prefix}.bam
```

Sort bam files, make pileup files
```
samtools sort ${prefix}.bam > ${prefix}_s.bam
samtools index ${prefix}_s.bam
samtools mpileup -f $ref ${prefix}_s.bam > ${prefix}.mpileup
```

Use VarScan for variant calling
```
java -jar VarScan.jar mpileup2snp ${prefix}.mpileup --min-var-freq 0.90 \
--variants --output-vcf 1 > ${prefix}_raw.vcf
```

Format vcf file and delete temporary files
```
awk '{if (NR>24) $1="Chromosome"; print}' ${prefix}_raw.vcf> ${prefix}.vcf
rm *.sam *.bam *.mpileup *raw.vcf
```

## Alternative Methods
### Galaxy Tool Variant Calling Pipeline
Another demo using the Galaxy tool for variant calling in different settings (diploid/haploid, somatic/germline).


  

## References:
[1] [Fun statistic](https://www.genome.gov/17516714/2006-release-about-whole-genome-association-studies)
[2] [Wikipedia for variant calling](https://en.wikipedia.org/wiki/SNV_calling_from_NGS_data)

[3] [Burrows-Wheeler Alignment](https://github.com/lh3/bwa) (bwa):
* Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID:  [19451168](http://www.ncbi.nlm.nih.gov/pubmed/19451168)]

[4] [VarScan 2](http://dkoboldt.github.io/varscan/):
* Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing  Genome Research  DOI:  [10.1101/gr.129684.111](http://dx.doi.org/10.1101/gr.129684.111)  

[x] [Basic pipeline](https://datacarpentry.org/wrangling-genomics/04-variant_calling/index.html)

[x] [Galaxy pipelines](https://galaxyproject.github.io/training-material/topics/variant-analysis/)





<!--stackedit_data:
eyJoaXN0b3J5IjpbMTUwODIxOTYzMSwxNzA2MTg1NjkyLC0xNj
AzNTY5MDEsLTE0MTY1MTI3NjAsLTcxOTIyOTY2OSwtNjA0NTYx
MTQyLDE3NzMwNTU3MDQsLTIxOTk0MjAyMSwtNTE0Mzk2NjM1LC
04NTA0MTA3MzYsMTYyMDg3MzIyNiwtNDU0MzEzNjQzXX0=
-->
