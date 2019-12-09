# Variant Calling
BENG 183 Group 1

Authors: Haoyin Xu, Hongru Yu, Ginny Wu

## Introduction
The genomes of individuals and overall populations are all incredibly similar; humans share 99.9% of our DNA while the remaining 0.1% of variations instruct our diversity.[<sup>[1]</sup>](https://www.genome.gov/17516714/2006-release-about-whole-genome-association-studies) These variations arise from random mutations as well as from gene recombinations in the germ line. When we compare genomes, the variations can be searched for and used to study anything from diseases to body development. For example, variations can be used to map the development of cell lineages in an embryo or the growth of a cancerous tumor.
<p float="left">
  <img src='/pictures/snp.png' width='200'/>
</p>

*Figure 1a: SNP Example*

<p float="left">
  <img src='/pictures/indel.png'  width='350'/>
</p>

*Figure 1b: Indel Example*


The specific ways DNA variants appear can be categorized into three groups: SNPs, indels, and structural variations. Single nucleotide polymorphisms (SNPs) represent differences of a single nucleotide. Indels are insertions and deletions of DNA segments and not as common as SNPs. Structural variations are much larger and typically characterized as more than 1 kb in length. These segments can be inverted, translocated, or copied redundantly within the genome. Variant calling is the process by which these variations are identified from sequence data.

## Strategies

Different variant calling methods rely on several kinds of general strategies, including probabilistic strategy, heuristic strategy, and machine learning.[<sup>[2]</sup>](https://en.wikipedia.org/wiki/SNV_calling_from_NGS_data) Each of these approaches has its own advantages and disadvantages, and researchers' choice depends on the actual data and sample type.[<sup>[3]</sup>](https://doi.org/10.1186/gm432)

#### Bayes' Method
The probabilistic approach takes a Bayesian perspective on the data. Researchers use the data to generate prior estimates for genotype probabilities (**P(G)**), create error models for data observations (**P(D|G)**), and combine these steps to calculate the probabilities of variants at certain loci. During these calculations, researchers have to consider the effects of linkage disequilibrium, which makes genotypes at adjacent loci not independent.

<img src='/pictures/bayes.svg'>

*Figure 2: Bayes’ Theorem*
<br>

#### Heuristic Method
Heuristic based algorithms serve as an alternative method. Instead of calculating genotype possibilities, researchers would use a list of heuristic factors to set the bounds for variant calling. Those factors might include minimum allele counts, read quality cut-offs, and depth levels of read coverage. Though a relatively unpopular approach, the method could robustly outlay data that violate the assumptions of probabilistic models.

#### Machine Learning Method
Machine learning represents researchers’ recent attempts to optimize the current variant calling methods. Relying on convolutional neural network (CNN), the method is able to magically output genotype likelihoods. Researchers currently have few practical ways to understand the nature of these neural networks, but try their best to make sure that the accuracy of input data meet their expectations.

## Procedure
After the whole genome or exome is sequenced, the raw reads in FASTQ files are quality checked and aligned to the reference genome, resulting in BAM or CRAM files. The alignment shows how the sequences are different from the reference genome, and these variants can be further analyzed to confirm their significance.

Because this is a comparative analysis, the algorithms can differ depending on sample type. There are many packages and pipelines that have been developed to accommodate for diploidy, somatic cells, and germline cells.

## Demo
We demonstrate an analysis pipeline for identifying tuberculosis related SNPs starting from analysis ready reads to final VCF file visualizations. The tools used:
*  **FastQC**: quality check and trimming of reads[<sup>[4]</sup>](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* **bwa mem**: Burrows-Wheeler Alignment alignment[<sup>[5]</sup>](http://www.ncbi.nlm.nih.gov/pubmed/19451168)
* **samtools**: file conversion sam to bam[<sup>[6]</sup>](http://samtools.sourceforge.net)
* **VarScan**: variant calling[<sup>[7]</sup>](http://dkoboldt.github.io/varscan/)

### Demo Files
Sequenced read files: Feliciano, Cinara S. et al. (2018)[<sup>[8]</sup>](https://doi.org/10.1016/j.tube.2018.04.003)

Reference Genome: [Mycobacterium tuberculosis H37Rv NCBI database](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3?report=fasta)

[Full pipeline script on demo files](https://github.com/g8wu/beng183/blob/master/run_variance.txt): implements loop for running multiple read files and sorting them for variant calling. See output file 'log.txt' for runtime status notes.

Run Fastqc to quality check reads. `-o .` outputs files to current directory. Multiple files can be checked using one command line.
```
fastqc -o . \path\to\read_1.fastq.gz \path\to\read_2.fastq.gz
```

### FastQC: Per Base Sequence Quality before and after trimming
<p float="left">
  <img src='/pictures/fastqc1.png' width='400'/>
  <img src='/pictures/fastqc1-trim.png'  width='400' />
</p>

*Figure 3a: Per base sequence quality of read file 1 before and after trimming*

<p float="left">
  <img src='/pictures/fastqc2.png' width='400'/>
  <img src='/pictures/fastqc2-trim.png'  width='400'/>
</p>

*Figure 3b: Per base sequence quality of read file 2 before and after trimming*
<br>
<br>

Using Sickle, trim ends with QC score threshold 30
```
sickle pe -q 30 -f \path\to\read_1.fastq.gz -r \path\to\read_2.fastq.gz -t sanger \
-o \trimmed\read_1.fastq -p \trimmed\read_2.fastq -s singletons.fastq \
```

Using FastQC again, quality check reads with trimmed ends
```
fastqc -o . \trimmed_read_1.fastq \trimmed_read_2.fastq
```

Using bwa, align tuberculosis sequences to the reference genome and output as a sam file.

```
bwa mem tuberculosis.fasta \trimmed_read_1.fastq \trimmed_read_2.fastq > output_align.sam
```

Using samtools, quality check the alignment sequence and convert it from sam to bam file format
```
samtools flagstat output.sam
samtools view -S -b output_align.sam > output_align.bam
```

Sort bam file and make pileup file
```
samtools sort output_align.bam > sorted_align.bam
samtools mpileup -f tuberculosis.fasta sorted_align.bam > output.mpileup
```

Use VarScan for variant calling
```
java -jar VarScan.jar mpileup2snp output.mpileup --min-var-freq 0.90 \
--variants --output-vcf 1 > output_raw.vcf
```

Format vcf file using 'awk' and delete temporary files
```
awk '{if (NR>24) $1="Chromosome"; print}' output_raw.vcf> formatted_output.vcf
rm *.sam *.bam *.mpileup *raw.vcf
```

Variant calling complete! You can view the final vcf file with variant locations and sequnces using any plain text reader or spreadsheet software (i.e. import to Microsoft Excel and select "Tab" when it prompts you what delimiter is used).

## Results
Upon inspecting the vcf files, the following significant mutations can be observed:

Number of mutations at locus:

|     |rpoB |katG |pncA |rpsl |gyrA |embA |embB |embC |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|ERR2432987|1|1|0|0|0|0|0|1|
|ERR2432988|2|1|1|0|0|0|0|1|
|ERR2432989|2|1|1|1|0|0|1|1|
|ERR2433004|0|1|0|1|4|0|0|1|
|ERR2433005|2|1|1|0|3|1|1|1|
|ERR2433006|0|0|0|1|4|0|0|1|

Type of mutations at locus:
(M = Missense, S = Synonymous, N = None)

|     |rpoB |katG |pncA |rpsl |gyrA |embA |embB |embC |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|ERR2432987|S|M|N|N|N|N|N|M|
|ERR2432988|MS|M|M|N|N|N|N|M|
|ERR2432989|MS|M|M|M|N|N|M|M|
|ERR2433004|N|M|N|M|MMMS|N|N|S|
|ERR2433005|MS|M|M|N|MMM|S|S|S|
|ERR2433006|N|N|N|M|MMMS|N|N|S|

## Significance
A patient goes to hospital after showing classical food poisoning symptoms. After rapid testing, an E. Coli infection was confirmed. However, the specific type of E. Coli he contracted is antibiotic resistant. Luckily, because of the growth rate of E. Coli, specific testing can be done with different antibiotic trials to figure out which one will improve his condition at a relatively fast rate. But if the patient is unlucky and contracted a slow-growing pathogen like Mycobacterium tuberculosis, the specific testing that need to be done will take as long as eight weeks to complete. During that process, doctors must take a leap of faith and guess which antibiotic mix he needs to use to prevent the patient's condition from worsening. [<sup>[9]</sup>](https://www.nature.com/news/health-care-bring-microbial-sequencing-to-hospitals-1.15282)

However, new technologies have emerged since the old days. Now we can sequence the entire genome of Mycobacterium Tuberculosis in a matter of days. There have already been multiple mutations known to cause tuberculosis to develop resistance.Bringing pathogen sequencing to hospitals will being several benefits. First, we can improve the survival rate of patients especially those with more rapid and lethal infections. Second, we will significantly decrease the amount of broad-spectrum antibiotics used when specific testing of the pathogen is being done. This will reduce the possibility of the pathogen to develop new resistances. Finally, with enough sequencing done on the same pathogen, there is a bigger likelihood that new mutations that causes antibiotic resistance will be discovered. [<sup>[9]</sup>](https://doi.org/10.1016/j.tube.2018.04.003)

Now focusing on the mutations detected and how they affect protein expression. We are particularly interested in the mutations that may cause the resistance:

1. **rpoB**: The rpoB gene encodes the β subunit of bacterial RNA polymerase. It is the second-largest polypeptide in the bacterial cell. It is the site of mutations that confer resistance to the rifamycin antibacterial agents, such as rifampin. Mutations in rpoB that confer resistance to rifamycin do so by altering residues of the rifamycin binding site on RNA polymerase, thereby reducing rifamycin binding affinity for rifamycin. [<sup>[10]</sup>](https://www.sciencedirect.com/science/article/pii/S0092867401002860?via%3Dihub)
2. **katG**: Aside from providing energy, proline metabolism influences oxidative stress resistance in different organisms. Proline increases oxidative stress tolerance in tuberculosis bacteria via a proadaptive effect involving endogenous hydrogen peroxide production and enhanced catalase-peroxidase activity. katG increases proline presence in tuberculosis bacteria, making them resistance to certain antibiotics. [<sup>[11]</sup>](https://jb.asm.org/content/197/3/431)
3. **pncA**: PncA is a gene encoding pyrazinamidase in Mycobacterium species. Pyrazinamidase converts the drug pyrazinamide to the active form pyrazinoic acid. There is a strong correlation between mutations in pncA and resistance of M. tuberculosis to pyrazinamide. [<sup>[12]</sup>](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2346646/)
4. **rpsl**: Mutations in this region result in streptomycin resistance. The S12 protein is a component of the 30S subunit of the ribosome and plays a role in translational accuracy. [<sup>[13]</sup>](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC163252/)
5. **gyrA**: Mutations in this region is known to inhibit intake of ciprofloxacin by altering the transport protein structure. [<sup>[14]</sup>](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC148966/)

Because of the specific proteins that the mutations affect, there can be concluded a guideline for ideal treatment for the patient contracted the specific *Mycobacterium Tuberculosis* strain that we mapped. We need to use antibiotics that targets the portion of the cell division pathway that takes place inside the cell. So, no matter how the cell membrane surface changes, it will not fail to target the pathogen. We also prefer antibiotics that has high osmolarity, so it can pass through the cell membrane more effectively. Last, we need to increase the dosage of the antibiotics for the mutation in gyrA may cause excessive effluxion of the antibodies thus decreasing their effectiveness. Doctors should not prescribe rifamycin, pyrazinamide, streptomycin or ciprofloxacin for there are mutations specifically responsible for resistance to those drugs. 

## Alternative Methods
### Galaxy Tool Variant Calling Pipeline
Other variations of the process can use different algorithms to identify variations. The bioinformatics webtool Galaxy includes several pipelines for variant analysis that are optimized for different sample types. The tutorials found [here](https://galaxyproject.github.io/training-material/topics/variant-analysis/) take into consideration sequence characteristics of diploid, haploid, somatic, and germline sample reads.[<sup>[15]</sup>]


## References:
[1] [Fun statistic](https://www.genome.gov/17516714/2006-release-about-whole-genome-association-studies)

[2] [Wikipedia for variant calling](https://en.wikipedia.org/wiki/SNV_calling_from_NGS_data)

[3] Pipeline Concordance
* O'Rawe, Jason et al. (2013). Low concordance of multiple variant-calling pipelines: practical implications for exome and genome sequencing. Genome Medicine, 5:28. DOI: [10.1186/gm432](https://doi.org/10.1186/gm432)

[4] [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[5] [Burrows-Wheeler Alignment](https://github.com/lh3/bwa) (bwa):
* Li H. and Durbin R. (2009). Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. PMID: [19451168](http://www.ncbi.nlm.nih.gov/pubmed/19451168)

[6] [Samtools](http://samtools.sourceforge.net)
* Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009). The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25:2078-9. PMID: [19505943](http://www.ncbi.nlm.nih.gov/pubmed/19505943)

* Li H. et al. (2011). Statistical Framework for SNP Calling, Mutation Discovery, Association Mapping and Population Genetical Parameter Estimation from Sequencing Data. Bioinformatics, 27-21:2987-93. PMID: [21903627](http://www.ncbi.nlm.nih.gov/pubmed/21903627)

[7] [VarScan](http://dkoboldt.github.io/varscan/):
* Koboldt, D. et al. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing  Genome Research. DOI: [10.1101/gr.129684.111](http://dx.doi.org/10.1101/gr.129684.111)  

[8] Whole Genome Sequencing Accuracy
* Feliciano, Cinara S. et al. (2018). Accuracy of whole genome sequencing versus phenotypic (MGIT) and commercial molecular tests for detection of drug-resistant Mycobacterium tuberculosis isolated from patients in Brazil and Mozambique. Tuberculosis, 110:59-67. DOI: [10.1016/j.tube.2018.04.003](https://doi.org/10.1016/j.tube.2018.04.003)

[9] Sharon Peacock. "Health care: Bring microbial sequencing to hospitals." Nature. 22 May 2014. Print.

[10] Campbell, E.A.; Korzheva, N.; Mustaev, A.; Murakami, K.; Nair, S.; Goldfarb, A.; Darst, S.A. (2001). "Structural mechanism for rifampicin inhibition of bacterial RNA polymerase". Cell. 104:101–12. doi:10.1016/S0092-8674(01)00286-0. PMID 11290327

[11] Zhang, Lu, James R. Alfano, and Donald F. Becker. “Proline Metabolism Increases katG Expression and Oxidative Stress Resistance in Escherichia Coli.” Ed. P. de Boer. Journal of Bacteriology 197.3 (2015): 431–440. PMC. Web. 11 June 2018.

[12] Juréen, Pontus; Werngren, Jim; Toro, Juan-Carlos; Hoffner, Sven (2016-12-15). "Pyrazinamide Resistance and pncA Gene Mutations in Mycobacterium tuberculosis". Antimicrobial Agents and Chemotherapy. 52 (5): 1852–1854. doi:10.1128/AAC.00110-08. ISSN 0066-4804. PMC 2346646 Freely accessible. PMID 18316515

[13] Sreevatsan, S et al. “Characterization of rpsL and Rrs Mutations in Streptomycin-Resistant Mycobacterium Tuberculosis Isolates from Diverse Geographic Localities.” Antimicrobial Agents and Chemotherapy 40.4 (1996): 1024–1026. Print.

[14] Minnick, Michael F. et al. “gyrA Mutations in Ciprofloxacin-Resistant Bartonella Bacilliformis Strains Obtained In Vitro.” Antimicrobial Agents and Chemotherapy 47.1 (2003): 383–386. PMC. Web. 11 June 2018.

[15] [Galaxy pipelines](https://galaxyproject.github.io/training-material/topics/variant-analysis/)

[x] High-Coverage Samples
* Li, Heng (2014). Toward better understanding of artifacts in variant calling from high-coverage samples. Bioinformatics, 30-20:2843-2451. DOI: [10.1093/bioinformatics/btu356](https://doi.org/10.1093/bioinformatics/btu356)
