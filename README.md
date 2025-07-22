# Population-scale sequencing resolves correlates and determinants of latent Epstein-Barr Virus infection

Repository of code and processed data supporting analyses in Nyeo et al. 2025.


## Consortia

For this work, we made extensive use of the [UK Biobank (UKB)](https://ukbiobank.dnanexus.com/login)
and [All of Us (AoU)](https://www.researchallofus.org/data-tools/workbench/). These data sources are available
for qualified researchers to access genetic and phenotypic data on ~750,000 total individuals. 

## Extracting reads per EBV contig

The most cost-efficient way to access the EBV data from the .cram files containing aligned
WGS data is to stream the reads via GATK. This allows you to only pay egress costs 
associated with the aligned sequences rather than the costs of moving the entire cram. 

An example of this command is shown below that was used in a `for` loop for processing 
the All of Us consortium:

```
./gatk-4.2.6.0/gatk PrintReads -I gs://fc-aou-datasets-controlled/pooled/wgs/cram/v7_delta/wgs_*****.cram \
	-L chrEBV -R Homo_sapiens_assembly38.fasta -O ID*****ebv.bam \
	--gcs-project-for-requester-pays $GOOGLE_PROJECT --cloud-prefetch-buffer 0 \
	--cloud-index-prefetch-buffer 0
```

A few notes:
- GATK 4.2.6 turned out to be a critical required version. Newer versions would give strange errors. 
- The reference genome is worth copying locally, which has to be indexed according to GATK's liking.
	- Said reference genome can [be found here](https://github.com/broadinstitute/gatk/blob/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz).
- The asterisks in place should be replaced with the individual ID. 

This takes ~1-2 hours on a single core machine and will occupy only a few gigabytes of disk space. 
Hence, the full analysis for extracting these data should only cost ~10 USD, which fits easily
within the free credits available to AoU researchers. 

For the UK Biobank data, .bam files were accessed on a local cluster. Hence, the 
`samtools view` command was sufficient. To the best of our knowledge, 
`samtools view` is not compatible with `.cram` files hosted on buckets; hence, the gatk command seems necessary.

## EBV reference

Many files in this repository are pulled from the nuccore of the main EBV contig contained in the 
hg38 reference genome (sequenced from the Raji cell line). The full data 
[is available here - nuccore NC_007605.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_007605.1).

## Reproducing downstream analyses

As many analyses require access to individual-level data, many of these are omitted from 
this repository in compliance with data sharing agreements with these consortia. 
What are included allow for more processed analyses to be replicated. 
[Contact Caleb Lareau](lareauc@mskcc.org) for assistance in engineering these features 
on your own workspace. 

All custom analyses were performed using the `R 4.4.0` software environment. 
Set your working directory in `code` in the base of any of these three
to run the R code line-by-line.

## all-of-us-notebooks
This folder contains all of the major notebooks used for assembling the EBV quantification,
performing the EBV GWAS validation analyses, and the phenotypic associations using EBV as
an exposure for complex trait associations. (Fig. 1-3)

## celltype-pathway-mapping
This folder contains code for pathway and cell type enrichment analyses of the 
147 genes with protein-altering variation from the ExWAS analyses (Figure 4). 

## epitope-scoring
This folder contains code needed to create the processed peptide files 
as well as results from running NetMHC (I+II). Additional files for enrichment
analyses of known IEDB epitopes are also contained here (Fig. 5).

## viral-sequences
This folder contains aggregated data of the EBV contig from both consortia, 
including workflows for integrating these allele frequencies with annotations 
of the EBV contig (see ED Fig. 6).



<br><br>
