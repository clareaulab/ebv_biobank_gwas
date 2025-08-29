library(rtracklayer) 
library(dplyr)
library(stringr) 

## import gtf 
# This file can be obtained here: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/
gtf <- import("/data1/lareauc/references/hg38_variant_calling/hg38.ncbiRefSeq.gtf")
# subset to chr6 
gtf2 <- gtf[seqnames(gtf) == "chr6"]
# Remove HLA genes 
catch <- gtf2[str_detect(gtf2$gene_name,"HLA")]
# subset to transcripts 
catch2 <- catch[catch$type == "transcript"]

## Include TAP1/2 and MICA/B 
catch <- gtf2[str_detect(gtf2$gene_name,"^TAP[1|2]$")]
toadd <- catch[catch$type == "transcript"]

catch <- gtf2[str_detect(gtf2$gene_name,"^MIC[A|B]$")]
toadd2 <- catch[catch$type == "transcript"]

joint <- c(catch2,toadd,toadd2) 

# simplify regions to only 100 kB windows for all HLA genes
try1 <- reduce(joint, min.gapwidth=100000)

## How many bases assayed?
width(try1) %>% sum() #1549980 bp total

#chr6:28510120-33480577
#4970457bps used, 4970457 / 1321933 === 3.2 decrease in size so that saves a lot 

out <- try1 %>% as.data.frame() %>% mutate(out = sprintf("%s:%s-%s",seqnames,start,end)) %>% pull(out)

writeLines(out,"small_hla_regions_nounmapped.interval")
