######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

library(TitanCNA)
library(doMC)
registerDoMC()
options(cores=12)

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
## Users should run TITAN once for each setting of the number of clonal clusters, ranging from 1 to 5.
######################################################################################################
ID_SAMPLE <- "TEST"
CHROM <- "ALL"
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
## chrom <- "chr19"

chrom <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
           "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
           "chr20", "chr21", "chr22", "chrX")
 
# chrom <- c(1:22, "X")

# in the script that is available on GITHUB :
# make_option(c("--chrs"), type = "character", default = "c(1:22, 'X')",
#            help = "Chromosomes to analyze; string [Default: %default"),
# make_option(c("--gender"), type = "character", default = "male", help = "User specified gender: male or female
# [Default: %default]"),

######################################################################################################
######################################################################################################
################################################# reading CENTROMERES :
######################################################################################################
######################################################################################################

centromeres <- read.table("hg38.IDEOGRAM.track.only.CENTROMERES.in.BED.for.TITAN.30sep2018.bed",
                          sep="\t", header=F, stringsAsFactors=F)  

centromeres.df <- as.data.frame(centromeres)

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
############################################## perhaps reading also the IDEOGRAM :
##################################################################################
##################################################################################
### shall we run it on all CHR ie on WHOLE-GENOME ?
###  testing it again in order to see how it works :
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

#!/usr/bin/env Rscript

# library("Rsamtools")
# library("TitanCNA")
# library(doMC)
# registerDoMC()
# options(cores=4)

# tumbamFile <- "SPCG-OS040_8T.GRCh38p5M.MD_intervals_hg38_male.interval_list.IR.RC.bam"
# tumIndexFile <- paste(tumbamFile,".bai",sep = "")
# vcfFile <- "SPCG-OS040_7G.GRCh38p5M.MD_samtools.germline.HETERO.vcf"
# outFile <- "SPCG-OS040_7G.GRCh38p5M.MD_samtools.germline.HETERO.output.Rsamtools.5april2017.out"

# pp <- PileupParam(min_base_quality = 10, min_mapq = 20,
#          min_nucleotide_depth = 10, max_depth = 1000,
#          distinguish_strands = FALSE,
#          distinguish_nucleotides = TRUE)
#
# countsDF <- extractAlleleReadCounts(tumbamFile, tumIndexFile,
#                                     vcfFile, outFile,
#                                     pileupParam = pp)

# data <- loadAlleleCounts(countsDF, symmetric = TRUE, genomeStyle = "UCSC")

######################################################################################################
######################################################################################################
##################################################################################
################################################################################## for a CHR : chr19 :
######################################################################################################
######################################################################################################

## genome.fa.gc.1kb.wig.extracting.chr19.wig
## genome.fa.map.150bp.bin.1kb.wig.extracting.chr19.wig

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

#!/usr/bin/env Rscript

#library(Rsamtools)
#library(TitanCNA)
#library(doMC)
#registerDoMC()
#options(cores=12)

######################################################################################################
######################################################################################################
############### to get the ALLELE COUNTS in the TUMOR :
######################################################################################################
######################################################################################################

# tumbamFile <- "SPCG-OS040_8T.GRCh38p5M.MD_intervals_chr19.IR.RC.bam"
# tumIndexFile <- paste(tumbamFile,".bai",sep = "")
# vcfFile <- "SPCG-OS040_7G.GRCh38p5M.MD_samtools.germline.HETERO.chr19.vcf"
# outFile <- "SPCG-OS040_7G.GRCh38p5M.MD_samtools.germline.HETERO.chr19.Rsamtools.new.june2017.out"

# pp <- PileupParam(min_base_quality = 10, min_mapq = 20,
#          min_nucleotide_depth = 10, max_depth = 1000,
#          distinguish_strands = FALSE,
#          distinguish_nucleotides = TRUE)

# countsDF <- extractAlleleReadCounts(tumbamFile, tumIndexFile,
#                                     vcfFile, outFile,
#                                     pileupParam = pp)


######################################################################################################
######################################################################################################
### the function extractAlleleReadCounts is deprecated; probably not need to use it anymore.
### in fact it does not work anymore ...

# data <- loadAlleleCounts(countsDF, genomeStyle = "UCSC")

### so if we upload the file with the heterozygous SNV and the counts that was computed
### by some other programs, we do :

# the file shall have : chromosome, position,
#                       reference base, reference read counts,
#                       non-reference base, non-reference read counts
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

### SHALL WE START everything from here, as the previous functions are not available ...
### now revising on the date of 13 and 14feb2018:

######################################################################################################
######################################################################################################

### library(Rsamtools)
### library(TitanCNA)
### library(doMC)
### registerDoMC()
### options(cores=12)

#### uploading the ALLELE counts : actually the SCRIPT can start HERE : 
#### LOAD DATA :
#### just to make sure that we have the HETEROZYG SNV and no INDELS..

# x <- read.delim("SPCG-OS040_7G.GRCh38p5M.MD_intervals_ALL.HAPLOTYPE_CALLER.ALL_CHRS.vcf.chr19.HETERO.SNP.BIALLELIC.30sep2018.vcf.v2.rtable", 
#                 header=T, sep="\t", stringsAsFactors=F)

#### including the version for ENTIRE GENOME : obtained with GATK :

x <- read.delim("SPCG-OS040_7G.GRCh38p5M.MD_intervals_ALL.HAPLOTYPE_CALLER.ALL_CHRS.HETERO.SNP.BIALLELIC.06sept2018.rtable", 
                 header=T, sep="\t", stringsAsFactors=F)

######################################################################################################
######################################################################################################
#### possibly we may NEED an additional check here in order to make sure that the ENTRIES are UNIQUE :)
#### reordering X :
######################################################################################################
######################################################################################################

xr <- data.frame(chr=as.character(x$contig),    
                 position=as.numeric(x$position), 
                 ref=as.character(x$refAllele), 
                 refCount=as.numeric(x$refCount),
                 Nref=as.character(x$altAllele),
                 NrefCount=as.numeric(x$altCount), stringsAsFactors=F)

#### the structure of the data :
######################################################################################################
######################################################################################################

head(xr)
str(xr)
dim(xr)  # [1] 57156 ; 2566949

######################################################################################################
######################################################################################################
######## function loadAlleleCounts :
######################################################################################################
######################################################################################################

data <- loadAlleleCounts(xr, genomeStyle = "UCSC", sep = "\t", header = TRUE)

str(data)
head(data) 

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

### write.table(xr, file="SPCG-OS040_7G.GRCh38p5M.MD_samtools.germline.HETERO.chr19.vcf.HET4.SNP.BIALLELIC.vcf.12june2017.rtable.simple.to.usee", sep="\t")

### y <- read.delim("SPCG-OS040_7G.GRCh38p5M.MD_samtools.germline.HETERO.chr19.vcf.HET4.SNP.BIALLELIC.vcf.12june2017.rtable.simple.to.use", 
###                  header=T, sep="\t", stringsAsFactors=F)

#### PREVIOUSLY, we had loaded this data :
#### chr    position    ref    refCount    Nref    NrefCount

#### data <- loadAlleleCounts("SPCG-OS040_7G.GRCh38p5M.MD_samtools.germline.HETERO.output.Rsamtools.v3.chr19.cluster.out",
####                          genomeStyle = "UCSC", sep = "\t", header = TRUE)


### if we want to FILTER the data : 

# data <- filterData(data, chrs="chr19",
#                    minDepth = 10, maxDepth = 1000,
#                    positionList = NULL, centromere = NULL,
#                    centromere.flankLength = 10000)

## or using another file that is produced by GATK :

######################################################################################################
######################################################################################################

# tumWig  <- ""
# normWig <- ""
# gc <-
# map <-

######################################################################################################
######################################################################################################
#########################################################################################
#########################################################################################
######################################################################################### the read depth files :

# tumWig <- "SPCG-OS040_8T.GRCh38p5M.MD_intervals_chr19.IR.RC.1kb.wig.extracting.chr19.wig"
tumWig <- "SPCG-OS040_8T.GRCh38p5M.MD.bam.1kb.wig"

# normWig <- "SPCG-OS040_7G.GRCh38p5M.MD_intervals_chr19.IR.RC.1kb.wig.extracting.chr19.wig"
normWig <- "SPCG-OS040_7G.GRCh38p5M.MD.bam.1kb.wig"

# gcWig <- "genome.fa.gc.1kb.wig.extracting.chr19.wig"               ### for the entire genome : genome.fa.gc.1kb.wig
gcWig <- "genome.fa.gc.1kb.wig"

# mapWig <- "genome.fa.map.150bp.bin.1kb.wig.extracting.chr19.wig"   ### for the entire genome : genome.fa.map.150bp.bin.1kb.wig
mapWig <- "genome.fa.map.150bp.bin.1kb.wig"

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################## GC AND MAPPABILITY CORRECTION
######################################################################################################
######################################################################################################

# genome-wide GC content and mappability score for 1kb windows3

cnData <- correctReadDepth(tumWig, normWig, gcWig, mapWig, genomeStyle = "UCSC")

# find the log RATIO at each heterozy germline SNP

logR <- getPositionOverlap(data$chr, data$posn, cnData)

data$logR <- log(2^logR)                              ####### #transform the log ratio to natural logs

# rm(logR,cnData)

######################################################################################################
######################################################################################################

# data <- filterData(data, "chr19",
#                   minDepth = 10, maxDepth = 1000,
#                   positionList = NULL, centromere = NULL, centromere.flankLength = 10000)

######################################################################################################
######################################################################################################
######################################################################################### FILTER DATA :
######################################################################################################
###  data <- filterData(data,c(1:22,"X"), minDepth=10, maxDepth=200)
###  when running for 1, 2, 3, 4, or 5 clonal cluster for real data
######################################################################################################
######################################################################################################
### it may need filtering for the algorithm to run :)

data <- filterData(data,
                   chrs=chrom,
                   minDepth = 10,
                   maxDepth = 1000,
                   positionList = NULL,
                   centromere = as.data.frame(centromeres),
                   centromere.flankLength = 10000)

######################################################################################################
######################################################################################################
######################################################################################## CLONALITY :
######################################################################################################
# numClusters  :
# '1' treats the sample as being clonal (no subclonality).
# '2' or higher treats the tumour data as being subclonal.
######################################################################################################
######################################################################################################

# If provided and 'symmetric=TRUE', then it will compute the
# median allelic ratio to use as the baseline for heterozygous
# genotypes; otherwise, the baseline will default to 0.55
# (reference/depth) if 'data=NULL'.

######################################################################################################
######################################################################################################
# to add in the name of the files : .numClusters.copyNumbers.
######################################################################################################
######################################################################################################
######################################################################################################

numClusters <- 2  ## the number of clusters will be specified
                  ## when running for 1, 2, 3, 4, or 5 clonal cluster for real data.

copyNumbers <- 5  ## Default (and recommended) is 5. It can go up to 8.

######################################################################################################
######################################################################################################
################################################################################# or using the R code :
######################################################################################################

number_clusters <- 1:2

for (j in 1:2)
{
    numClusters <- number_clusters[j]
   
    params <- loadDefaultParameters(copyNumber = copyNumbers,
                                    numberClonalClusters = numClusters,
                                    symmetric = TRUE,
                                    data = data)

    convergeParams <- runEMclonalCN(data,
                                     gParams=params$genotypeParams,
                                     nParams=params$normalParams,
                                     pParams=params$ploidyParams,
                                     sParams=params$cellPrevParams,
                                     maxiter=20, maxiterUpdate=1500,
                                     txnExpLen=1e15,txnZstrength=1e5,
                                     useOutlierState=FALSE,
                                     normalEstimateMethod="map",
                                     estimateS=TRUE, estimatePloidy=TRUE)

    optimalPath <- viterbiClonalCN(data, convergeParams)

######################################################################################################
######################################################################################################
############################################# here using the CHROMOSOME number for the OUTPUT FIGURES :
######################################################################################################
######################################################################################################
#####################################################################################  chrom <- "chr19"
#######################################################################################################
########################################################################## OUTPUT THE FORMATED RESULTS :

###  Position-specific results that includes genotype, clonal cluster, and cellular prevalence estimates.

out_results <- paste("30sep.", CHROM, ".",numClusters, ".numClusters.", copyNumbers,
                     ".copyNumbers.TITAN.RESULTS.txt", sep="")

results <- outputTitanResults(data,  convergeParams,
                                     optimalPath,
                                     filename = out_results,
                                     posteriorProbs = FALSE,
                                     subcloneProfiles = TRUE)

######################################################################################################
######################################################################################################
str(results)
######################################################################################################
######################################################################################################

### Model parameters and summary values such as the number of clonal clusters and 
### their corresponding cellular prevalence, normal contamination, ploidy, and 
### the S_Dbw validity index for model selection.

######################################################################################################
######################################################################################################
##################### if we have used ALL the CHROMOSOMES we may say chrom=ALL
######################################################################################################
######################################################################################################

out_parameters <- paste("30sep.", CHROM, ".", numClusters, ".numClusters.", copyNumbers,
                                              ".copyNumbers.TITAN.PARAMETERS.txt", sep = "")
                                              
######################################################################################################
######################################################################################################
outputModelParameters(convergeParams,
                         results,  
                         file=out_parameters)
######################################################################################################
######################################################################################################

out_segments <- paste("30sep.", CHROM, ".", numClusters, ".numClusters.", copyNumbers,
                                            ".copyNumbers.TITAN.SEGMENTS.txt", sep = "")

out_segments_IGV <- paste("30sep.", CHROM, ".", numClusters, ".numClusters.", copyNumbers,
                                             ".copyNumbers.TITAN.SEGMENTS.for.IGV.txt", sep = "")

segs <- outputTitanSegments(results,
                            id = ID_SAMPLE,
                            convergeParams,
                            filename = out_segments,
                            igvfilename = out_segments_IGV)

######################################################################################################
######################################################################################################
####################################################################################
#################################################################################### display per CHR
#### ############## ################################################################
######################################################################################################
######################################################################################################

#### ploidy <- tail(convergeParams$phi, 1)
#### normal <- tail(convergeParams$n, 1)

for (chrom in unique(results$Chr))
   {

######################################################################################################
######################################################################################################
### 1. COPY NUMBER ALTERATIONS (LOG RATIO) :
######################################################################################################
######################################################################################################

    ploidy <- tail(convergeParams$phi, 1)
### ploidy
    normal <- tail(convergeParams$n, 1)
### normal

    png(filename=paste("30sep.", chrom, ".", numClusters, ".numClusters.", copyNumbers,
                       ".copyNumbers.TITAN.COPY_NUMBER_alterations.png", sep=""),
                        width=1800, height=600, res=200)

    plotCNlogRByChr(results,
                    chr = chrom,
                    ploidy = ploidy,
                    normal = normal,
                    ylim = c(-2, 2), cex = 0.25,
                    xlab = "",
                    main = paste("CNA on ", chrom, sep=""))
     dev.off()

# The Y-axis is based on log ratios. Log ratios are computed ratios between normalized tumour and normal
# read depths. Data points close to 0 represent diploid, above 0 are copy gains, below 0 are deletions.
# Bright Green - HOMD Green - DLOH Blue - HET, NLOH Dark Red - GAIN Red - ASCNA, UBCNA, BCNA

######################################################################################################
######################################################################################################
#### 2. LOSS OF HETEROZYGOSITY (ALLELIC RATIO) :
######################################################################################################
######################################################################################################

    png(filename=paste("30sep.", chrom, ".", numClusters, ".numClusters.", copyNumbers,
                       ".copyNumbers.TITAN.ALLELIC_RATIOS.png", sep=""),
                       width=1800, height=600, res=200)

    plotAllelicRatio(results,
                          chr = chrom,
                          ylim = c(0, 1), cex = 0.25,
                          xlab = "",
                          main = paste("ALLELIC RATIOS on ", chrom, sep=""))
    dev.off()

# Allelic ratios are computed as RefCount/Depth. Data points close to 1 represent homozygous reference base,
# close to 0 represent homozygous non-reference base, and close to 0.5 represent heterozygous.
# Normal contamination influences the divergence away from 0.5 for LOH events.
# Grey - HET, BCNA Bright Green - HOMD Green - DLOH, ALOH Blue - NLOH Dark Red - GAIN Red - ASCNA, UBCNA
######################################################################################################
######################################################################################################
##### 3. Cellular prevalence and clonal clusters :
######################################################################################################
######################################################################################################

norm <- tail(convergeParams$n, 1)
## norm # estimated normal contamination
### ... [1] 0.1969555
##1 - convergeParams$s[, ncol(convergeParams$s)] # estimated cellular prevalence
### ... [1] 0.9927419 0.5791513
######################################################################################################
######################################################################################################

    png(filename=paste("30sep.", chrom, ".", numClusters, ".numClusters.", copyNumbers,
                 ".copyNumbers.TITAN.CLONAL_CLUSTERS.png", sep=""),
                  width=1800, height=600, res=200)

    plotClonalFrequency(results,
                    chr = chrom,
                    normal = norm,
                    ylim = c(0, 1), cex = 0.25, xlab = "",
                    main = paste("CLONAL CLUSTERS on ", chrom, sep=""))
    dev.off()

## The black horizontal line represents the tumour content labeled as "T".
## Each horizontal grey line
## represents the cellular prevalence of the clonal clusters labeled as Z1, Z2, etc.
## Colours are the sames for allelic ratio plots.

######################################################################################################
######################################################################################################
###### 4. Subclone profiles :
######################################################################################################
######################################################################################################

## Plotting the copy number profiles for the predicted subclones
## (containing 1 or 2 clonal clusters).
## Colours have the same definition as for the allelic ratio plots

    png(filename=paste("30sep.", chrom, ".", numClusters, ".numClusters.",
                 copyNumbers,".copyNumbers.TITAN.SUBCLONE_PROFILES.png", sep=""),
                 width=1800, height=600, res=200)
    plotSubcloneProfiles(results,
                         chr = chrom,
                         cex = 1, spacing = 2,
                         main = paste("SUBCLONE PROFILES on ", chrom, sep="")
                        )
    dev.off()

######################################################################################################
######################################################################################################
###############################################################
#### 5. and if we can combine ALL the figures together :
###############################################################
######################################################################################################
######################################################################################################

    png(filename=paste("30sep.",  chrom, ".", numClusters, ".numClusters.",
                       copyNumbers,".copyNumbers.TITAN.ALL_these_4figures.png", sep=""),
                       width=2000, height=1800, res=200)

    par(mfrow=c(4, 1))   
    plotCNlogRByChr(results,
                         chr = chrom,
                         ploidy = ploidy,
                         normal = normal,
                         geneAnnot = NULL, segs = NULL,
                         ylim = c(-2, 2), cex = 0.25,
                         xlab = "", main = paste("CNA logR on ", chrom, sep=""))

    plotAllelicRatio(results,
                          chr = chrom,
                          geneAnnot = NULL,
                          ylim = c(0, 1), cex = 0.25,
                          xlab = "", main = paste("ALLELIC RATIO on ", chrom, sep=""))

    plotClonalFrequency(results,
                    chr = chrom,
                    normal = norm,
                    geneAnnot = NULL,
                    ylim = c(0, 1), cex = 0.25,
                    xlab = "", main = paste("CLONAL CLUSTERS on ", chrom, sep=""))

    plotSubcloneProfiles(results,
                              chr = chrom,
                              cex = 1, spacing = 2,
                              main = paste("SUBCLONE PROFILES on ", chrom, sep=""))


    dev.off()

   }

## closing the LOOP if we display per CHROMOSOME ..

}

######################################################################################################
######################################################################################################
## closing the LOOP regarding the number of clusters ...
## possibly shall we run only 2 clusters ?
######################################################################################################
######################################################################################################
######################################################################################################
###################################################################################################### HETEROZYG :
######################################################################################################
###################################################################################################### from other file

### re-checking the number of HETEROZYG variants :
### grep -e "#" -e "0/1" VCF

##################################################################
##################################################################
##################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

# SOME comments about the HETEROZYGOUS VARIANTs:

######################################################################################################
######################################################################################################
################################# using GATK : possibly continuing the pipeline with HAPLOTYPE_CALLER : 

VCF="SPCG-OS040_7G.GRCh38p5M.MD_intervals_ALL.HAPLOTYPE_CALLER.ALL_CHRS.vcf"
BAM_TUMOR_chr19="SPCG-OS040_8T.GRCh38p5M.MD_intervals_chr19.IR.RC.bam"

######################################################################################################
######################################################################################################
######################################################
######################################################
################################ using GATK : possibly continuing the pipeline with HAPLOTYPE_CALLER : 
######################################################################################################
######################################################################################################

VCF="SPCG-OS040_7G.GRCh38p5M.MD_intervals_ALL.HAPLOTYPE_CALLER.ALL_CHRS.vcf"
BAM_TUMOR_chr19="SPCG-OS040_8T.GRCh38p5M.MD_intervals_chr19.IR.RC.bam"

### if we would like to include a specific chromosome : -L chr19 \
### as it has been a limitation here to a specific chromosome ...

$GATK \
-T SelectVariants \
-R $REFERENCE_HG38_MALE \
-L chr19 \
-select 'vc.getGenotype("SPCG-OS040_7G").isHet()' \
-selectType SNP \
-restrictAllelesTo BIALLELIC \
-V $VCF \
-o "${VCF}.chr19.HETERO.SNP.BIALLELIC.18aug2018.vcf"

$GATK -T ASEReadCounter \
-R $REFERENCE_HG38_MALE \
-U ALLOW_N_CIGAR_READS \
-I $BAM_TUMOR_chr19 \
-minDepth 2 \
--minBaseQuality 10 \
--minMappingQuality 10 \
-sites "${VCF}.chr19.HETERO.SNP.BIALLELIC.18aug2018.vcf" \
-o "${VCF}.chr19.HETERO.SNP.BIALLELIC.18aug2018.vcf.v2.rtable"

######################################################################################################
######################################################################################################
##################################################################
######################################################################################################
######################################################################################################

more "${VCF}.chr19.HETERO.SNP.BIALLELIC.18aug2018.vcf"
more "${VCF}.chr19.HETERO.SNP.BIALLELIC.18aug2018.vcf.rtable" 

######################################################################################################
######################################################################################################
############## after inspection on 18AUG2018, this code works OK : RSAMTOOLS
######################################################################################################
######################################################################################################

library("VariantAnnotation")
library("Rsamtools")
library("TitanCNA")

# vcf <-VariantAnnotation::readVcf("SPCG-OS040_7G.GRCh38p5M.MD_samtools.germline.HETERO.chr19.vcf.HET4.SNP.BIALLELIC.12june2017.vcf", 
#                                  genome="hg38")

vcf <- VariantAnnotation::readVcf("SPCG-OS040_7G.GRCh38p5M.MD_intervals_ALL.HAPLOTYPE_CALLER.ALL_CHRS.vcf.chr19.HETERO.SNP.BIALLELIC.18aug2018.vcf",
                                   genome="hg38")

vcf_granges <- rowRanges(vcf)

# vr <- readVcfAsVRanges("SPCG-OS040_7G.GRCh38p5M.MD_samtools.germline.HETERO.chr19.vcf.HET3.vcf.BI.ALLELIC.vcf", 
#                        genome="hg38")
# df <- as.data.frame(vr)
# sbp_vr = ScanBamParam(which=vr)

bam_file <- "SPCG-OS040_8T.GRCh38p5M.MD_intervals_chr19.IR.RC.bam"

# BAM_file <- BamFile(bam_file)
# seqinfo(BAM_file)
# ALN <- scanBam(BAM_file)
# names(ALN)
# quickBamFlagSummary(BAM_file)

bam_param = ScanBamParam(which=vcf_granges)

pileup_param = PileupParam(max_depth=1000, 
                   min_base_quality=13, 
                   min_mapq=0,
                   min_nucleotide_depth=1, 
                   min_minor_allele_depth=0,
                   distinguish_strands=FALSE, 
                   distinguish_nucleotides=TRUE,
                   ignore_query_Ns=TRUE, 
                   include_deletions=FALSE, 
                   include_insertions=FALSE)

snp_counts <- pileup(bam_file, scanBamParam=bam_param, pileupParam=pileup_param)

# head(snp_counts)
#  seqnames   pos nucleotide count         which_label
#1    chr19 61766          A    35 chr19:61766-61766.1
#2    chr19 61766          C     9 chr19:61766-61766.1
#3    chr19 61774          A    40 chr19:61774-61774.2
#4    chr19 61774          G     9 chr19:61774-61774.2
#5    chr19 61799          G     7 chr19:61799-61799.3
#6    chr19 61799          T    40 chr19:61799-61799.3

write.table(snp_counts, 
            file="SPCG-OS040_8T.GRCh38p5M.MD_intervals_chr19.IR.RC.bam _ table with allele counts in R v3",
            sep="\t", row.names=F, col.names=T)   

######################################################################################################
######################################################################################################
### if we would like to include a specific chromosome : -L chr19 \
######################################################################################################
######################################################################################################

$GATK \
-T SelectVariants \
-R $REFERENCE_HG38_MALE \
-L chr19 \
-select 'vc.getGenotype("SPCG-OS040_7G").isHet()' \
-selectType SNP \
-restrictAllelesTo BIALLELIC \
-V $VCF \
-o "${VCF}.chr19.HETERO.SNP.BIALLELIC.18aug2018.vcf"

$GATK -T ASEReadCounter \
-R $REFERENCE_HG38_MALE \
-U ALLOW_N_CIGAR_READS \
-I $BAM_TUMOR_chr19 \
-minDepth 2 \
--minBaseQuality 10 \
--minMappingQuality 10 \
-sites "${VCF}.chr19.HETERO.SNP.BIALLELIC.18aug2018.vcf" \
-o "${VCF}.chr19.HETERO.SNP.BIALLELIC.18aug2018.vcf.v2.rtable"

##################################################################

more "${VCF}.chr19.HETERO.SNP.BIALLELIC.18aug2018.vcf"
more "${VCF}.chr19.HETERO.SNP.BIALLELIC.18aug2018.vcf.rtable" 

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
