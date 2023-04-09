#################################################################
#################################################################
################################################################# HMMcopy utils and binaries
# https://github.com/shahcompbio/HMMcopy
# https://github.com/shahcompbio/hmmcopy_utils
# mapCounter
# gcCounter
# readCounter
#################################################################
#################################################################
################################################################# folder BINS

HMMcopy_util_bin="./HMMcopy/bin/"

HMMcopy_util_bin_readCounter="./HMMcopy/bin/readCounter"
HMMcopy_util_bin_gcCounter="./HMMcopy/bin/gcCounter"
HMMcopy_util_bin_mapCounter="./HMMcopy/bin/mapCounter"

#################################################################
#################################################################
################################################################# folder UTILS
#################################################################
#################################################################
#################################################################

HMMcopy_util_wig=".//HMMcopy/util/bigwig"
HMMcopy_util_mappability="./HMMcopy/util/mappability"
HMMcopy_util_seg="./HMMcopy/util/seg"

HMMcopy_util_wig_bigwigtowig="./HMMcopy/util/bigwig/bigWigToWig"
HMMcopy_util_mappability_map="./HMMcopy/util/mappability/generateMap.pl"

#################################################################
#################################################################
################################################################# PRE-COMPUTED

### PRECOMPUTED these files for the GENOMES we are using :
### if we have a mapability file already: $HMMcopy/bin/mapCounter bw.file
### http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeMapability/
### if we DO NOT have a mappability file : and need to generate it from SCRATCH :
### we do : ./HMMcopy/util/mappability/generateMap.pl (the default seems to have been -w, --window 35)

./mappability/generateMap.pl -i ./genome-bowtie/genome -w 150 -o genome.fa.map.150.bw genome.fa
./mappability/generateMap.pl -i ./genome-bowtie/genome -w 100 -o genome.fa.map.100.bw genome.fa
./mappability/generateMap.pl -i ./genome-bowtie/genome -w 50 -o genome.fa.map.50.bw   genome.fa

#################################################################
#################################################################
#################################################################
### some scripts : ./script_mappability_100bp.sh
### fixedStep chrom=chr1 start=1 step=1 span=1
#################################################################
#################################################################
#################################################################

$HMMcopy_util_wig_bigwigtowig genome.fa.map.100.bw genome.fa.map.100.wig

### to make a wig file : ./bigwig/bigWigToWig genomeM.fa.map.bw genomeM.fa.map.wig

$HMMcopy_util_bin_mapCounter -w 1000 genome.fa.map.100.bw > genome.fa.map.100.bin.1kb.wig
$HMMcopy_util_bin_mapCounter -w 1000 genomeM.fa.map.100.bw > genomeM.fa.map.100.bin.1kb.wig
$HMMcopy_util_bin_mapCounter -w 1000 genomeF.fa.map.100.bw > genomeF.fa.map.100.bin.1kb.wig

$HMMcopy_util_bin_mapCounter -w 1000 genomeM.fa.map.50.bw > genomeM.fa.map.50.bin.1kb.wig
$HMMcopy_util_bin_mapCounter -w 1000 genomeF.fa.map.50.bw > genomeF.fa.map.50.bin.1kb.wig

$HMMcopy_util_bin_mapCounter -w 1000 genome.fa.map.150bp.bw > genome.fa.map.150bp.bin.1kb.wig
$HMMcopy_util_bin_mapCounter -w 1000 genomeM.fa.map.150bp.bw > genomeM.fa.map.150bp.bin.1kb.wig
$HMMcopy_util_bin_mapCounter -w 1000 genomeF.fa.map.150bp.bw > genomeF.fa.map.150bp.bin.1kb.wig

###########################################################################################
###########################################################################################
#################################################################
#################################################################
################################################################# PRE-COMPUTED 

$HMMcopy_util_bin_gcCounter -w 1000 genome.fa > genome.fa.gc.1kb.wig
$HMMcopy_util_bin_gcCounter -w 1000 genomeF.fa > genomeF.fa.gc.1kb.wig
$HMMcopy_util_bin_gcCounter -w 1000 genomeM.fa > genomeM.fa.gc.1kb.wig

###########################################################################################
###########################################################################################

## to use these files :

## genome.fa.gc.1kb.wig                (fixedStep chrom=chr1 start=1 step=1000 span=1000)
## genome.fa.map.150bp.bin.1kb.wig     (fixedStep chrom=chr1 start=1 step=1000 span=1000)

## genomeF.fa.gc.1kb.wig               (fixedStep chrom=chr1 start=1 step=1000 span=1000) 
## genomeF.fa.map.150bp.bin.1kb.wig    (fixedStep chrom=chr1 start=1 step=1000 span=1000)

## genomeM.fa.gc.1kb.wig               (fixedStep chrom=chr1 start=1 step=1000 span=1000)
## genomeM.fa.map.150bp.bin.1kb.wig    (fixedStep chrom=chr1 start=1 step=1000 span=1000)

###########################################################################################
###########################################################################################
###########################################################################################
### important, if we choose to work per CHROMOSOME, all the work to be done per CHROMOSOME
#################################################################
#################################################################
###########################################################################################
###########################################################################################

NORMAL="SPCG-OS040_7G.GRCh38p5M.MD_intervals_chr19.IR.RC.bam"
TUMOR="SPCG-OS040_8T.GRCh38p5M.MD_intervals_chr19.IR.RC.bam"

# the indexes of these files :

NORMAL="SPCG-OS040_7G.GRCh38p5M.MD.bam"
TUMOR="SPCG-OS040_8T.GRCh38p5M.MD.bam"

###########################################################################################
###########################################################################################
# samtools index \
# SPCG-OS040_7G.GRCh38p5M.MD_intervals_chr19.IR.RC.bam \
# SPCG-OS040_7G.GRCh38p5M.MD_intervals_chr19.IR.RC.bam.bai
###########################################################################################
###########################################################################################

###########################################################################################
###########################################################################################
# samtools index \
# SPCG-OS040_8T.GRCh38p5M.MD_intervals_chr19.IR.RC.bam \
# SPCG-OS040_8T.GRCh38p5M.MD_intervals_chr19.IR.RC.bam.bai
###########################################################################################
###########################################################################################

$SAMTOOLS index \
$NORMAL \
"${NORMAL}.bai"

$SAMTOOLS index \
$TUMOR \
"${TUMOR}.bai"

###########################################################################################
###########################################################################################
# $HMMcopy_util_bin_readCounter -w 1000 $NORMAL  > "${NORMAL%.bam}.1kb.wig"
# $HMMcopy_util_bin_readCounter -w 1000 $TUMOR  > "${TUMOR%.bam}.1kb.wig"
###########################################################################################
###########################################################################################

### to change the names of the output files :

$HMMcopy_util_bin_readCounter -w 1000 $NORMAL > "${NORMAL}.1kb.wig"
$HMMcopy_util_bin_readCounter -w 1000 $TUMOR  > "${TUMOR}.1kb.wig"

###########################################################################################
###########################################################################################
########################################################################################### HMMCOPY
########### to continue with HMMcopy code from the R section :
########### to use the following files 
###########################################################################################
###########################################################################################
## genome.fa.gc.1kb.wig
## genome.fa.map.150bp.bin.1kb.wig

## genomeF.fa.gc.1kb.wig
## genomeF.fa.map.150bp.bin.1kb.wig

## genomeM.fa.gc.1kb.wig
## genomeM.fa.map.150bp.bin.1kb.wig

## SPCG-OS040_7G.GRCh38p5M.MD_intervals_chr19.IR.RC.1kb.wig
## SPCG-OS040_8T.GRCh38p5M.MD_intervals_chr19.IR.RC.1kb.wig

# library(HMMcopy)
# rfile <- system.file("extdata", "normal.wig", package = "HMMcopy") ## germline
# gfile <- system.file("extdata", "gc.wig", package = "HMMcopy")     ## genome-gc
# mfile <- system.file("extdata", "map.wig", package = "HMMcopy")    ## genome-map
# tfile <- system.file("extdata", "tumour.wig", package = "HMMcopy") ## tumor

## ch19, to work with these files :

#  SPCG-OS040_7G.GRCh38p5M.MD_intervals_chr19.IR.RC.1kb.wig.extracting.chr19.wig
#  SPCG-OS040_7G.GRCh38p5M.MD_intervals_hg38_male.interval_list.IR.RC.1kb.wig      ## WHOLE GENOME
#  SPCG-OS040_7G.GRCh38p5M.MD_intervals_hg38_male.interval_list.IR.RC.1kb.wig.extracting.chr19.wig


#  SPCG-OS040_8T.GRCh38p5M.MD_intervals_chr19.IR.RC.1kb.wig.extracting.chr19.wig
#  SPCG-OS040_8T.GRCh38p5M.MD_intervals_hg38_male.interval_list.IR.RC.1kb.wig      ## WHOLE GENOME
#  SPCG-OS040_8T.GRCh38p5M.MD_intervals_hg38_male.interval_list.IR.RC.1kb.wig.extracting.chr19.wig

###########################################################################################
###########################################################################################
###########################################################################################

library(HMMcopy)
library(TitanCNA)
library(doMC)

### rfile <- "SPCG-OS040_7G.GRCh38p5M.MD_intervals_hg38_male.interval_list.IR.RC.1kb.wig"
rfile <- "SPCG-OS040_7G.GRCh38p5M.MD_intervals_chr19.IR.RC.1kb.wig.extracting.chr19.wig"
### in order to subset a chromosome : rfile["chr19"]

### tfile <- "SPCG-OS040_8T.GRCh38p5M.MD_intervals_hg38_male.interval_list.IR.RC.1kb.wig"
tfile <- "SPCG-OS040_8T.GRCh38p5M.MD_intervals_chr19.IR.RC.1kb.wig.extracting.chr19.wig"
### in order to subset a chromosome : tfile["chr19"]

### gfile <- "genome.fa.gc.1kb.wig"
gfile <- "genome.fa.gc.1kb.wig.extracting.chr19.wig"

###  mfile <- "genome.fa.map.150bp.bin.1kb.wig"
mfile < -"genome.fa.map.150bp.bin.1kb.wig.extracting.chr19.wig"

### mfile <- "hg38.M.map.from_KS.wig"


###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
########################################################################################### VISUALIZE THE NORMAL SAMPLE :

normal_reads <- wigsToRangedData(rfile, gfile, mfile)
normal_reads[1000:1010,]

### in order to subset a chromosome : normal_reads <- normal_reads["chr19"]

normal_copy <- correctReadcount(normal_reads)
normal_copy[1000:1010, ]

tumor_reads <- wigsToRangedData(tfile, gfile, mfile)
tumor_reads[1000:1010, ]
### in order to subset a chromosome : tumor_reads <- tumor_reads["chr19"]

tumor_copy <- correctReadcount(tumor_reads)
tumor_copy[1000:1010, ]

###########################################################################################
###########################################################################################
#### correctReadcount requires at least about 1000 bins to work properly.
###########################################################################################
########################################################################################### visualize the corrections : in GERMLINE sample :

#### to visualize THE EFFECT OF CORRECTION :

par(cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7, mar = c(4, 4, 2, 0.5))
png("06april_visualize_effect_correction_germline.png")
plotBias(normal_copy, pch = 20, cex = 0.5)
dev.off()

#### to visualize the CNA PROFILES :

to par(mar = c(4, 4, 2, 0))
png("06april_visualize_profiles_germline.png")
plotCorrection(normal_copy, pch = ".")
dev.off()

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
########################################################################################### visualize the corrections : in TUMOR sample :

#### to visualize THE EFFECT OF CORRECTION :

par(cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7, mar = c(4, 4, 2, 0.5))
png("06april_visualize_effect_correction_tumor.png")
plotBias(tumor_copy, pch = 20, cex = 0.5)
dev.off()

#### visualize the CNA PROFILES :

par(mar = c(4, 4, 2, 0))
png("06april_visualize_profiles_tumor.png")
plotCorrection(tumor_copy, pch = ".")
dev.off()

###########################################################################################
###########################################################################################
###########################################################################################

## In these files, COPY fields is used, which is simply the base 2 log values of cor.map.
## cor.gc Readcounts after the first GC correction step
## cor.map cor.gc readcounts after a furthur mappability correction
## copy cor.map transformed into log2 space

###########################################################################################
###########################################################################################
###########################################################################################
######################################################################
######################################################################
######################################################################
##################################### we are doing MATCHED TUMOR-NORMAL SAMPLE CORRECTION :
######################################################################
######################################################################
###########################################################################################
###########################################################################################
###########################################################################################

somatic_copy <- tumor_copy
head(somatic_copy,100)
somatic_copy$copy <- tumor_copy$copy - normal_copy$copy
head(somatic_copy,100)

### We can do the SAME SEGMENTATION and VISUALIZATION  as it follows:
### in order to get the parameters :

# param <- HMMsegment(somatic_copy, getparam = TRUE) # retrieve converged parameters via EM
# param$mu <- log(c(1, 1.4, 2, 2.7, 3, 4.5) / 2, 2)
# param$m <- param$mu
# somatic_segments <- HMMsegment(somatic_copy, param)
# there are some adjustments to increase or decrease the number of segments :
# https://raw.githubusercontent.com/Bioconductor/copy-number-analysis/master/inst/script/HMMcopy.R
# https://github.com/Bioconductor/copy-number-analysis/blob/master/inst/script/HMMcopy.R

### in the manual, it is : somatic_segments <- HMMsegment(somatic_copy, newmu_param, verbose = FALSE)

somatic_segments <- HMMsegment(somatic_copy, verbose = TRUE)

### VISUALIZATION :

par(cex.main = 0.5, cex.lab = 0.5, cex.axis = 0.5,
    mar = c(2, 1.5, 0, 0), mgp = c(1, 0.5, 0))

png("06april_visualize_after_HMMsegment_on_somatic_copy.png")
plotSegments(somatic_copy, somatic_segments, pch = ".",
             ylab = "Tumour Copy Number",
             xlab = "Chromosome Position")
dev.off()

######################################################
######################################################

### https://academic.oup.com/bioinformatics/article/26/24/3051/288761/CNAseg-a-novel-framework-for-identification-of
### export to SEG format for CNAseq segmentation, or for inspection of the data :

rangedDataToSeg(somatic_copy, file = "06april_for_CNAseg_tumor_corrected_copy.seg")
rangedDataToWig(somatic_copy, file = "06april_for_CNAseg_tumor_corrected_copy.wig")

##### VISUALIZATION :

png("06april_visualize_plotBias_somatic_copy.png")
plotBias(somatic_copy)
dev.off()

png("06april_visualize_plotCorrection_somatic_copy.png")
plotCorrection(somatic_copy)
dev.off()

png("06april_visualize_plotSegments_somatic_copy.png")
plotSegments(somatic_copy, somatic_segments)
dev.off()

###########################################################################################
########################################################################################### messages from the authors (DL) :

# util/mappability/generateMap.pl -o hg38.mappability.bw -w <READ_LENGTH> hg38.fa# 

# where -o specifies the output file, -w is the read length, and the input to the script is the reference genome.
# it created a 1 basepair resolution BigWig (.bw) mappability file specific to your read length.  
# this is the time consuming step.

# You can then use this BigWig mappability file as input into bin/mapCounter script as specified on the guide.  
# this should be fairly quick.
# when calling bin/mapCounter, you use the -w option to specify the “bin” width.

# As long as the -w option matches between your calls for gcCounter (on the hg38.fa), 
# mapCounter (on the .bw file you’re creating), and the readCounter (for your GERMLINE + TUMOR BAMs), 
# you should be fine.  Do you not need to create a unique .bw file for each bin size.  
# Changing window sizes should be relatively quick, so you can play around with readCounter to see what bin size you prefer.

###########################################################################################
########################################################################################### messages from the authors (DL) :

# The bin width depends on your depth of coverage.
# to test out a few bin width with the readCounter function to get bins of around 100 to 200 reads.
# In practice a 30x human genome does well with 1k windows, but if your coverage is high enough, 
# you can confidently lower the bin size even smaller to detect small-scale CNA events.

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
### to continue the pipeline for TitanCNA 
###########################################################################################
###########################################################################################
