## We retreived all D solani(14-01-22)
[download_Dsolani_fna.sh](download_Dsolani_fna.sh)

## We run the phylogenomics pipeline
[BUSCO_phylogenomics.sh](BUSCO_phylogenomics.sh)


## Run Snippy on the strains


## Plot the snp using

````R
library(karyoploteR)
library(GenomicRanges)
library(seqinr)
library(Biostrinsg)# use BiocManager to install the package if not installed


# read the genome file
ref <- readDNAStringSet("Dsl3337_17.fasta")

# show the names of the chromosomes :
names(ref)
# show the length of every chromosome
width(ref)
# so we will use theses data to build a custom genome for our sequence using toGranges

# to rename the sequence
names(ref) <-c("D. solani RNS 08.23.3.1.A")

custom.genome <- toGRanges(data.frame(chr=c(names(ref)), start=c(1), end=c(width(ref))))


# read the snp tab file output of snippy for evry strain
rawsnp13481A <- read.table("13-48.1A.snp.tab", sep="\t", header=T)

# converte the raw snp to GRanges object
snp13481A <- toGRanges(data.frame(chr="D. solani RNS 08.23.3.1.A", start=rawsnp13481A$POS, end=rawsnp13481A$POS+1)

# plot the chromosome
kp <- plotKaryotype(genome=custom.genome, plot.type =2)

# plot the snp density for each strain (SNP here is for the strain 13-31...)
kpPlotDensity(kp, data=SNP, window.size=1000, r0=0,r1=0.45, col="#FFAACB")
kpPlotDensity(kp, data=snp13481A, window.size=1000, r0=0.55,r1=1, col="red")


