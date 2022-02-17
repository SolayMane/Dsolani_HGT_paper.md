## We retreived all D solani(14-01-22)
[download_Dsolani_fna.sh](download_Dsolani_fna.sh)

## We run the phylogenomics pipeline
[BUSCO_phylogenomics.sh](BUSCO_phylogenomics.sh)

## We alsoe used orthofinder to cluster the proteomes and infere the Spceies Tree
````bash
orthofinder -t 56 -M msa -A mafft -T iqtree -f data/faa/
````
## Run Snippy on the strains





## Get stats about vatriants

````bash

echo -e "Strain\tNumVariant\tGenicVariants\tIntergenicVariant\tNumGeneInfected\tSynoV\tMissSens\tComplexType\tsnpType\tmnpType\tdelType\tinsType" > deepStat.tab


for file in $(find . -type f -name "snps.tab")
do
        strain=$(echo "$file" | cut -d"/" -f 2)
        NumVariant=$(cat $file | grep -v "CHROM" | wc -l)
        GenicVariants=$(cat $file | grep -v "CHROM" | cut -f 12 | grep "^D" | wc -l)
        IntergenicVariant=$(($NumVariant - $GenicVariants))
        NumGeneInfected=$(cat $file | grep -v "CHROM" | cut -f 12 | grep "^D" | sort | uniq | wc -l)
        SynoV=$(cat $file | grep -v "CHROM" | cut -f 11 | grep "^syn" | wc -l)
        MissSens=$(cat $file | grep -v "CHROM" | cut -f 11 | grep "^miss" | wc -l)
        ComplexType=$(cat $file | grep -v "CHROM" | cut -f 3 | grep "complex" | wc -l)
        snpType=$(cat $file | grep -v "CHROM" | cut -f 3 | grep "snp" | wc -l)
        mnpType=$(cat $file | grep -v "CHROM" | cut -f 3 | grep "mnp" | wc -l)
        delType=$(cat $file | grep -v "CHROM" | cut -f 3 | grep "del" | wc -l)
        insType=$(cat $file | grep -v "CHROM" | cut -f 3 | grep "ins" | wc -l)
        #genesCount=$(cat $file | grep -v "CHROM" | cut -f 12 | grep "^D" | sort | uniq)
        varPerGene=$(cat $file | grep -v "CHROM" | cut -f 12 | grep "^D" | sort | uniq -c)
        #totalgeneBP=$(cat $file | grep -v "CHROM" | cut -f 12 | grep "^D" | sort | uniq | grep -f - geneSize.tab |  awk '{ sum+= $2} END { if (NR > 0) print sum}')
        #AvrageVariantPerGene=$(cat $file | grep -v "CHROM" | cut -f 12 | grep "^D" | sort | uniq -c | awk -v var="$totalgeneBP" '{ sum+= $1} END { print (sum/var)/1000 }')

        echo -e "$strain\t$NumVariant\t$GenicVariants\t$IntergenicVariant\t$NumGeneInfected\t$SynoV\t$MissSens\t$ComplexType\t$snpType\t$mnpType\t$delType\t$insType" >>deepStat.tab
        echo -e "$varPerGene" > ${strain}.varPerGene.txt
        
        
````

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
````

## Extact the genes that are affected by snps in the referecen strain :
````python
from Bio import SeqIO

def get_annotation_with_ID(gbk_file, annoType, qualifier, value):

        #read the gbk file
        recorde = SeqIO.read(gbk_file,"genbank")

        #loop over the features
        for feature in recorde.features:
                if feature.type == annoType and value in feature.qualifiers.get(qualifier, []):
                        return feature.extract(recorde.seq).translate(table=11,cds=True)
        return None





fasta_of_liste_gene = open("liste.gene.fasta","w")

with open("liste.gene") as file:
        geneid = file.readlines()
        geneid = [ids.rstrip() for ids in geneid]

        for id in geneid:

                translation = get_annotation_with_ID("Dsl3337_17.gbk","CDS","locus_tag", id)
                fasta_of_liste_gene.write(">" + id + "\n" + str(translation) + "\n")
                print ("Gene number" + id +" was added")

