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
library(Biostrings)


ref <- readDNAStringSet("Ds0432.1/reference/ref.fa")

#names(ref)
#width(ref)

custom.genome <- toGRanges(data.frame(chr=c(names(ref)), start=c(1), end=c(width(ref))))


# Load raw snp file
rawsnpRNS07 <- read.table("RNS07.7.3B/snps.tab", sep="\t", header=T, quote="")
rawsnp13301B <- read.table("13301B/snps.tab", sep="\t", header=T, quote="")
rawsnp13311A <- read.table("13311A//snps.tab", sep="\t", header=T, quote="")
rawsnp13481A <- read.table("13481A/snps.tab", sep="\t", header=T, quote="")
rawsnp101051A <- read.table("101051A/snps.tab", sep="\t", header=T, quote="")
rawsnp151021A <- read.table("151021A/snps.tab", sep="\t", header=T, quote="")
rawsnpDs0432 <- read.table("Ds0432.1/snps.tab", sep="\t", header=T, quote="")
rawsnpIPO2222 <- read.table("IPO2222/snps.tab", sep="\t", header=T, quote="")
rawsnpPPO9019 <- read.table("PPO9019/snps.tab", sep="\t", header=T, quote="")
rawsnpPPO9134 <- read.table("PPO9134/snps.tab", sep="\t", header=T, quote="")
rawsnpRNS0512A <- read.table("RNS05.1.2A/snps.tab", sep="\t", header=T, quote="")
rawsnpA62 <- read.table("A623S20A17/snps.tab", sep="\t", header=T, quote="")

# convert raw snp file to GRanges file
snpRNS07 <- toGRanges(data.frame(chr=c(names(ref)), start=rawsnpRNS07$POS, end=rawsnpRNS07$POS))
snp13311A <- toGRanges(data.frame(chr=c(names(ref)), start=rawsnp13311A$POS, end=rawsnp13311A$POS))
snp13301B <- toGRanges(data.frame(chr=c(names(ref)), start=rawsnp13301B$POS, end=rawsnp13301B$POS))
snp13481A <- toGRanges(data.frame(chr=c(names(ref)), start=rawsnp13481A$POS, end=rawsnp13481A$POS))
snp101051A <- toGRanges(data.frame(chr=c(names(ref)), start=rawsnp101051A$POS, end=rawsnp101051A$POS))
snp151021A <- toGRanges(data.frame(chr=c(names(ref)), start=rawsnp151021A$POS, end=rawsnp151021A$POS))
snpDs0432 <- toGRanges(data.frame(chr=c(names(ref)), start=rawsnpDs0432$POS, end=rawsnpDs0432$POS))
snpIPO2222 <- toGRanges(data.frame(chr=c(names(ref)), start=rawsnpIPO2222$POS, end=rawsnpIPO2222$POS))
snpIPPO9019 <- toGRanges(data.frame(chr=c(names(ref)), start=rawsnpPPO9019$POS, end=rawsnpPPO9019$POS))
snpIPPO9134 <- toGRanges(data.frame(chr=c(names(ref)), start=rawsnpPPO9134$POS, end=rawsnpPPO9134$POS))
snpRNS0512A <- toGRanges(data.frame(chr=c(names(ref)), start=rawsnpRNS0512A$POS, end=rawsnpRNS0512A$POS))
snpA62 <- toGRanges(data.frame(chr=c(names(ref)), start=rawsnpA62$POS, end=rawsnpA62$POS))

# create a pdf file to store the figure
pdf(file="snp_plot_density.pdf")

# Plot the reference sequence
kp <- plotKaryotype(genome=custom.genome, plot.type =4, labels.plotter=NULL)

# Add unites on the reference sequence
kpAddBaseNumbers(kp, tick.dist = 1000000, tick.col="red", cex=0.5,minor.tick.dist = 500000, minor.tick.col = "gray", units ="Mb", add.units=TRUE)

# Plot SNP densities
ka<-kpPlotDensity(kp,data=snpIPO2222, window.size=1000,r0=0,r1=0.08, col="#32CD32", border="#32CD32")
kpAxis(kp, ymax=ka$latest.plot$computed.values$max.density, r0=0, r1=0.08, cex=0.25)
kpAddLabels(kp, labels="IPO2222", r0=0,r1=0.08, side="right", cex=0.35)

kb<-kpPlotDensity(kp,data=snpDs0432, window.size=1000,r0=0.09,r1=0.17, col="#32CD32", border="#32CD32")
kpAxis(kp, ymax=kb$latest.plot$computed.values$max.density, r0=0.09, r1=0.17, cex=0.25)
kpAddLabels(kp, labels="Ds0432-1", r0=0.09,r1=0.17, side="right", cex=0.35)

kc <-kpPlotDensity(kp,data=snp13301B, window.size=1000,r0=0.18,r1=0.23, col="#228B22", border="#228B22")
kpAxis(kp, ymax=kc$latest.plot$computed.values$max.density, r0=0.18, r1=0.23, cex=0.25)
kpAddLabels(kp, labels="13-30-1B", r0=0.18,r1=0.23, side="right", cex=0.35)

kd <-kpPlotDensity(kp,data=snp13311A, window.size=1000,r0=0.24,r1=0.32, col="#228B22", border="#228B22")
kpAxis(kp, ymax=kd$latest.plot$computed.values$max.density, r0=0.24,r1=0.32, cex=0.25)
kpAddLabels(kp, labels="13-31-1A", r0=0.24,r1=0.32, side="right", cex=0.35)

ke <-kpPlotDensity(kp,data=snp13481A, window.size=1000,r0=0.33,r1=0.41, col="#228B22", border="#228B22")
kpAxis(kp, ymax=ke$latest.plot$computed.values$max.density, r0=0.33,r1=0.41, cex=0.25)
kpAddLabels(kp, labels="13-48-1A", r0=0.33,r1=0.41, side="right", cex=0.35)

kf <-kpPlotDensity(kp,data=snp151021A, window.size=1000,r0=0.42,r1=0.5, col="#228B22", border="#228B22")
kpAxis(kp, ymax=kf$latest.plot$computed.values$max.density, r0=0.42,r1=0.5, cex=0.25)
kpAddLabels(kp, labels="15-102-1A", r0=0.42,r1=0.5, side="right", cex=0.35)

kg<-kpPlotDensity(kp,data=snpRNS07, window.size=1000,r0=0.51,r1=0.59, col="#228B22", border="#228B22")
kpAxis(kp, ymax=kg$latest.plot$computed.values$max.density, r0=0.51,r1=0.59, cex=0.25)
kpAddLabels(kp, labels="07.7.3B", r0=0.51,r1=0.59, side="right", cex=0.35)

kh <-kpPlotDensity(kp,data=snp101051A, window.size=1000,r0=0.6,r1=0.68, col="#FF4500", border="#FF4500")
kpAxis(kp, ymax=kh$latest.plot$computed.values$max.density, r0=0.6,r1=0.68, cex=0.25)
kpAddLabels(kp, labels="101-05-1A", r0=0.6,r1=0.68, side="right", cex=0.35)

ki<-kpPlotDensity(kp,data=snpIRNS0512A, window.size=1000,r0=0.69,r1=0.77, col="#FF4500", border="#FF4500")
kpAxis(kp, ymax=ki$latest.plot$computed.values$max.density, r0=0.69,r1=0.77, cex=0.25)
kpAddLabels(kp, labels="05.1.2A", r0=0.69,r1=0.77, side="right", cex=0.35)

kj<-kpPlotDensity(kp,data=snpA62, window.size=1000,r0=0.78,r1=0.83, col="#FF4500", border="#FF4500")
kpAxis(kp, ymax=kj$latest.plot$computed.values$max.density, r0=0.78,r1=0.83, cex=0.25)
kpAddLabels(kp, labels="A623.S20.A17", r0=0.78,r1=0.83, side="right", cex=0.35)

kk<-kpPlotDensity(kp,data=snpIPPO9019, window.size=1000,r0=0.84,r1=0.92, col="#0000FF", border="#0000FF")
kpAxis(kp, ymax=kk$latest.plot$computed.values$max.density, r0=0.84,r1=0.92, cex=0.25)
kpAddLabels(kp, labels="PPO9019", r0=0.84,r1=0.92, side="right", cex=0.35)

kl<-kpPlotDensity(kp,data=snpIPPO9134, window.size=1000,r0=0.93,r1=1.01, col="#0000FF", border="#0000FF")
kpAxis(kp, ymax=kl$latest.plot$computed.values$max.density, r0=0.93,r1=1.01, cex=0.3)
kpAddLabels(kp, labels="PPO9134", r0=0.93,r1=1.01, side="right", cex=0.35)

dev.off()
````

## Extact the genes that are affected by snps in the referecen strain :
````python
from Bio import SeqIO
gbk_file_2017V = open("Dsl3337_17.gbk","r")

# read the gbk file
recorde = SeqIO.read(gbk_file_2017V,"genbank")

# look for cds features and get infos
for feature in recorde.features:
    if feature.type == 'CDS':        
        product = feature.qualifiers["product"][0]
        id = feature.qualifiers["locus_tag"][0]
        prot_seq = feature.extract(recorde.seq).translate(table=11, cds=True)
        #rawtr = feature.qualifiers.get('translation',[0])
        #tr = str(rawtr).replace("[","").replace("]","").replace("'","")
        # make a dict with keys that have multiple values as list
        if id not in idDict:
             idDict[id] = list()
             idDict[id].extend([prot_seq,product])  
             
# output file to store fasta
fasta_of_liste_gene = open("uniqIDgeneAffected.fasta","w")

# open the the file that contains the gene IDs, one ids per line 
with open("genesAffected.tab") as file:

        # erad the lines
        geneid = file.readlines()
        
        # make a liste of the ids, and remove the the backends
        geneid = [ids.rstrip() for ids in geneid]
        #print(geneid)

        # loop over the ids and gets information from the dicetionary generated above
        for ID in geneid:

        

            
            id_seq = idDict[ID][0]
            id_product = idDict[ID][1]
            print("Gene : " + str(ID) + " was added.." + id_product)# this for checking whta happen
            
            # write the fasta file
            fasta_of_liste_gene.write(">" + ID + "#" + str(id_product) + "\n" + str(id_seq) + "\n")
````

 ## Raun ANI calculations and plot pylogeny with ANI matrix
 
 ```` bash
 # reformat boostraps values
 sed 's:\([0-9]\{2\}\:\):B=\1:g' Dsolani.tree
 
 
 # Prepare the genomes 
 ls data/*.fna | while read file; do echo "$file" >> liste.genomes;done
 
 # Run fastANi
 fastANI --ql liste.genome --rl liste.genomes -o aniout -t 56
 
 ````
 ````R
 library(phytools)
 library(reshape2)
 
 data <- read.table("aniout.txt", sep ="\t")
 tree <- read.tree("Dsolani.tree")
 
 
 # make 2 D matrix
 mat <- acast(data, V1~V2, value.var="V3")
 
 
 # rename row and col names
 gsub(".fna","",rownames(mat)) -> rnames
 gsub(".fna","",colnames(mat)) -> clnames
 rownames(mat) <- rnames
 colnames(mat) <- clnames
 
 # rename tip labels
 gsub(".fna","",tree$tip) -> tipsL
 tree$tip.label <- tipsL
 
 # create pdf file
 pdf(file="final_tree.pdf")
 
 # Plot the tree with matrix
 phylo.heatmap(tree, mat, fsize=0.5, pts=FALSE, lwd=1, colors=colors<-colorRampPalette(colors=c("green","red"))(100))
 
 # node BP value support
 nodelabels(tree$node.label,node=2:tree$Nnode+Ntip(tree), adj=c(1,-0.2),frame="none")
 
 dev.off()
 ````
 ## Blast the genes affected with snp
 
 ````bash
 # to make the headers shorter
 sed -i 's/_all_R1_001_cutadapt_(paired)_trimmed_(paired)//g' data/*.fa
 
 
 # make blast db
 cat data/Contig_27_S26.fa data/Contig_28_S27.fa data/Contig_36_S35.fa data/Contig_55_S53.fa data/Contig_8_S8.fa | makeblastdb -dbtype nucl -title DslStrains -out db/dslstrains -parse_seqids
 
 
 # run tblastx
 
 ````
 ## Plot phylogenomics tree using ggtree with meta data 
 ````R
pdf(file="tree.pdf", width=10)
p <- ggtree(tree) %<+% meta + geom_tiplab(size =4, align=TRUE, linesize=0.1, offset=0.005) +geom_tippoint(aes(color = coutry), size =1.5) +  scale_color_brewer(name="Country", palette ="Spectral",na.value="grey")+ geom_tiplab(aes(label= date), color ="blue", offset = 0.012, size =3,align=TRUE,linetype="blank") + geom_tiplab(aes(label= host), color ="red", offset = 0.014, size =3,align=TRUE,linetype="blank")+ geom_text2(aes(subset=(as.numeric(label)> 80), label=label),size =3, hjust =1, vjust =-1)
gheatmap(p, mat, offset=0.016, width = 0.2, legend_title="ANI value", font.size=2, colnames=FALSE)
dev.off()
````
