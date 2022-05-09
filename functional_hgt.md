## We will analyse functionally the hgt regions
````bash
#extrcat the genes that intersect with the hgt regions
#one region could intersects with several genes!
bedtools intersect -a Dsl3337.gff -wa -wb -b reg.bed > inter.tab

# get the liste of the genes
awk -F"\t" '{ print $9}' inter.tab  | cut -d"=" -f3 > liste.genes

# get the genes in aa that intersect with hgt regions using python tool
./getGenes_aa.py -h

usage: getGenes_aa.py [-h] -g GBK [-l LISTE] [-o OUT]

Tools to retreive the aa from gbk given a liste of gene IDs

optional arguments:
  -h, --help            show this help message and exit
  -g GBK, --gbk GBK     Input gbk file (default: None)
  -l LISTE, --liste LISTE
                        Input the gene liest (default: None)
  -o OUT, --out OUT     Output file in fasta format (default: None)

# run the script
[./getGenes_aa.py](/getGenes_aa.py) -g Dsl3337_17.gbk -o genes_in_hgt.fasta -l liste.genes 
````

