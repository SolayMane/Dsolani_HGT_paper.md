## We will analyse functionally the hgt regions
````bash
#extrcat the genes that intersect with the hgt regions
#one region could intersects with several genes!
bedtools intersect -a Dsl3337.gff -wa -wb -b reg.bed > inter.tab

# get the liste of the genes
awk -F"\t" '{ print $9}' inter.tab  | cut -d"=" -f3 > liste.genes

# get the genes in aa that intersect with hgt regions using python tool
./getGenes_aa.py -h


````

