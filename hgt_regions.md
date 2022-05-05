### This steps decribe how we analyzed the regions with variants hotspot
- we performed a snippy analysie on the whole genome sequences of D. solani
- We obtained  an alignment of "core SNPs" called core SNP genome of all strains based on the reference Dsl 3337
- We have calculated the SNP density for the five strains include in this study (window size = 1000 bp)
- We extracted all the windows that have SNP density > 9 varinat/1Kbp
- The coordinates were used to extract those regions from the Core SNP genome alignment using extractalign from EMBOSS
- The phylogeneies were infered from those subalignments
- We used ete3 to root the trees using D. dadantii as outgroup 
