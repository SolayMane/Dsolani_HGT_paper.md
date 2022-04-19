## Download the genome sequences for all D. solani plus 2 D. dadantii and 2 D. dianthicola including type strains:
## Run KSNP3 as fellow :
````bash
# create input file contianing the paths for the genomes sequences
for file in data/*.fasta; do name=$(basename $file);echo -e "$(pwd)/$file\t${name%%.*}";done >> input.txt
# Run Kchooser to determine the best kmer to use
MakeFasta input.txt in.fas
Kchooser in.fasta
$ cat Kchooser.report
Initial value of k is 13.
When k is 13 0.818904598267137 of the kmers from the median length sequence are unique.
When k is 15 0.966633122511766 of the kmers from the median length sequence are unique.
When k is 17 0.992283005184264 of the kmers from the median length sequence are unique.
The optimum value of K is 19.
When k is 19 0.995847056367441 of the kmers from the median length sequence are unique.

There were 43 genomes.
The median length genome was 4879104 bases.
The time used was 802 seconds

From a sample of 995 unique kmers  239 are core kmers.
0.240201005025126 of the kmers are present in all genomes.

FCK = 0.240.
 # we run Kchooser again
 Kchooser in.fasta 0.995
 $ cat Kchooser.report
 Initial value of k is 13.
When k is 13 0.818904598267137 of the kmers from the median length sequence are unique.
When k is 15 0.966633122511766 of the kmers from the median length sequence are unique.
When k is 17 0.992283005184264 of the kmers from the median length sequence are unique.
When k is 19 0.996029984406822 of the kmers from the median length sequence are unique.
The optimum value of K is 21.
When k is 21 0.99660221536075 of the kmers from the median length sequence are unique.

There were 43 genomes.
The median length genome was 4879104 bases.
The time used was 829 seconds

From a sample of 997 unique kmers  228 are core kmers.
0.228686058174524 of the kmers are present in all genomes.

FCK = 0.229.
# so we will select K = 21 as length
# Run kSNP3 as fellow :
kSNP3 -in input.txt -outdir run1 -k 21 -ML -all_annotations
