#!/bin/bash

# create a work direcotry
if [ -d "phylo" ]
then
    rm -r phylo
else
    mkdir phylo
fi

cd phylo

#check if the output directory for storing busco runs exists if not create it
if [ -d "outBusco" ]
then
    rm -r outBusco/*
else
    mkdir outBusco
fi


cd ../
# running the busco on the genome sequence
for file in data/*.fna
do
        output=$(basename ${file%%.fna*})


        busco -i $file -l enterobacterales_odb10 -o $output -m genome --out_path phylo/outBusco --cpu 56
done

cd phylo
#check the folder where to move busco results for phylogenomics analysis
if [ -d "phyloInput" ]
then
    rm -r phyloInput/*
else
    mkdir phyloInput
fi

find outBusco/ -type d -name "run_*" | while read path
do
        new_name=$(basename ${path%/*})
        cp -r $path phyloInput/run_${new_name}
done


#run the phylogenomics pipeline on the busco runs and infere the tree

#check if the output folder exists and remove it . the pipeline will create a new one!
if [ -d "phylOutput" ]
then
    rm -r phylOutput && echo "The phylogenomics pipeline will start..."
else
    echo "The phylogenomics pipeline will start..."
fi

cd ../

python BUSCO_phylogenomics/BUSCO_phylogenomics.py -t 56 -l enterobacterales_odb10 --supermatrix -o phylo/phylOutput -d phylo/phyloInput/

