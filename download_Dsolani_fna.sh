# check if the file summary file exist
if [ -s assembly_summary.txt ];then
        echo "The file is already there"
else
        echo -e " --\n-- Downloading the assembly summary file --\n--" && \
        curl https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -o assembly_summary.txt
fi

# get the ftp path from the table and also the strain name
count=0
awk -F"\t" ' ( $8 ~ /Dickeya solani/) {print $8";"$9";"$20";"$12}' assembly_summary.txt | while IFS=";" read -a field


do

        strain=$(echo "${field[1]##*=}" | sed 's/ //g' | sed 's/\///g')
        species=$(echo "${field[0]}"| sed 's/ /_/g' | sed 's/\///g')
        #echo $output
        #echo $species
       level=$(echo "${field[3]}"| sed 's/ /_/g')
       tcount=$(awk -F"\t" ' ( $8 ~ /Dickeya solani/) {print $12}' assembly_summary.txt | wc -l)
        ((count++))
       echo -e "Donwloading ${species}_${strain} $count/$tcount files ..... \n"
       echo -n "-----------------"
       curl --request GET "${field[2]}/${field[2]##*/}_genomic.fna.gz" > ${species}_${strain}_${level}.fna.gz
       gunzip ${species}_${strain}_${level}.fna.gz
done
# check if the file summaru file exist
if [ -s assembly_summary.txt ];then
        echo "The file is already there"
else
        echo -e " --\n-- Downloading the assembly summary file --\n--" && \
        curl https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -o assembly_summary.txt
fi

# get the ftp path from the table and also the strain name
count=0
awk -F"\t" ' ( $8 ~ /Dickeya solani/) {print $8";"$9";"$20";"$12}' assembly_summary.txt | while IFS=";" read -a field


do

        strain=$(echo "${field[1]##*=}" | sed 's/ //g' | sed 's/\///g')
        species=$(echo "${field[0]}"| sed 's/ /_/g' | sed 's/\///g')
        #echo $output
        #echo $species
       level=$(echo "${field[3]}"| sed 's/ /_/g')
       tcount=$(awk -F"\t" ' ( $8 ~ /Dickeya solani/) {print $12}' assembly_summary.txt | wc -l)
        ((count++))
       echo -e "Donwloading ${species}_${strain} $count/$tcount files ..... \n"
       echo -n "-----------------"
       curl --request GET "${field[2]}/${field[2]##*/}_genomic.fna.gz" > ${species}_${strain}_${level}.fna.gz
       gunzip ${species}_${strain}_${level}.fna.gz
done
