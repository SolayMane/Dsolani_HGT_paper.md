#!/usr/bin/env python3
#Created on 09-05-2022
#Author: Slimane Khayi; slimane.khayi@inra.ma



from Bio import SeqIO
import sys
import argparse

def main():
    """get the command line arguments """
    parser = argparse.ArgumentParser(
	description="Tools to retreive the aa from gbk given a liste of gene IDs",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-g","--gbk",help="Input gbk file", required=True)
    parser.add_argument("-l","--liste", help="Input the gene liest")
    parser.add_argument("-o","--out",help="Output file in fasta format")

    return parser.parse_args()

if __name__ == "__main__":
    main()



#inliste = sys.argv[1]
#ingbk = sys.argv[0]
args = main()


gbk_file_2017V = open(args.gbk,"r")

# read the gbk file
recorde = SeqIO.read(gbk_file_2017V,"genbank")


# a dict to store info
idDict ={}

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
fasta_of_liste_gene = open(args.out,"w")

# open the the file that contains the gene IDs, one ids per line 
with open(args.liste) as file:

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

