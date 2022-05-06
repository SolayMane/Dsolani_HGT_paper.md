### This steps decribe how we analyzed the regions with variants hotspot
- we performed a snippy analysie on the whole genome sequences of D. solani
- We obtained  an alignment of "core SNPs" called core SNP genome of all strains based on the reference Dsl 3337
- We have calculated the SNP density for the five strains include in this study (window size = 1000 bp)
- We extracted all the windows that have SNP density > 9 varinat/1Kbp
- The coordinates were used to extract those regions from the Core SNP genome alignment using extractalign from EMBOSS
- The phylogeneies were infered from those subalignments
- We used ete3 to root the trees using D. dadantii as outgroup 
### Here is the pipeline used
````bash
#This script aims to extract subalignement of the regions containing SNP hotspots from the whole genomes alignment of 40 sequences of D solani strains btained using snippy software
# The subalignment will be used to infere trees

# create a folder where to store the sequences
mkdir Regions_rooted

# rename the reference sequence
sed 's/Reference/Dickeya_solani_RNS08.23.3.1A/' core.full.aln > core.full.aln.fasta

# get the locations of the regions from the bed file and read is as an array
cat reg.txt | awk -v region="1" '{print $1"-"$2"#""region"region++}' | while IFS="#" read -a reg

# use extractalign tool from EMBOSS
do extractalign -sequence core.full.aln.fasta -regions ${reg[0]}  -outseq Regions_rooted/${reg[1]}.fasta
done

#infere the tree from the subalignments using nextstrain pipeline
for alignment in Regions_rooted/*.fasta

do
name=$(basename $alignment)
FastTree -nt -boot 1000 -gtr $alignment > Regions_rooted/${name%%.*}.tree

#augur tree -a $alignment --method fastree --tree-builder-args="--root-seq $alignment,Dickeya_dadantii_3937" -o Regions_rooted/${name%%.*}.tree.nwk --nthreads 56

done
````
````python
import ete3
from ete3 import TextFace, Tree,faces, AttrFace, TreeStyle, NodeStyle
import sys
import glob


#read the file
#inputdire = sys.argv[1]

# read the tree files
mytrees = [f for f in glob.glob("Regions_rooted/*.tree")]


#print(mytrees)

def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=10)
        faces.add_face_to_node(N, node, 0)
#        F = AttrFace("support", fsize=30)
#        faces.add_face_to_node(F, node, 0, position="branch-top")


# iterate over the tree file in the list
for tree in mytrees:

    #read the tree file
    t = Tree(tree)
  #set the outgroup
    t.set_outgroup("Dickeya_dadantii_3937")

    # write the tree to new file
#    st = NodeStyle()
#    st["vt_line_width"] = 5
#    st["hz_line_width"] = 5
#    st["size"] = 6
    nst1 = NodeStyle()

    nst1["bgcolor"] = "LightSteelBlue"
#    nst1["vt_line_width"] = 8
#    nst1["hz_line_width"] = 8
    nst2 = NodeStyle()
    nst2["bgcolor"] = "Moccasin"
#    nst2["vt_line_width"] = 8
#    nst2["hz_line_width"] = 8
#    for n in t.traverse():
#        n.dist = 0
    n1 = t.get_common_ancestor("Dickeya_solani_Am3a","Dickeya_solani_CH05026","Dickeya_solani_CH07044","Dickeya_solani_CH9635-1","Dickeya_solani_CH9918-774","Dickeya_solani_D12","Dickeya_solani_Ds0432-1","Dickeya_solani_F012","Dickeya_solani_GBBC2040","Dickeya_solani_IFB0099","Dickeya_solani_IFB0167","Dickeya_solani_IFB0212","Dickeya_solani_IFB0223","Dickeya_solani_IFB0231","Dickeya_solani_IFB0311","Dickeya_solani_IFB0417","Dickeya_solani_IFB0421","Dickeya_solani_IFB0487","Dickeya_solani_IFB0695","Dickeya_solani_IFB_0158","Dickeya_solani_IFB_0221","Dickeya_solani_IPO2019","Dickeya_solani_IPO_2222","Dickeya_solani_M21a","Dickeya_solani_MIE35","Dickeya_solani_MK10","Dickeya_solani_MK16_MK16","Dickeya_solani_PPO9019","Dickeya_solani_PPO9134","Dickeya_solani_RNS10-27-2A","Dickeya_solani_RNS08.23.3.1A","Dickeya_solani_Sp1a")
    n1.set_style(nst1)
        n2 = t.get_common_ancestor("101051A","13301B","13311A","13481A","151021A","Dickeya_solani_A623-S20-A17","RNS05.1.2A","Dickeya_solani_RNS0773B")
    n2.set_style(nst2)
#    n2.add_face(TextFace("0512 sub-clade "), column=2, position = "branch-right")
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False
    ts.show_branch_support = True
    ts.mode = "r"
    ts.root_opening_factor = 1



#    t.render(format=1, outfile= tree + ".newick")
    t.render(tree+ ".png", w=400, units="px",tree_style=ts)
````


