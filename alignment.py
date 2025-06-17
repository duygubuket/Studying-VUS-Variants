# Importing necessary libraries from BioPython
import Bio
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import SeqIO,  SearchIO
from Bio import Seq
from matplotlib import pylab
from Bio.Blast import NCBIWWW


import os
from collections import Counter


#read the protein seq fasta file of homo_sapiens ALDH3A2 gene
seq_file_read = SeqIO.read('homo_sapiens_ALDH3A2_seq.fasta', 'fasta')
seqfromfile = seq_file_read.seq

#the protein seq 
protein_seq = "MELEVRRVRQAFLSGRSRPLRFRLQQLEALRRMVQEREKDILTAIAADLCKSEFNVYSQEVITVLGEIDFMLENLPEWVTAKPVKKNVLTMLDEAYIQPQPLGVVLIIGAWNYPFVLTIQPLIGAIAAGNAVIIKPSELSENTAKILAKLLPQYLDQDLYIVINGGVEETTELLKQRFDHIFYTGNTAVGKIVMEAAAKHLTPVTLELGGKSPCYIDKDCDLDIVCRRITWGKYMNCGQTCIAPDYILCEASLQNQIVWKIKETVKEFYGENIKESPDYERIINLRHFKRILSLLEGQKIAFGGETDEATRYIAPTVLTDVDPKTKVMQEEIFGPILPIVPVKNVDEAINFINEREKPLALYVFSHNHKLIKRMIDETSSGGVTGNDVIMHFTLNSFPFGGVGSSGMGAYHGKHSFDTFSHQRPCLLKSLKREGANKLRYPPNSQSKVDWGKFFLLKRFNKEKLGLLLLTFLGIVAAVLVKAEYY"



def blast_query(protein_seq):
    #query this protein from NCBI
    result_handle = NCBIWWW.qblast('blastp', 'pdb', protein_seq)
    blast_qresult = SearchIO.read(result_handle, 'blast-xml')
    return blast_qresult[0:20]


#This function will make the vus sequences based on the aa changes in protein seq
def Protein_changer (protein, location, aa , replace):
  threetoone = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N', 
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', 
     'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}
  location -= 1
  if protein[location] == threetoone[aa]:
    VUS_Protein = protein[:location] + threetoone[replace] + protein[location + 1:]
  else:
    print(protein[location])
    return "Something is wrong"
  
  return VUS_Protein


#This function will return our vus variants in a new obtained file
def vus_file():
    new_file = open("VUS_protein.fasta", "w")

    Vus_dict = {"Variant #0000325290 (NC_000017.10:g.19552298T>C, ALDH3A2(NM_000382.2):c.14T>C)" : "Val5Ala",
                "Variant #0000255886 (NC_000017.10:g.19552403A>G, ALDH3A2(NM_000382.2):c.119A>G)" : "Asp40Gly",
                "Variant #0000262212 (NC_000017.10:g.19559858T>G, ALDH3A2(NM_000382.2):c.651T>G)" : "Asp217Glu",
                "Variant #0000560653 (NC_000017.10:g.19564581G>C, ALDH3A2(NM_000382.2):c.940G>C)" : "Ala314Pro",
                "Variant #0000560656 (NC_000017.10:g.19566694A>G, ALDH3A2(NM_000382.2):c.989A>G)" : "Glu330Gly",
                "Variant #0000807853 (NC_000017.10:g.19566720A>G, ALDH3A2(NM_000382.2):c.1015A>G)" : "Ile339Val",
                "Variant #0000256069 (NC_000017.10:g.19566759A>G, ALDH3A2(NM_000382.2):c.1054A>G)" : "Ile352Val",
                "Variant #0000325293 (NC_000017.10:g.19566768C>T, ALDH3A2(NM_000382.2):c.1063C>T)" : "Arg355Cys",
                "Variant #0000560657 (NC_000017.10:g.19566799C>T, ALDH3A2(NM_000382.2):c.1094C>T)" : "Ser365Leu"}

    for gene , variant in Vus_dict.items():

        number = ""

        for i in range(len(variant)):
            if variant[i].isdigit():
                number += variant[i]

        first = variant[:3]
        second = variant[-3:]
        number = int(number)

        code = Protein_changer(protein_seq,number,first,second)

        new_file.write(">" + gene + "\n" + code + "\n")

    new_file.close()

    return new_file


#This function is for implement multiple alignment, at first, the function look at the lengths, do proper changes if the lengths are not the same and later do multiple sequence alignment
def same_length_multiple_align(input_file):
    records = SeqIO.parse(input_file, 'fasta')
    records = list(records) # make a copy, otherwise our generator
                            # is exhausted after calculating maxlen
    maxlen = max(len(record.seq) for record in records)

    # pad sequences so that they all have the same length
    for record in records:
        if len(record.seq) != maxlen:
            sequence = str(record.seq).ljust(maxlen, '.')
            record.seq = Seq.Seq(sequence)
    assert all(len(record.seq) == maxlen for record in records)

    # write to temporary file and do alignment
    output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')

    return output_file



# Calculate the distance matrix
def dist_matrix(alignment):
    calculator = DistanceCalculator('identity')
    distMatrix = calculator.get_distance(alignment)
    #print(distMatrix)
    return distMatrix
    
#This function will return the phylogenetic tree that we constructed.
def construct_phylo_tree(distMatrix):
    # Create a DistanceTreeConstructor object
    constructor = DistanceTreeConstructor()
    # Construct the phlyogenetic tree using UPGMA algorithm
    UGMATree = constructor.upgma(distMatrix)
    # Construct the phlyogenetic tree using NJ algorithm
    NJTree = constructor.nj(distMatrix)

    # Draw the phlyogenetic tree
    Phylo.draw_ascii(UGMATree)
    # Draw the phlyogenetic tree using terminal
    Phylo.draw_ascii(NJTree)

 
#return the variants file
vus_file=vus_file()

#blast the ALDH3A gene
blast_results=blast_query(protein_seq)
print(blast_results)

#Read the sequences and align
alignment_file=same_length_multiple_align('blast1000.fasta')
alignment = AlignIO.read(alignment_file, "fasta")
distMatrix=dist_matrix(alignment)
print(construct_phylo_tree(distMatrix))


