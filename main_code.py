#REFERENCES 
#Gado, Japheth E., 2021, "Machine learning and bioinformatic insights into key enzymes for a bio-based circular economy". Theses and Dissertations:
# Chemical and Materials Engineering, 129. DOI: 10.13023/etd.2021.009.


#This code is adopted from the above article and modified. 


from pycanal import Canal

# Create an instance of the Canal class
canal = Canal(fastafile='aligned_new.fasta', #Multiple sequence alignment (MSA) of homologous sequences
              ref=0, #Position of reference sequence in MSA, use first sequence
              startcount=-16, #Position label of first residue in reference sequence
              verbose=True # Print out progress
              )

# Compute conservation scores for each site in reference sequence with relative entropy method
cons_scores = canal.analysis(include=None, method='relative')
print(cons_scores)

# Plot the distribution of amino acids in at position 77 and save image as position77.png
canal.plotSiteDistribution(site=352, saveplot='position352')

# Determine consensus sequence from the alignment and save as fasta file
consensus_sequence = canal.getConsensusSequence(savefasta='consensus_sequence.fasta')