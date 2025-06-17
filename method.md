In this code, the detailed information of the codes in main_code.py could be found:

At the beginning we tried to analyze the nucleotide sequence and then the protein sequences and vus variants.
The "homo_sapiens_ALDH3A2_seq.fasta" has the nucleotide seq for ALDH3A2 gene in Homo sapiens.

After we read it from file and transcribe and translate it into protein sequence, we looked at the aa freqeuencies in our sequence.

##alignment.py

We did blast search (blastp) to compare our protein sequence with a database of protein sequences in the blast_query function.

```ruby
def blast_query(protein_seq):
    #query this protein from NCBI
    result_handle = NCBIWWW.qblast('blastp', 'pdb', protein_seq)
    blast_qresult = SearchIO.read(result_handle, 'blast-xml')
    return blast_qresult[0:5]
```

The blast results are shown in below:
Program: blastp (2.13.0+)
  Query: unnamed (485)
         protein product
 Target: pdb
   Hits: ----  -----  ----------------------------------------------------------
            #  # HSP  ID + description
         ----  -----  ----------------------------------------------------------
            0      1  pdb|4QGK|A  Structure of the Human Sjogren Larsson Synd...
            1      1  pdb|3SZA|A  Crystal structure of human ALDH3A1 - apo fo...
            2      1  pdb|1AD3|A  CLASS 3 ALDEHYDE DEHYDROGENASE COMPLEX WITH...
            3      1  pdb|6K0Z|A  Chain A, Aldehyde dehydrogenase [Staphyloco...
            4      1  pdb|5MYP|A  Structure of apo-TbALDH3 [Trypanosoma brucei]

The Protein_changer function made the proper substitution changes.

```ruby
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

```


The vus_file function will return the variants protein sequences that we obtained from the previous Protein_change function in a new file called "VUS_proteins.fasta"

```ruby
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

```

The same_length_multiple_alig function will return the aligned homologs protein sequences, we read the unaligned seqeunces from "blast1000.fas, than check the lengths for multiple sequence alignments and if they are not in the same length, proper changes were made. After that multiple sequnce alignment took place and aligned sequences were recorded in the new file called "blast1000_padded.fas"

```ruby
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

```

The dist_matrix function is calculated the distance matrix adn than return it.
```ruby
def dist_matrix(alignment):
    calculator = DistanceCalculator('identity')
    distMatrix = calculator.get_distance(alignment)
    #print(distMatrix)
    return distMatrix

```


The construct_phylo_tree will construct the phylogenetic tree based on UPGMA and Neighbor joining tree methods. 
```ruby
  def construct_phylo_tree(distMatrix):
    # Create a DistanceTreeConstructor object
    constructor = DistanceTreeConstructor()
    # Construct the phlyogenetic tree using UPGMA algorithm
    UGMATree = constructor.upgma(distMatrix)
    # Construct the phlyogenetic tree using NJ algorithm
    NJTree = constructor.nj(distMatrix)

    # Draw the phlyogenetic tree
    Phylo.draw(UGMATree)
    # Draw the phlyogenetic tree using terminal
    Phylo.draw_ascii(NJTree)

```

<<<<<<< HEAD
#blast1000.fasta

The fasta file obtained from blastp, 1000 hits.

#aligned_blast1000.fas

Aligned homologs fasta file.

The constructed phylogenetics trees are recorded in the files were printed.


The utils.py include the utils functions such as reading a fasta file and writing in a fasta file. 

### scoringchanges.py

To make comparison in cases of the changes in aminoacid bases and DNA bases, a code that scores transition and transversion for bases and aminoacid chemical property changes is written.

```ruby
  def check_transition (previous_base, new_base):
    if previous_base == "T" or previous_base == "C":
      if new_base == "T" or new_base == "C":
        return True
    else:
      if new_base == "G" or new_base == "A":
        return True
    return False
```

```ruby
  def chemical (previousaa, newaa):
    chemical = {"nonpolar": ["G","A","V","C","P","L","I","M","W","F"], "polar" : ["S","T","Y","N","Q"], "neg" : ["D","E"], "pos" : ["K","R","H"]}
    if previousaa in chemical["nonpolar"]:
      if newaa in chemical["nonpolar"]:
        return True
    elif previousaa in chemical["polar"]:
      if newaa in chemical["polar"]:
        return True
    elif previousaa in chemical["neg"]:
      if newaa in chemical["neg"]:
        return True
    elif previousaa in chemical["pos"]:
      if newaa in chemical["pos"]:
        return True
    return False
```
### rates.py

A function is written that takes dictionary of blast results and position of interest and counts the amount of amino acids. 
It sends the found amino acids with their rate in whole sequences

```ruby
  def checker(position,dict):
      position = position - 1
      aminoacids = {}
      total = 0
      for header,sequence in dict.items():
          aa = sequence[position]
          if aa == "-":
              aa = "empty"
          if aa in aminoacids.keys():
              aminoacids[aa] += 1
          else:
              aminoacids[aa] = 1
          total += 1
      
      for i in aminoacids.keys():
          aminoacids[i] = round((int(aminoacids[i])/total)*100,2)
      return aminoacids
```
A secondary finder is used to pinpoint exact sequences that has VUS in their genome.

```ruby
  def find(position, dict , whichaa):
      position = position - 1
      selected = []
      for header,sequence in dict.items():
          aa = sequence[position]
          if aa == whichaa:
              selected.append(header)
      return selected
```

Lastly, using loops, each frequency for each input of position is detected and written in a new file.

```ruby
  while True:

      user_input = input("Enter position or quit: ")
      if user_input == "quit":
          break
      inputs.append(int(user_input))
  for i in inputs:
      aawithrate = checker(i,dict_of_blast)
      rates.write(">" + str(i) + "\n")
      for aa,rate in aawithrate.items():
          rates.write(aa + '\t' + str(rate) + '\n')
      rates.write("\n")

```