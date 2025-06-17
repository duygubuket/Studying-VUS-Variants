def fastareader(filename):
    seqDict = {}
    with open(filename, "r") as filein:
        for line in filein:
            if line.startswith(">"):
                header = line[1:].strip()
                seqDict[header] = ""
            else:
                seqDict[header] += line.strip()
    return seqDict

#Saves aa and their rates throughout all sequences in dictionary

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

#Finds exact sequences that has selected aa in selected position

def find(position, dict , whichaa):
    position = position - 1
    selected = []
    for header,sequence in dict.items():
        aa = sequence[position]
        if aa == whichaa:
            selected.append(header)
    return selected


#Opens a fasta file, reads aligned fasta and saves rates of aa in entered positions

rates = open('ratesofaa.faa', 'w')
dict_of_blast = fastareader("no_empty_sude_blast_aligned.fas")
inputs = []

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

rates.close()

# Used to find special sequences that have VUS aa change for the comparison using Neighborjoining tree
print(find(355,dict_of_blast,"C"))