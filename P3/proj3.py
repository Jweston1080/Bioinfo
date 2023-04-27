
import time 
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
#------------Task 1 ------------------#
start_time = time.time()
print(f"----------------Task 1 start ----------------")
# function to calculate Hamming distance between two sequences
def hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2)
    return sum(s1 != s2 for s1, s2 in zip(seq1, seq2))

# read in fasta file
with open("input_file.fasta", "r") as f:
    seq_records = f.readlines()

# create list of sequences from fasta file
sequences = []
current_seq = ""
for line in seq_records:
    if line.startswith(">"):
        if current_seq:
            sequences.append(current_seq)
            current_seq = ""
    else:
        current_seq += line.strip()
if current_seq:
    sequences.append(current_seq)

# initialize empty distance matrix
num_seqs = len(sequences)
distance_matrix = [[0] * num_seqs for _ in range(num_seqs)]

# iterate through pairs of sequences and calculate Hamming distance
for i in range(num_seqs):
    for j in range(i+1, num_seqs):
        dist = hamming_distance(sequences[i], sequences[j])
        distance_matrix[i][j] = dist
        distance_matrix[j][i] = dist

# write distance matrix to output file in CSV format
with open("output_file.fasta", "w") as f:
    for row in distance_matrix:
        f.write(",".join(str(x) for x in row) + "\n")

end_time = time.time()
task1_time = end_time - start_time
print(f"----------------Task 1 runtime: {task1_time:.2f} seconds----------------")
#------------Task 2 ------------------#

start_time = time.time()
print(f"----------------Task 2 start ----------------")
# Replace "sequences.fasta" with the path to your FASTA file
records = AlignIO.read("input_file.fasta", "fasta")

# Create a MultipleSeqAlignment object
alignment = MultipleSeqAlignment(records)


# Calculate pairwise distances using the identity metric
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)

# Construct the phylogenetic tree
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# Draw and show the tree using Phylo module
Phylo.draw(tree)


end_time = time.time()
task2_time = end_time - start_time
print(f"----------------Task 2 runtime: {task2_time:.2f} seconds----------------")

#------------Task 3 ------------------#
start_time = time.time()
print(f"----------------Task 3 start ----------------")
records = AlignIO.read("input_file.fasta", "fasta")

# Create a MultipleSeqAlignment object
alignment = MultipleSeqAlignment(records)

# Calculate pairwise distances using the identity metric
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)

# Construct the phylogenetic tree
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# Get the ancestral sequences at each node using Phylo module
for clade in tree.find_clades(order='postorder'):
    if clade.is_terminal():
        seq_record = next((record for record in records if record.id == clade.name), None)
        clade.sequence = seq_record.seq
    else:
        seqs = [child.sequence for child in clade.clades]
        combined_seq = SeqRecord(Seq(''))
        for i in range(len(seqs[0])):
            most_common = max([seq[i] for seq in seqs], key=[seq[i] for seq in seqs].count)
            combined_seq.seq += Seq(most_common)
        clade.sequence = combined_seq.seq
    print(f"Node {clade.name}: {clade.sequence}")

end_time = time.time()
task3_time = end_time - start_time
print(f"----------------Task 3 runtime: {task3_time:.2f} seconds----------------")
#------------Task 4 ------------------#
start_time = time.time()
print(f"----------------Task 4 start ----------------")
# genomic sequence files of multiple coronaviruses

# Task 2 from P1 : Frequent words by hashing to identify 9 mers. 
def FrequencyMap(Word, k):
  freq_map = {}
  location_map={}
  n = len(Word)
  for i in range(n - k + 1):
      pattern = tuple(Word[i:i+k])
      if pattern in freq_map:
          freq_map[pattern] += 1
          location_map[pattern].append([i+1,i+k+1])
      else:
          freq_map[pattern] = 1
          location_map[pattern]=[[i+1,i+k+1]]
  return freq_map, location_map
   
def BetterFrequentWords(Word, k):
  freq_map = FrequencyMap(Word, k)[0]
  location_map=FrequencyMap(Word, k)[1]
  max_count = max(freq_map.values())
  frequent_patterns = [[pattern,location_map[pattern]] for pattern, count in freq_map.items() if count == max_count]
  return frequent_patterns[0][0]


with open ('data/COVID2.fasta') as f:
 genome=f.read().strip()
genome=genome[:-33] 

with open ('data/Bat.fasta') as f:
 genomebat=f.read().strip()
genomebat=genomebat[:-27] 

with open ('data/Middle_east.fasta') as f:
 genomeme=f.read().strip()

with open ('data/Pangolin.fasta') as f:
 genomep=f.read().strip()


with open ('data/SARS.fasta') as f:
 genomes=f.read().strip()

with open ('data/Severe_acute.fasta') as f:
 genomesa=f.read().strip()
genomesa=genomesa[:-33] 

print("frequent 9-mers in Covid 2 genome: ", BetterFrequentWords(genome,9))
print("\n\n")

print("frequent 9-mers in Bat genome: ", BetterFrequentWords(genomebat,9))
print("\n\n")

print("frequent 9-mers in Middle East genome: ", BetterFrequentWords(genomeme,9))
print("\n\n")
 
print("frequent 9-mers in Panglin genome: ", BetterFrequentWords(genomep,9))
print("\n\n")
 
print("frequent 9-mers in SARS genome: ", BetterFrequentWords(genomes,9))
print("\n\n")

print("frequent 9-mers in Sevre acute genome: ", BetterFrequentWords(genomesa,9))


#Task 4: The evolution of the most frequent 9-mer repeat in coronaviruses. 
# Define the 9-mers for each genome
bat = SeqRecord(Seq('CAGCTGGTA'), id='Bat')
me = SeqRecord(Seq('TTAACGAAc'), id='Middle East')
pangolin = SeqRecord(Seq('TAATGGTAA'), id='Pangolin')
sars = SeqRecord(Seq('TAAACGAAc'), id='SARS coronavirus')
severe = SeqRecord(Seq('GATGGTTTT'), id='Severe acute')

# Construct a MultipleSeqAlignment object from the sequences
alignment = MultipleSeqAlignment([bat, me, pangolin, sars, severe])

# Compute the distances between the genomes based on the alignment
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)

# Construct the phylogenetic tree
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# Draw and show the tree using Phylo module
Phylo.draw(tree)




end_time = time.time()
task4_time = end_time - start_time
print(f"----------------Task 4 runtime: {task4_time:.2f} seconds----------------")


#----------------------------------------
#Task E1: 
start_time = time.time()
print(f"----------------Task E1 start ----------------")

with open ('data/HIV1.fasta') as f:
 genomehiv=f.read().strip()

with open ('data/Hepatitis_B.fasta') as f:
 genomehep=f.read().strip()

with open ('data/adenovirus.fasta') as f:
 genomead=f.read().strip()


with open('data/Ebola.fasta') as f:
    genomeb = f.read().strip()

# Print the most frequent 9-mers in each genome sequence
print("Most frequent 9-mers in HIV genome: ", BetterFrequentWords(genomehiv,9))
print("Most frequent 9-mers in Hepatitis B genome: ", BetterFrequentWords(genomehep,9))
print("Most frequent 9-mers in Adenovirus genome: ", BetterFrequentWords(genomead,9))
print("Most frequent 9-mers in Ebola genome: ", BetterFrequentWords(genomeb,9))

#Task E1:
# Define the 9-mers for each genome
hiv = SeqRecord(Seq('AAA GAAAAA'), id='HIV')
hep_b = SeqRecord(Seq('TTC TTGTTG'), id='Hepatitis B')
adeno = SeqRecord(Seq('GGC GGCGGC'), id='Adenovirus')
ebola = SeqRecord(Seq('GAA GATTAA'), id='Ebola')

# Construct a MultipleSeqAlignment object from the sequences
alignment = MultipleSeqAlignment([hiv, hep_b, adeno, ebola])

# Compute the distances between the genomes based on the alignment
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)

# Construct the phylogenetic tree
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# Draw and show the tree using Phylo module
Phylo.draw(tree)

end_time = time.time()
taskE1_time = end_time - start_time
print(f"----------------Task E1 runtime: {taskE1_time:.2f} seconds----------------")