
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
print(f"Task 1 runtime: {task1_time:.2f} seconds")
#------------Task 2 ------------------#

start_time = time.time()

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
print(f"Task 2 runtime: {task2_time:.2f} seconds")

#------------Task 3 ------------------#
#start_time = time.time()




# Read the input file
'''
end_time = time.time()
task3_time = end_time - start_time
print(f"Task 3 runtime: {task3_time:.2f} seconds")'''


#------------Task 4 ------------------#
#start_time = time.time()

# genomic sequence files of multiple coronaviruses
# ... rest of Task 4 code ...

#end_time = time.time()
#task4_time = end_time - start_time
#print(f"Task 4 runtime: {task3_time:.2f} seconds")


