# Load necessary packages
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo import read, draw_ascii
from Bio.Phylo.PAML import baseml
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

# Replace "sequences.fasta" with the path to your FASTA file
records = AlignIO.read("input_file.fasta", "fasta")

# Create a MultipleSeqAlignment object
alignment = MultipleSeqAlignment(records)


# Calculate pairwise distances using the identity metric
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)

# Use NJ to construct the tree
nj_constructor = DistanceTreeConstructor(calculator, 'nj')
nj_tree = nj_constructor.build_tree(alignment)

# Print the tree in ASCII format
print("Neighbor joining tree:")
draw_ascii(nj_tree)


# Infer the ancestral sequences
if baseml is not None:
    ancestral_seqs = {}
    for node in nj_tree.get_nonterminals():
        if node.name:
            # Get the branch length from the tree
            branch_length = nj_tree.distance(node)
            
            # Get the inferred ancestral sequence using baseml
            ancestral_seq = baseml.get_local_alignment(node.clades[0].name, node.clades[1].name, branch_length)[0]
            
            # Store the inferred ancestral sequence
            ancestral_seqs[node] = ancestral_seq
else:
    print("Error: baseml is None")

# Print the inferred ancestral sequences
print("\nInferred ancestral sequences:")
for node in nj_tree.get_nonterminals():
    if node.name:
        print("Node {}: {}".format(node.name, ancestral_seqs[node]))

# Write the inferred ancestral sequences to a file
with open("ancestral_sequences.txt", "w") as handle:
    for node in nj_tree.get_nonterminals():
        if node.name:
            handle.write("Node {}: {}\n".format(node.name, ancestral_seqs[node]))
