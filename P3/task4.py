import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Define the list of genome sequences
try:
    record = SeqIO.read('/content/COVID2.fasta', 'fasta')
    genome = record.seq
except ValueError:
    print('No records found in file.')

genomebat = SeqIO.read('data/Bat.fasta', 'fasta').seq
genomeme = SeqIO.read('data/Middle_east.fasta', 'fasta').seq
genomep = SeqIO.read('data/Pangolin.fasta', 'fasta').seq
genomes = SeqIO.read('data/SARS.fasta', 'fasta').seq
genomesa = SeqIO.read('data/Severe_acute.fasta', 'fasta').seq

# Define the list of sequences
sequences = [genome, genomebat, genomeme, genomep, genomes, genomesa]

# Calculate pairwise distances
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(sequences)

# Construct the tree
constructor = DistanceTreeConstructor()
upgma_tree = constructor.upgma(dm)

# Plot the tree
fig, axes = plt.subplots(1, 1, figsize=(10, 10), dpi=100)
Phylo.draw(upgma_tree, axes=axes)
plt.show()
