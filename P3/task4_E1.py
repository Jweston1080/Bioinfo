
#-----------------------
import pandas as pd
import numpy as np
import timeit
import math
import random
import os
import pathlib
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

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

#----------------------------------------
#Task E1: 
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


