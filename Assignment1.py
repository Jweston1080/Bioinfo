# Assignment: P1
# Authors: Jason Weston, Raghda Kailany, Calicia Perea
#%%
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt2
import pandas as pd
import numpy as np
import timeit
import math
import random
import os
import pathlib

# Task 1: Frequent words by counting. Implement the FrequentWords() algorithm in class based on the PatternCount() method.
def PatternCount(Word, Pattern):
  """"""
  count = 0
  for i in range(len(Word) - len(Pattern) + 1):
      if Word[i:i+len(Pattern)] == Pattern:
          count += 1
  return count
   
def FrequentWords(Word, k):
  frequent_patterns = []
  counts = []
  for i in range(len(Word) - k + 1):
      pattern = Word[i:i+k]
      count = PatternCount(Word, pattern)
      counts.append(count)
  max_count = max(counts)
  for i in range(len(Word) - k + 1):
      if counts[i] == max_count and Word[i:i+k] not in frequent_patterns:
          frequent_patterns.append(Word[i:i+k])
  return frequent_patterns
   
   
   # Task 2: Frequent words by hashing. Implement the BetterFrequentWords() algorithm in class based on the FrequencyMap() method.
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
  return frequent_patterns
   
   
# Task 3: Write a program to generate random DNA sequences of various lengths
def generate_random_dna(length):
  n = ['A', 'C', 'G', 'T']
  p = [0.25, 0.25, 0.25, 0.25]
  sequence = np.random.choice(n, size=length, p=p)
  return ''.join(sequence)

# Task 4:Visualize the empirical runtime of the two algorithms by counting versus by hashing
# Graph for FrequentWords()
ns = np.linspace(100, 2000, 15, dtype=int)
ts = [timeit.timeit('FrequentWords(lst, 0)',
                    setup='lst=list(range({})); random.shuffle(lst)'.format(n),
                    globals=globals(),
                    number=1)
         for n in ns]
plt.plot(ns, ts, 'or');

# Graph for BetterFrequentWords()
bf = np.linspace(100, 2000, 15, dtype=int)
bfr = [timeit.timeit('BetterFrequentWords(lst2, 0)',
                    setup='lst2=list(range({})); random.shuffle(lst2)'.format(n),
                    globals=globals(),
                    number=1)
         for n in ns]
plt2.plot(bf, bfr, 'ob');

#Task 5: Task 5: Report the most frequent k-mers for each k in 3, 6, 9, 12, 15. 
with open ('sequence.txt') as f:
  genome=f.read().strip()
  genome=genome[:-33] # Frequent k-mers in the SARS-Cov-2 genome
print("frequent 3-mers in Covid genome: ", BetterFrequentWords(genome,3)) 
print("\n\n")
print("frequent 6-mers in Covid genome: ", BetterFrequentWords(genome,6))
print("\n\n")
print("frequent 9-mers in Covid genome: ", BetterFrequentWords(genome,9))
print("\n\n")
print("frequent 12-mers in Covid genome: ", BetterFrequentWords(genome,12))
print("\n\n")
print("frequent 15-mers in Covid genome: ", BetterFrequentWords(genome,15))
print("\n\n")

# Testing of functions
def test():
  r = generate_random_dna(20)
  print('Testing FrequentWords(Task1) with string: ' + r)
  print(FrequentWords(r, 3))
  print("\n\n")
  f = generate_random_dna(20);
  print('Testing BetterFrequentWords(Task2) with string: ' + f)
  print(BetterFrequentWords(f, 3))
  print("\n\n")

  print("This graph displays the runtime of the FrequentWords() method in red")
  print("and the BetterFrequentWords() method in blue")
def main():
  test()

if __name__ == "__main__":
  main()

# %%
