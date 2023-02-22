import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import timeit
import math
import random


# Task 1: Frequent words by counting. Implement the FrequentWords() algorithm in class based on the PatternCount() method.
    
def PatternCount(Word, Pattern):
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
    n = len(Word)
    for i in range(n - k + 1):
        pattern = Word[i:i+k]
        if pattern in freq_map:
            freq_map[pattern] += 1
        else:
            freq_map[pattern] = 1
    return freq_map
   
def BetterFrequentWords(Word, k):
    freq_map = FrequencyMap(Word, k)
    max_count = max(freq_map.values())
    frequent_patterns = [pattern for pattern, count in freq_map.items() if count == max_count]
    return frequent_patterns
   
   
# Task 3: Write a program to generate random DNA sequences of various lengths
   
def generate_random_dna(length):
    n = ['A', 'C', 'G', 'T']
    p = [0.25, 0.25, 0.25, 0.25]
    sequence = np.random.choice(n, size=length, p=p)
    return ''.join(sequence)

# Task 4:Visualize the empirical runtime of the two algorithms by counting versus by hashing
ns = np.linspace(100, 2000, 15, dtype=int)
ts = [timeit.timeit('FrequentWords(lst, 0)',
                    setup='lst=list(range({})); random.shuffle(lst)'.format(n),
                    globals=globals(),
                    number=1)
         for n in ns]
plt.plot(ns, ts, 'or');

degree = 4
coeffs = np.polyfit(ns, ts, degree)
p = np.poly1d(coeffs)
plt.plot(ns, [p(n) for n in ns], '-r')

# implement plot for BetterFrequentWords

#Task 5: 
with open ('/content/sample_data/sequence.txt') as f:
  genome=f.read().strip()
  genome=genome[:-33] #Frequent k-mers in the SARS-Cov-2 genome

print("frequent 3-mers in Covid genome: ", BetterFrequentWords(genome,3)) 
print("frequent 6-mers in Covid genome: ", BetterFrequentWords(genome,6))
print("frequent 9-mers in Covid genome: ", BetterFrequentWords(genome,9))
print("frequent 12-mers in Covid genome: ", BetterFrequentWords(genome,12))
print("frequent 15-mers in Covid genome: ", BetterFrequentWords(genome,15))

# Testing of functions

def test():
  r = generate_random_dna(10)
  print('Testing FrequentWords(Task1) with string: ' + r)
  print(FrequentWords(r, 3))

  f = generate_random_dna(20);
  print('Testing BetterFrequentWords(Task2) with string: ' + f)
  print(BetterFrequentWords(f, 3))

def main():
  test()

if __name__ == "__main__":
  main()

   