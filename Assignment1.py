# Task: Frequent words by counting. Implement the FrequentWords() algorithm in class based on the PatternCount() method.
    
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count
   
def FrequentWords(Text, k):
    frequent_patterns = []
    counts = []
    for i in range(len(Text) - k + 1):
        pattern = Text[i:i+k]
        count = PatternCount(Text, pattern)
        counts.append(count)
    max_count = max(counts)
    for i in range(len(Text) - k + 1):
        if counts[i] == max_count and Text[i:i+k] not in frequent_patterns:
            frequent_patterns.append(Text[i:i+k])
    return frequent_patterns
   
   
   # Task 2: Frequent words by hashing. Implement the BetterFrequentWords() algorithm in class based on the FrequencyMap() method.
   
def FrequencyMap(Text, k):
    freq_map = {}
    n = len(Text)
    for i in range(n - k + 1):
        pattern = Text[i:i+k]
        if pattern in freq_map:
            freq_map[pattern] += 1
        else:
            freq_map[pattern] = 1
    return freq_map
   
def BetterFrequentWords(Text, k):
    freq_map = FrequencyMap(Text, k)
    max_count = max(freq_map.values())
    frequent_patterns = [pattern for pattern, count in freq_map.items() if count == max_count]
    return frequent_patterns
   
   
   # Task 3: Write a program to generate random DNA sequences of various lengths
   
import numpy as np
   
def generate_random_dna(length):
    nucleotides = ['A', 'C', 'G', 'T']
    probabilities = [0.25, 0.25, 0.25, 0.25]
    sequence = np.random.choice(nucleotides, size=length, p=probabilities)
    return ''.join(sequence)
   