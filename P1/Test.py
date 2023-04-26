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

with open ('sequence.txt') as f:
  genome=f.read().strip()
  genome=genome[:-33] # Frequent k-mers in the SARS-Cov-2 genome





print("frequent 9-mers in Covid genome: ", BetterFrequentWords(genome,9))
