# function to calculate Hamming distance between two sequences
def hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2)
    return sum(s1 != s2 for s1, s2 in zip(seq1, seq2))

# read in fastq file
with open("input_file.fastq", "r") as f:
    seq_records = f.readlines()

# create list of sequences from fastq file
sequences = [seq_records[i+1].strip() for i in range(0, len(seq_records), 4)]

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
with open("output_file.csv", "w") as f:
    for row in distance_matrix:
        f.write(",".join(str(x) for x in row) + "\n")