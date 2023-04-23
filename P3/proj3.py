from ete3 import Tree,TreeStyle,TextFace, PhyloTree

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import average, to_tree
import numpy as np
import time 

#------------Task 1 ------------------#
start_time = time.time()

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

end_time = time.time()
task1_time = end_time - start_time
print(f"Task 1 runtime: {task1_time:.2f} seconds")
#------------Task 2 ------------------#

start_time = time.time()
# read in distance matrix from output file of Task 1
with open("output_file.csv", "r") as f:
    lines = f.readlines()
    matrix = [list(map(int, line.strip().split(","))) for line in lines]

# convert distance matrix to condensed distance matrix format
condensed_matrix = squareform(matrix)

# compute average linkage hierarchical clustering
linkage_matrix = average(condensed_matrix)

# convert linkage matrix to tree structure
tree = to_tree(linkage_matrix, False)

# helper function to convert tree structure to Newick format
def tree_to_newick(tree):
    if tree.is_leaf():
        return f"{tree.id}"
    else:
        left_subtree = tree_to_newick(tree.left)
        right_subtree = tree_to_newick(tree.right)
        branch_length = "{:.2f}".format(tree.dist)
        return f"({left_subtree}:{branch_length},{right_subtree}:{branch_length}){tree.id}"

# convert tree structure to Newick format
newick_str = tree_to_newick(tree)

# write tree to output file in Newick format
with open("output_file.nwk", "w") as f:
    f.write(newick_str)



# read in tree from output file of Task 2
with open("output_file.nwk", "r") as f:
    newick_str = f.read().strip() + ";"  # add semicolon to the end

# try to create tree object from Newick string
try:
    tree = Tree(newick_str)
except:
    print("Error: unable to create tree object from Newick string.")
    exit()

# check if the tree has at least one child
if len(tree.children) > 0:
    # set the first child as the root and add a branch length of 0
    tree.set_outgroup(tree.children[0])
    tree.children[0].dist = 0.0
else:
    print("Error: Tree has no children.")

# visualize tree
ts = TreeStyle()
ts.show_leaf_name = True
ts.scale = 10000 
ts.branch_vertical_margin = 10
ts.layout_fn = "phylogeny"
ts.title.add_face(TextFace("Hierarchical Clustering Tree", fsize=20), column=0)
for n in tree.traverse():
    if not n.is_leaf():
        nstyle = TextFace("{:.2f}".format(n.dist), fsize=10, fgcolor="black")
        n.add_face(nstyle, column=0)

# attempts to render the tree to a PDF file and if there is an error it prints a message and exits.
try:
    tree.render("tree_task2.pdf", tree_style=ts)
except:
    print("Error: unable to render tree to PDF file.")
    exit()


# save tree object for use in Task 3
tree.write(format=1, outfile="tree_task2.nw")

end_time = time.time()
task2_time = end_time - start_time
print(f"Task 2 runtime: {task2_time:.2f} seconds")

#------------Task 3 ------------------#
#start_time = time.time()

# read in tree from output file of Task 2
# ... rest of Task 3 code ...

#end_time = time.time()
#task3_time = end_time - start_time
#print(f"Task 3 runtime: {task3_time:.2f} seconds")


#------------Task 4 ------------------#
#start_time = time.time()

# genomic sequence files of multiple coronaviruses
# ... rest of Task 4 code ...

#end_time = time.time()
#task4_time = end_time - start_time
#print(f"Task 4 runtime: {task3_time:.2f} seconds")


