library(phangorn)
library(ape)
library(seqinr)

# Task 2: Building a rooted phylogenetic tree

# load the distance matrix from Task 1
dist_matrix <- as.dist(read.csv("output_file.csv", header=FALSE)$V1)

# compute the NJ tree using phangorn package
tree <- NJ(dist_matrix)


# set the first child as the root and add a branch length of 0
tree$edge.length[1] <- 0
tree$root.edge <- tree$edge.length[1]

# Visualize
# create a pdf file to save the plot
pdf("tree_nj.pdf")
plot(tree, cex = 0.7, main = "Neighbor Joining Tree")


# save the tree to a file in Newick format
write.tree(tree, file="tree_nj.nwk")

# close the pdf file
dev.off()