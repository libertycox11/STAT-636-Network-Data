

setwd("C:/GitHub/STAT-636-Multivariate-Network-Data")


library(igraph)

# Read the edge list, skipping the 4 comment lines
edges <- read.table("Amazon0302.txt", 
                    skip = 4, 
                    col.names = c("from", "to"))



# Create an UNDIRECTED graph and simplify to remove duplicates/self-loops
g <- graph_from_data_frame(edges, directed = FALSE)
g <- simplify(g)  # removes duplicate edges and self-loops

# Quick summary
vcount(g)  # number of nodes
ecount(g)  # number of edges