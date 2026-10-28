# ===================================================
# TITLE: STAT 636 Multivariate & Network Data Analysis
# DATE: Fall 2026
# AUTHORS: Liberty Cox, Aysegul Buber, Rhea Perumalil
# ===================================================

setwd("C:/GitHub/STAT-636-Multivariate-Network-Data")


library(igraph)

# Read the edge list, skipping the 4 comment lines
edges <- read.table("Amazon0302.txt", skip = 4, col.names = c("from", "to"))


# Create an UNDIRECTED graph and simplify to remove duplicates/self-loops
g <- graph_from_data_frame(edges, directed = FALSE)
g <- simplify(g)  # removes duplicate edges and self-loops

# Quick summary
vcount(g)  # number of nodes
ecount(g)  # number of edges

# ===================================================
# METHOD 1:
# Exploratory Network Characterization
# ===================================================

# mean node degree (the average number of connections a product has)


# network density (the proportion of possible connections that are present)


# check the network’s connectivity by identifying 
# product connections
# product communities
# the size of each community


# examine the degree distribution for right skewness 
# to determine whether a small number of products hold a 
# disproportionate number of connections

# products with the most connections ranked by degree centrality


# ===================================================
# METHOD 2:
# Community Detection using Stochastic Block Modeling
# ===================================================

# reference class code notes 3-4 SBM.R

# ===================================================
# METHOD 3:
# Community Detection using Louvain Method
# ===================================================

