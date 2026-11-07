# ===================================================
# TITLE: STAT 636 Multivariate & Network Data Analysis
# DATE: Fall 2026
# AUTHORS: Liberty Cox, Aysegul Buber, Rhea Perumalil
# ===================================================

setwd("C:/GitHub/STAT-636-Network-Data")

# ===================================================
# DATA CLEANING
# ===================================================

library(igraph)
library(Matrix)
library(lattice)

# Read the edge list, skipping the 4 comment lines
edges_raw <- read.table("Amazon0302.txt", skip = 4, col.names = c("from", "to"))

# Build directed first (as given in the raw file), then convert
g_directed <- graph_from_edgelist(as.matrix(edges_raw), directed = TRUE)

# Convert to undirected, collapse duplicate edges, remove self-loops
g <- as.undirected(g_directed, mode = "collapse")
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

cat("Network loaded and cleaned.\n")
cat("Nodes:", vcount(g), "\n")
cat("Edges:", ecount(g), "\n")

# ===================================================
#  ANALYSIS METHOD 1:
# Exploratory Network Characterization
# ===================================================

# Adjacency Matrix (spot-check; full matrix is too large to plot)
top50_ids  <- order(degree(g), decreasing = TRUE)[1:50]
g_sub50    <- induced_subgraph(g, top50_ids)
A_sub      <- get.adjacency(g_sub50)
levelplot(as.matrix(A_sub),
          main = "Adjacency Matrix: Top 50 Highest-Degree Products")

# ------------------------------------------------------------
# A. mean node degree (avg. # of connections a product has)
#    & network density (proportion of possible connections that are present)
# ------------------------------------------------------------

num_nodes   <- vcount(g)
num_edges   <- ecount(g)
net_density <- edge_density(g, loops = FALSE)   # edge_density() -- 2-1
mean_deg    <- mean(degree(g))                  # degree()       -- 1-1 / 2-1
med_deg     <- median(degree(g))

cat("\n--- Global Network Properties ---\n")
cat("Nodes             :", format(num_nodes, big.mark = ","), "\n")
cat("Edges             :", format(num_edges, big.mark = ","), "\n")
cat("Network density   :", round(net_density, 8), "\n")
cat("Mean degree       :", round(mean_deg, 4), "\n")
cat("Median degree     :", med_deg, "\n")


# ------------------------------------------------------------
# B. check the network’s connectivity by identifying 
#     i. product connections
#     ii. product communities
#     iii. the size of each community
# ------------------------------------------------------------

G_comp    <- components(g)
num_comps <- G_comp$no
comp_sizes <- sort(G_comp$csize, decreasing = TRUE)

cat("\n--- Connected Components ---\n")
cat("Number of components          :", num_comps, "\n")
cat("Top 10 component sizes        :", head(comp_sizes, 10), "\n")
cat("Largest component (nodes)     :", comp_sizes[1],
    sprintf("(%.2f%% of network)\n", 100 * comp_sizes[1] / num_nodes))

# Component size distribution table
cat("\nComponent size distribution (size : count):\n")
print(table(G_comp$csize))

# Extract the Largest Connected Component (LCC)
# Pattern: membership == which.max(csize), induced_subgraph()
lcc_id    <- which.max(G_comp$csize)
vert_lcc  <- V(g)[G_comp$membership == lcc_id]  # V(G_stat)[membership==1]
G_lcc     <- induced_subgraph(g, vert_lcc)

cat("\nLCC nodes :", vcount(G_lcc), "\n")
cat("LCC edges :", ecount(G_lcc), "\n")

# Summary stats ON the LCC
cat("\n--- LCC Summary Statistics ---\n")
cat("Edge density:", round(edge_density(G_lcc),  8), "\n")
cat("Transitivity:", round(transitivity(G_lcc),  4), "\n")
cat("Diameter:", diameter(G_lcc), "\n")  
cat("Radius:", radius(G_lcc), "\n")

# ------------------------------------------------------------
# C. examine the degree distribution for right skewness 
#     to determine whether a small number of products hold a 
#     disproportionate number of connections
#            from 1-1 / 2-1 class notes
# ------------------------------------------------------------
deg_lcc <- degree(G_lcc) 

# Skewness
deg_skew <- mean((deg_lcc - mean(deg_lcc))^3) / (sd(deg_lcc)^3)
cat("\n--- Degree Distribution (LCC) ---\n")
cat("Min    :", min(deg_lcc),              "\n")
cat("Median :", median(deg_lcc),           "\n")
cat("Mean   :", round(mean(deg_lcc), 4),   "\n")
cat("Max    :", max(deg_lcc),              "\n")
cat("Skewness:", round(deg_skew, 4), " --> Right-skewed if > 0\n")

# Plot 1: Degree distribution bar plot 
f_d   <- degree_distribution(G_lcc) 
max_d <- length(f_d) - 1
barplot(f_d,
        names.arg = 0:max_d,
        main      = "Degree Distribution (LCC)",
        xlab      = "Degree",
        ylab      = "Proportion",
        col       = "steelblue",
        border    = "white",
        las       = 2)

# Plot 2: Log-log degree distribution
plot(0:max_d, f_d,
     log  = "xy",     
     pch  = 19,
     cex  = 0.6,
     col  = "steelblue",
     main = "Degree Distribution: Log-Log Scale (LCC)",
     xlab = "Degree (log scale)",
     ylab = "Proportion (log scale)")

# Plot 3: Standard degree distribution plot 
degree_vals <- 0:max(deg_lcc)
plot(degree_vals,
     degree_distribution(G_lcc),
     xlab = "Degree",
     ylab = "Proportion",
     main = "Degree Distribution (LCC)",
     pch  = 19,
     cex  = 0.5,
     col  = "steelblue")


# ------------------------------------------------------------
# D. products w/ most connections ranked by degree centrality
#       from 2-1 notes
# ------------------------------------------------------------

cat("\n--- Top 20 Products by Degree (LCC) ---\n")
top_deg <- sort(degree(G_lcc), decreasing = TRUE)[1:20]
print(top_deg)

# Eigenvector centrality
eigencent <- eigen_centrality(G_lcc)
cat("\n--- Top 20 Products by Eigenvector Centrality (LCC) ---\n")
top_eig <- sort(eigencent$vector, decreasing = TRUE)[1:20]
print(round(top_eig, 6))

# Combined summary table
top20_ids   <- names(top_deg)
top20_table <- data.frame(
  Rank             = 1:20,
  Node_ID          = top20_ids,
  Degree           = as.integer(top_deg),
  Eigenvec_Centrality = round(eigencent$vector[top20_ids], 8)
)

cat("\n--- Combined Top 20 Centrality Table ---\n")
print(top20_table, row.names = FALSE)

# Bar chart of top 20 by degree
barplot(top_deg,
        las       = 2,
        col       = "steelblue",
        main      = "Top 20 Products by Degree (LCC)",
        xlab      = "Product Node ID",
        ylab      = "Degree",
        cex.names = 0.75)

# ------------------------------------------------------------
# E. Final Console Summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("    EXPLORATORY NETWORK SUMMARY\n")
cat("========================================\n")
cat(sprintf("  Nodes (full graph)   : %s\n",   format(num_nodes,big.mark=",")))
cat(sprintf("  Edges (full graph)   : %s\n",   format(num_edges,big.mark=",")))
cat(sprintf("  Network Density      : %.2e\n", net_density))
cat(sprintf("  Mean Degree          : %.4f\n", mean_deg))
cat(sprintf("  Median Degree        : %d\n",   as.integer(med_deg)))
cat(sprintf("  # Components         : %d\n",   num_comps))
cat(sprintf("  LCC Size (nodes)     : %s (%.2f%%)\n",
            format(vcount(G_lcc), big.mark=","),
            100 * vcount(G_lcc) / num_nodes))
cat(sprintf("  LCC Density          : %.2e\n", edge_density(G_lcc)))
cat(sprintf("  LCC Transitivity     : %.4f\n", transitivity(G_lcc)))
cat(sprintf("  LCC Diameter         : %d\n",   diameter(G_lcc)))
cat(sprintf("  Degree Skewness      : %.4f\n", deg_skew))
cat("========================================\n")

# ============================================================
# Amazon Co-Purchasing Network Analysis
# Part 1: Exploratory Network Characterization
# ============================================================
# Class code sources used:
#   1-1 Introduction         --> graph building, degree, components, subgraph
#   1-2 Eigenvalues          --> adjacency matrix, levelplot
#   2-1 Descriptive Stats    --> degree_distribution, components, LCC,
#                                edge_density, transitivity, centrality,
#                                top-N ranking pattern
# ============================================================





# ===================================================
# ANALYSIS METHOD 2:
# Community Detection using Stochastic Block Modeling
# ===================================================

# reference class code notes 3-4 SBM.R and 2-2 Graph paritioning.R

# ===================================================
# ANALYSIS METHOD 3:
# Community Detection using Louvain Method
# ===================================================

