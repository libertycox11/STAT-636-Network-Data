# ===================================================
# TITLE: STAT 636 Multivariate & Network Data Analysis
# DATE: Fall 2026
# AUTHORS: Liberty Cox, Aysegul Buber, Rhea Perumalil
# ===================================================

setwd("C:/GitHub/STAT-636-Network-Data")

# ===================================================
# DATA CLEANING
# ===================================================
# this is a random comment

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





# ============================================================
# ANALYSIS METHOD 2:
# Community Detection using Stochastic Block Modeling
# ============================================================
# Class code sources used:
#   2.2 Graph Partitioning  --> Laplacian, eigen, kmeans clustering
#   3.4 Estimation for SBMs --> eigs() for large sparse matrices,
#                               kmeans(V, K, nstart=100), table(cluster,z)
#   3.5 DCSBM & MMSBM       --> levelplot on reordered adjacency matrix,
#                               index_order = order() reordering pattern
# ============================================================

library(igraph)
library(Matrix)
library(lattice)
library(rARPACK)   # efficient top-k eigenvectors for large sparse matrices

# NOTE: This section assumes g and G_lcc were already built in Part 1.
#       All community detection is performed on G_lcc (the largest
#       connected component), which is standard practice -- the LCC
#       ensures the graph is fully connected, which is required for
#       the normalized Laplacian to be well-defined.

# ============================================================
# STEP 1: Build the Adjacency and Laplacian Matrices
# ============================================================
# Pattern: get.adjacency(), D = diag(degrees), L = D - A  -- (2.2)

A_lcc <- get.adjacency(G_lcc)          # get.adjacency()     -- 2.2 / 1.2
n_lcc <- vcount(G_lcc)                 # number of nodes in LCC

# Degree vector and diagonal degree matrix
degrees_lcc <- degree(G_lcc)           # degree()            -- 1.1 / 2.2
D_lcc       <- Diagonal(x = degrees_lcc)  # sparse diagonal  -- efficient for large n

# Unnormalized Laplacian: L = D - A
L_lcc <- D_lcc - A_lcc                 # L = D - A           -- 2.2

# Normalized Laplacian: L_norm = D^{-1/2} L D^{-1/2}
# Pattern: D_m12 = diag(colSums(A)^(-1/2)), L_norm = D_m12 %*% L %*% D_m12 -- (2.2)
# Using sparse version for efficiency on 260k-node graph
D_m12_lcc  <- Diagonal(x = degrees_lcc^(-1/2))   # D^{-1/2}  -- 2.2
L_norm_lcc <- D_m12_lcc %*% L_lcc %*% D_m12_lcc  # L_norm    -- 2.2

cat("Adjacency and Laplacian matrices built.\n")

# ============================================================
# STEP 2: Choose Number of Communities K via Scree Plot
# ============================================================
# Pattern: plot(eig.L$values, main = "Scree plot of L")    -- (2.2)
# For large sparse matrices, use eigs() NOT eigen()
# Pattern: eigs(A, k) for top-k eigenvectors               -- (3.4)
# For the Laplacian we want the SMALLEST eigenvalues,
# so we use eigs_sym() with which = "SM" (smallest magnitude)

K_max <- 10   # inspect top K_max smallest eigenvalues to choose K

eig_L_small <- eigs_sym(
  L_norm_lcc,
  k     = K_max,
  which = "SM"       # Smallest Magnitude -- for Laplacian community detection
)

# Eigenvalues come back in descending order from eigs_sym;
# reverse so they read smallest -> largest for the scree plot
eig_vals_sorted <- rev(eig_L_small$values)
eig_vecs_sorted <- eig_L_small$vectors[, K_max:1]

# Scree plot of smallest eigenvalues -- plot(eig.L$values) pattern from 2.2
plot(
  1:K_max, eig_vals_sorted,
  type = "b",
  pch  = 19,
  col  = "steelblue",
  main = "Scree Plot: Smallest Eigenvalues of Normalized Laplacian (LCC)",
  xlab = "Index (smallest to largest)",
  ylab = "Eigenvalue"
)
abline(v = 4, col = "red", lty = 2)   # adjust after inspecting the plot
text(4, max(eig_vals_sorted) * 0.5,
     "K = 4", col = "red", adj = c(-0.2, 0))

# ============================================================
# STEP 3: Spectral Clustering with Chosen K
# ============================================================
# Pattern: V = eigs(A, 2)$vectors
#          clusters = kmeans(V, 2, nstart = 100)$cluster    -- (3.4)

K <- 4   # set based on scree plot above;
# Amazon metadata has 4 categories: Books, Music, DVD, Video

# Top K eigenvectors of normalized Laplacian (smallest eigenvalues)
eig_K <- eigs_sym(
  L_norm_lcc,
  k     = K,
  which = "SM"
)

# Reorder columns smallest -> largest eigenvalue
V_K <- eig_K$vectors[, K:1]           # eigenvector matrix, n x K

# Spectral clustering: kmeans on the K-dim eigenvector embedding
# Pattern: kmeans(V, K, nstart = 100)$cluster              -- (3.4)
set.seed(636)
spectral_clusters <- kmeans(
  V_K,
  centers = K,
  nstart  = 100,     # nstart = 100 for stability                -- 3.4
  iter.max = 500
)$cluster

cat("\n--- Spectral Clustering (K =", K, ") ---\n")
cat("Cluster sizes:\n")
print(sort(table(spectral_clusters), decreasing = TRUE))

# ============================================================
# STEP 4: Visualize Eigenvector Embedding (2D)
# ============================================================
# Pattern: plot(V, col = z)                                 -- (3.4 / 3.5)

# Plot first two eigenvectors colored by cluster assignment
plot(
  V_K[, 1], V_K[, 2],
  col  = spectral_clusters,
  pch  = 19,
  cex  = 0.4,
  main = "Spectral Embedding: Top 2 Eigenvectors of L_norm\nColored by Cluster",
  xlab = "Eigenvector 1 (Fiedler)",
  ylab = "Eigenvector 2"
)
legend(
  "topright",
  legend = paste("Cluster", 1:K),
  col    = 1:K,
  pch    = 19,
  cex    = 0.8
)

# ============================================================
# STEP 5: Reorder Adjacency Matrix by Cluster (Block Structure)
# ============================================================
# Pattern: index_order = order(Z[,2])
#          Anew = A[index_order, index_order]
#          levelplot(as.matrix(Anew))                       -- (3.5)
# NOTE: plotting the full 260k x 260k matrix is not feasible;
#       we use a subgraph of the top 500 highest-degree nodes
#       split across clusters, same subgraph approach as Part 1

top500_ids   <- order(degree(G_lcc), decreasing = TRUE)[1:500]
A_sub500     <- A_lcc[top500_ids, top500_ids]      # submatrix   -- 1.2
cluster_sub  <- spectral_clusters[top500_ids]      # their labels

# Reorder rows/cols by cluster -- index_order pattern from 3.5
index_order  <- order(cluster_sub)
A_reordered  <- A_sub500[index_order, index_order]

# Levelplot of reordered adjacency (block structure should appear)
levelplot(
  as.matrix(A_reordered),
  main   = "Reordered Adjacency Matrix: Top 500 Nodes by Cluster",
  xlab   = "Node (sorted by cluster)",
  ylab   = "Node (sorted by cluster)",
  col.regions = colorRampPalette(c("white", "steelblue"))(100)
)

# ============================================================
# STEP 6: Cluster Summary Statistics
# ============================================================
# Pattern: sort(degree(), decreasing=T)[1:10]
#          table(cluster, z)                                -- (2.2 / 2.1)

# Assign cluster labels back to LCC vertex attributes
V(G_lcc)$cluster <- spectral_clusters          # attach to graph  -- 1.1 pattern

# Degree statistics within each cluster
cat("\n--- Per-Cluster Degree Summary ---\n")
deg_lcc <- degree(G_lcc)

cluster_summary <- do.call(rbind, lapply(1:K, function(k) {
  idx <- which(spectral_clusters == k)
  data.frame(
    Cluster    = k,
    Size       = length(idx),
    Pct_Network = round(100 * length(idx) / n_lcc, 2),
    Mean_Degree = round(mean(deg_lcc[idx]), 4),
    Median_Degree = median(deg_lcc[idx]),
    Max_Degree  = max(deg_lcc[idx])
  )
}))

print(cluster_summary, row.names = FALSE)

# Top 10 nodes by degree within each cluster
# Pattern: sort(degree(G_lcc), decreasing = T)[1:10]       -- (2.1)
cat("\n--- Top 10 Nodes by Degree Per Cluster ---\n")
for (k in 1:K) {
  idx      <- which(spectral_clusters == k)
  top_k    <- sort(deg_lcc[idx], decreasing = TRUE)[1:10]
  cat(sprintf("\nCluster %d:\n", k))
  print(top_k)
}

# ============================================================
# STEP 7: Estimate SBM Block Probability Matrix B
# ============================================================
# Pattern: B_matrix = crossprod(tau, crossprod(A, tau)) /
#                     tcrossprod(apply(tau, 2, sum))        -- (3.4)
# Here we use hard cluster assignments (tau = indicator matrix)
# instead of soft variational assignments, following the
# spectral clustering result.

# Build hard assignment matrix tau (n x K indicator)       -- 3.4 pattern
tau <- matrix(0, nrow = n_lcc, ncol = K)
tau[cbind(1:n_lcc, spectral_clusters)] <- 1

# Estimate B: within/between cluster edge probabilities
cluster_sizes <- apply(tau, 2, sum)                        # cluster sizes

# Numerator: observed edges between clusters
B_num   <- crossprod(tau, crossprod(A_lcc, tau))           # 3.4 pattern

# Denominator: possible edges between clusters
B_denom <- tcrossprod(cluster_sizes)                       # 3.4 pattern
# Correct diagonal: within-cluster has n_k*(n_k-1)/2 possible edges (no self-loops)
diag(B_denom) <- cluster_sizes * (cluster_sizes - 1)

# Estimated B matrix
B_hat <- as.matrix(B_num / B_denom)
rownames(B_hat) <- paste0("C", 1:K)
colnames(B_hat) <- paste0("C", 1:K)

cat("\n--- Estimated SBM Block Probability Matrix (B_hat) ---\n")
print(round(B_hat, 6))

# Heatmap of B_hat -- levelplot(P) pattern from 3.4 / 3.5
levelplot(
  B_hat,
  main      = "Estimated SBM Block Probability Matrix",
  xlab      = "Community",
  ylab      = "Community",
  col.regions = colorRampPalette(c("white", "steelblue"))(100)
)

# ============================================================
# STEP 8: Modularity of the Spectral Partition
# ============================================================
# Modularity quantifies how well the partition separates
# within-community from between-community edges.
# A value close to 1 = strong community structure.

modularity_val <- modularity(
  G_lcc,
  membership = spectral_clusters
)
cat("\n--- Modularity of Spectral Partition ---\n")
cat("Modularity (K =", K, "):", round(modularity_val, 6), "\n")
cat("(Range: -1 to 1; values > 0.3 indicate meaningful structure)\n")

# ============================================================
# STEP 9: Final Summary
# ============================================================

cat("\n========================================\n")
cat("     SBM COMMUNITY DETECTION SUMMARY\n")
cat("========================================\n")
cat(sprintf("  Method              : Spectral Clustering\n"))
cat(sprintf("  Matrix used         : Normalized Laplacian\n"))
cat(sprintf("  Num communities (K) : %d\n",   K))
cat(sprintf("  LCC nodes           : %s\n",   format(n_lcc, big.mark=",")))
cat(sprintf("  Modularity          : %.6f\n", modularity_val))
cat("  Cluster sizes:\n")
for (k in 1:K) {
  cat(sprintf("    Cluster %d : %s nodes (%.2f%%)\n",
              k,
              format(cluster_summary$Size[k], big.mark=","),
              cluster_summary$Pct_Network[k]))
}
cat("========================================\n")



# ===================================================
# ANALYSIS METHOD 3:
# Community Detection using Louvain Method
# ===================================================

