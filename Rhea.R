# Rhea
# set your working directory

#testing this again hopefully


#Exploratory Characterization
library(igraph); library(Matrix)

# Load and clean 
edges <- read.table("Amazon0302.txt", comment.char = "#")
G <- graph_from_edgelist(as.matrix(edges) + 1, directed = FALSE)  # SNAP IDs start at 0
G <- simplify(G)         # removes reciprocal duplicates and self-loops
vcount(G); ecount(G)     # 262,111 nodes, 899,792 edges
mean(degree(G))          # 6.87

# Connectivity (2-1) 
components(G)            # 1 component, so the graph is connected

# Summary statistics (2-1) 
#transitivity = clustering coefficient
transitivity(G)          # 0.236 (proportion of “two-step” connections that are also direct connections)
edge_density(G)          # 2.6e-05

# Comparison to ER graph with the same density (3-3)
G_np = sample_gnp(vcount(G), p = edge_density(G))
transitivity(G_np)       # near 0, far below the real network

# Degree distribution (1-1, 2-1) 
degree_vals = 0:max(degree(G))
plot(degree_vals, degree_distribution(G), log = "xy",
     xlab = "Degree", ylab = "Proportion", main = "Degree distribution")

# Explains the spike at degree 5: crawl recorded at most 5 co-purchases per product
table(table(edges[, 1]))

# Hill estimator and Hill plot (3-1) 
deg = degree(G)

t_min = 10
Y = deg[which(deg >= t_min)]
alpha = 1 + 1/mean(log(Y/(min(Y)-0.5)))
alpha

sort.deg = sort(deg, decreasing = T)
alpha_t = sapply(1:5000, function(k) {
  Y = sort.deg[1:k]
  alpha = 1 + 1/mean(log(Y/sort.deg[k+1]))
  alpha
})
plot(alpha_t)

# Centrality (2-1)
sort(degree(G), decreasing = T)[1:10]
eigencent = eigen_centrality(G) #not necessarily needed but in notes
sort(eigencent$vector, decreasing = T)[1:10] #same here