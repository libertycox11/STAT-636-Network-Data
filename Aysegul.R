##########################
#Line 44 of the RMD: change to...
##########################
# On average, each product is linked to about 7 other products in the
# co-purchasing network, counting both the products on its own
# "Customers Also Bought" list and the products whose lists include it.
# so not necessarily "its 'Customers Also Bought' list
# contains roughly 7 items"
# since degree counts both directions (undirected)
# products on its own “Also Bought” list (out-links), and
# products whose lists include it (in-links).



##########################
# Line 91 of the RMD...
# Degree distribution
##########################
# degree_distribution() returns a vector of proportions whose FIRST entry is
# degree 0, the second is degree 1, and so on:
#   f_d[1] = P(degree = 0), f_d[2] = P(degree = 1), ..., f_d[k] = P(degree = k - 1)
f_d = degree_distribution(G)

# The true degree for each entry of f_d: 0, 1, 2, ..., max degree.
# Plotting f_d alone uses positions 1, 2, 3, ... as x, which shifts every point
# one unit to the right (degree d drawn at x = d + 1). Using d fixes that.
d = 0:(length(f_d) - 1)

# ---- Consistency checks ----
# 1. The last degree should equal the largest degree in the graph
max(d) == max(degree(G))              # should be TRUE
# 2. The proportions should add up to 1
sum(f_d)                              # should be 1 (up to rounding)
# 3. No isolated products (degree 0) should remain after removing product 0
f_d[1]                                # proportion with degree 0: should be 0
sum(degree(G) == 0)                   # number of isolated products: should be 0
# 4. The spike from the 5-item cap should now sit at degree 5
d[which.max(f_d)]                     # should be 5


# check out the highest degree and understand the degree distribution
max(degree(G))        # 420; the highest degree
sum(degree(G) == 420) # 1; how many products have exactly 420 links
sum(degree(G) > 10)   # number of products with more than 10 links
mean(degree(G) > 10)  # 0.09742856; their share of all products
median(degree(G))     # 6; 50% of the products have degree <= 6! (skewness sign)

barplot(f_d)
#zoom in on the degree distribution barplot
max_show = 50 
barplot(f_d[1:(max_show + 1)],
        names.arg = 0:max_show, 
        xlab = "Degree", ylab = "Proportion",
        main = "Degree distribution (degrees 0-50)")
# Because the degree distribution is right-skewed, this mean is pulled up
# by a few very popular products; a typical product has fewer links
# (median = 6). Each product's own "Also Bought" list has at most 5 items.


##########################
# Line 97 of the RMD...
##########################
# The crawler recorded at most 5 "Customers Also Bought" items per product page,
# so every product has at most 5 OUTGOING links (out-degree <= 5), as this table shows.
# Degree in our undirected graph = outgoing + incoming links, and incoming links
# (other products listing this one) are not capped. So:
#  - degrees above 5 come from products that appear on other products' lists;
#  - the spike at degree 5 comes from products that list 5 items but are rarely
#    or never listed by others.
# The spike is therefore a data-collection artifact, not consumer behavior.
# edges[, 1] is the first column of the raw directed edge list, the product whose 
# page lists the other one. So table(edges[, 1]) counts how many products each 
# product lists, its out-degree. The crawl recorded at most 5, so the table stops at 5.
# The plot counts all links, in both directions. 
# Your degree distribution comes from G, the undirected graph. In an undirected graph, 
# a product’s degree counts:
#  - the products it lists (at most 5), plus
#  - the products that list it (no limit).

# Explains the spike at degree 5: crawl recorded at most 5 co-purchases per product
table(table(edges[, 1]))


##########################
# Line 138 of the RMD...
# power-law
##########################
# ---- Corrected log-log plot ---
# A log scale can't show 0, so drop degree 0 and any degree with proportion 0
keep = d > 0 & f_d > 0
plot(d[keep], f_d[keep], log = "xy",
     xlab = "Degree", ylab = "Proportion of products",
     main = "Log-log plot of degree distribution")
abline(v = 5, lty = 2, col = "grey50")   # cap artifact at degree 5 (spike)



# alpha_MLE = 3.5 here
# scale-free degree distr have alpha values ranging from 2 to 3
# we can conduct a GoF test with Kolmogorov-Smirnow (a HW had this) 
# if we want to but alpha_MLE = 3.5 is a good justification of
# our degree distr not being scale-free
plot(alpha_t, ylim=c(3, 4.5),
     xlab = "k(number of top products used)", ylab="Estimated alpha", main = "Hill plot")






