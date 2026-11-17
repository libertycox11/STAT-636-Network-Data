##########################
#Line 44 of the RMD: change to...
##########################
# On average, each product is linked to about 7 other products in the
# co-purchasing network, counting both the products on its own
# "Customers Also Bought" list and the products whose lists include it.
# so not necessarily "its 'Customers Also Bought' list
# contains roughly 7 items"
# since degree counts both directions:
# products on its own “Also Bought” list (out-links), and
# products whose lists include it (in-links).

##########################
# Line 91 of the RMD...
# Degree distribution
##########################
f_d = degree_distribution(G)
max_d = length(f_d)
x_upper = range(degree(G))
barplot(f_d)
plot(f_d,  log='xy')

# check out the highest degree and understand the degree distribution
max(degree(G))        # the highest degree
sum(degree(G) == 420) # how many products have exactly 420 links
sum(degree(G) > 10)   # number of products with more than 10 links
mean(degree(G) > 10)  # their share of all products
median(degree(G))     # 50% of the products have degree <= 6! (skewness sign)

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
# Line 138 of the RMD...
# power-law
##########################
# alpha_MLE = 3.5 here
# scale-free degree distr have alpha values ranging from 2 to 3
# we can conduct a GoF test with Kolmogorov-Smirnow (a HW had this) 
# if we want to but alpha_MLE = 3.5 is a good justification of
# our degree distr not being scale-free
plot(alpha_t, ylim=c(3, 4.5),
     xlab = "k(number of top products used)", ylab="Estimated alpha", main = "Hill plot")







