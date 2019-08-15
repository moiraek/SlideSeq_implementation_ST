#!/usr/bin/Rscript
#-----------------------------------------------------------------------
#
# Moira Ek
# 2019 07 17
# Attempt to implement the SlideSeq DGE method
# 1. Direct, without any adaptations to ST, initially working with raw
#    counts. Using spot coordinates for distance calculations.
# 2. Attempt to optimize the code in order to decrease computational 
#    time.
# 
#-----------------------------------------------------------------------

#----------------------------------------------------------------------
# Preparations, one sample is analyzed at a time. Initially looking
#  at a very small number of genes, due to long computation times
#----------------------------------------------------------------------
setwd("/home/moiraek/summerp19/SlideSeq_etc/Till_git")

data <- as.data.frame(t(read.table("Rep1_MOB_count_matrix-1.tsv", check.names=FALSE)))

testdata <- data[1:10,]

#--------------------------------------------------------------------
# Create a distance matrix
#--------------------------------------------------------------------

spotnames <- colnames(testdata)
euk <- matrix(0, ncol=ncol(testdata), nrow=ncol(testdata))
rownames(euk) <- spotnames
colnames(euk) <-spotnames

xcoords <- as.numeric(sapply(strsplit(colnames(testdata), "x"), "[[", 1))
ycoords <- as.numeric(sapply(strsplit(colnames(testdata), "x"), "[[", 2))

for (j in 1:(ncol(testdata) - 1)){
  rj <- c(xcoords[j], ycoords[j])
  for (k in (j + 1):ncol(testdata)){
    rk <- c(xcoords[k], ycoords[k])
    dist <- sqrt(sum((rj-rk)^2))
    euk[spotnames[j],spotnames[k]]<-dist
    euk[spotnames[k],spotnames[j]]<-dist
  }
}


#--------------------------------------------------------------------
# Calculate probabilities of sampling the spots
#--------------------------------------------------------------------

tot_spot <- colSums(testdata)
tot <- sum(tot_spot)
P <- vector(mode="numeric", length=ncol(testdata))
for (val in 1:length(P)){
  P[val] <- tot_spot[val]/tot
}
indices <- 1:length(P)


#--------------------------------------------------------------------
# Attempt
#--------------------------------------------------------------------

# Calculate the maximum distance between included spots, in order to
#  set appropriate breakpoints for hist. Before, the true distribution 
#  was used to set the bin breakpoints (set a number of bins for it,
#  used $breaks as input to hist for the random distributions), but
#  this entails a risk that all distances in the random samples aren't
#  covered, e.g. if a gene is expressed in a limited region. (Since
#  those spots are sampled from the whole array.)
#  0.5 due to the actual distances in the arrays.
maxdist <- max(euk)
breaks <- seq(0,maxdist+1, by=0.5)



# For each gene, take the spots with non-zero counts ("values"), 
#  calculate the distances. eukl contains the distances between the 
#  spots withexpression of each gene, followed by zeroes. First create
#  a matrix for the counts of (number of bins) columns, (n+2) rows - 
#  one distribution for the true distances, n random distributions, and
#  one mean random distribution. n= the number of random sets of spots 
#  to sample, SlideSeq uses n = 1000.
# Create eukl - a vector for the distances, of length (x-1)+(x-2)+...+1
#  where x is the number of spots with expression of the gene in 
#  question. p - to save the p value for each gene.
n <- 1000
p <- vector(mode="numeric", length=nrow(testdata))

start_time <- Sys.time()
for (i in 1:nrow(testdata)){
  non_zero <- which(testdata[i,]!=0)
  
  eukl <- euk[non_zero,non_zero]
  eukl <- eukl[lower.tri(eukl, diag=FALSE)]
  
  eukdistr <- hist(eukl, breaks, plot=FALSE)
  counts_matrix <- matrix(0, nrow = (n + 1), ncol = (length(breaks)-1))
  counts_matrix[1,] <- eukdistr$counts
  
  # The same number of distances for each random set of spots for each
  #  gene as for the actual spots where it is expressed. n random 
  #  distance distributions are saved (as the number of counts in each
  #  bin) on row 2 to n+1 in count_matrix. In eukl_for_rand the actual
  #  distances are saved.
  
  eukl_for_rand <- vector(mode="numeric", length=length(eukl))
  
  n_random <- length(non_zero)
  
  for (k in 1:n){
    random_indices <- sort(sample(indices, size=n_random, replace=FALSE, prob=P))
    interm_eukl <- euk[random_indices,random_indices]
    eukl_for_rand <- interm_eukl[lower.tri(interm_eukl, diag=FALSE)]
    counts_matrix[(k+1),] <- hist(eukl_for_rand, breaks, 
                                  plot=FALSE)$counts
  }
  
  # Interpreting element-wise mean as taking the mean of the counts in 
  #  each bin
  mean_counts <- colMeans(counts_matrix[(2:(n+1)),])
  
  # Element-wise difference for the distance distribution of the GOI
  #  and the mean distribution, as well as for the random distributions
  #  and the mean.
  diff_matr <- matrix(0, (n+1), ncol(counts_matrix))
  diff_matr <- t(apply(counts_matrix,1,'-',mean_counts))
  
  # Take the absolute value of these, and sum over each row --> the
  #  L1-norm of the distance between the true distribution (element 1)
  #  and the random distributions (element 2 to (n+1)), respectively, 
  #  and the mean random distribution. 
  abs_diff <- abs(diff_matr)
  L1_norms <- rowSums(abs_diff)
  #hist(L1_norms[2:(n+1)], n/2, xlim=c(min(L1_norms),max(L1_norms)))
  #abline(v=L1_norms[1], col="red")
  
  # The smallest possible p value is p < 1/n, with n being the number
  #  of random sets of spots.In general, p = (no. random samples with
  #  L1>L1(true sample))/(total no. random samples). If the numerator
  #  is 0, it can however only be said that p < 1/n. In the supplementary
  #  material to the paper it is said that p = (no. random samples with
  #  L1<L1(true sample))/(total no. random samples), but this contradicts
  #  the authors' argument, and the example in figure S10. For p < 1/n:
  #  for the time being, the value is saved as p = (1/10n), to set it apart
  #  from actual p values of 1/n. Can this be done in a nicer way?
  
  over_L1 <- L1_norms[which(L1_norms>L1_norms[1])]
  if (length(over_L1)!=0){
    p[i] <- length(over_L1)/(length(L1_norms)-1)
  } else {
    p[i] <- 1/(10*n)
  }
  
  
}

p <- as.data.frame(p, row.names=rownames(testdata), col.names="p-value")


diff_expr <- rownames(testdata[which(p<0.005),])
print(Sys.time() - start_time)

end_time <- Sys.time()


# --------------------------------------------------------------------
# Quick plotting, test
# --------------------------------------------------------------------
library(ggplot2)

gene <- "Chrm1"
col <- as.numeric(as.vector(testdata[which(rownames(testdata)==gene),]))
# Create a colour gradient
rbPal <- colorRampPalette(c('yellow','red'))
color_vector <- rbPal(10)[as.numeric(cut(col,breaks = 10))]

coords <- GetCoords(colnames(data), delim="x")
xcoords <- coords$x
ycoords <- coords$y

# Plot the array with spot colours according to the gradient
plot(x=xcoords, y=ycoords, col=alpha(color_vector, 1), lwd=1, asp=1,
     ylab="", xlab="", main=paste(gene), pch=19, cex.main=1.5, 
     xaxt="n", yaxt="n", bty="n", col.main="black")
