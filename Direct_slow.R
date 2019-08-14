#!/bin/Rscript
#-----------------------------------------------------------------------
#
# Moira Ek
# 2019 07 12
# Attempt to implement the SlideSeq DGE method.
# 1. Direct, without any adaptations to ST, initially working with
#    raw count data. Using spot coordinates for distance calculations.
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
# Functions:
#--------------------------------------------------------------------


# Create a function to give a random set of spots, with probabilities
#  according to the number of transcripts at each spot. Without
#  replacement. STdata is a dataframe with spots as columns, genes as
#  rows (the raw data), nr is the desired number of spots. Nr can never
#  exceed the total number of spots in the intended application --> no
#  if statement is needed.

SpotProb <- function(STdata, nr){
  tot_spot <- colSums(STdata)
  tot <- sum(tot_spot)
  P <- vector(mode="numeric", length=ncol(STdata))
  for (val in 1:length(P)){
    P[val] <- tot_spot[val]/tot
  }
  out <- sample(colnames(STdata), size=nr, replace=FALSE, prob=P)
  return(out)
}


# Function to give a distribution of the pairwise Euclidean distances 
#  between a given set of spots, spot_coords is a vector of spot names
#  of the form "x coordinate x y coordinate". For each spot the 
#  Euclidean distance to all spots to its right in the vector is 
#  calculated (--> the distance for each pair of spots is only included
#  once). Output: the distances.

distrib <- function(spot_coords){
  len <- 0
  for (val in 1:(length(spot_coords)-1)){
    len <- len + val
  }
  euk <- vector(mode="numeric", length=len)
  xcoords <- as.numeric(sapply(strsplit(colnames(spot_coords), "x"), "[[", 1))
  ycoords <- as.numeric(sapply(strsplit(colnames(spot_coords), "x"), "[[", 2))
  count <- 1
  for (j in 1:(length(spot_coords)-1)){
    rj <- c(xcoords[j], ycoords[j])
    for (k in (j+1):length(spot_coords)){
      rk <- c(xcoords[k], ycoords[k])
      euk[count] <- sqrt(sum((rj-rk)^2)) 
      count <- count+1
    }
  }
  return(euk)
}



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
# NB - what is an appropriate number of bins?
maxdist <- max(distrib(colnames(testdata)))
breaks <- seq(0,maxdist, by=(maxdist/200))


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
  values <- testdata[i,which(testdata[i,]!=0)]
  
  len <- 0
  for (val in 1:(ncol(values)-1)){
    len <- len + val
  }
  eukl <- vector(mode="numeric", length=len)
  
  spot_names <- colnames(values)
  eukl <- distrib(spot_names)
  
  eukdistr <- hist(eukl, breaks, plot=FALSE)
  # breaks <- eukdistr$breaks
  counts_matrix <- matrix(0, nrow = (n + 2), ncol = (length(breaks)-1))
  counts_matrix[1,] <- eukdistr$counts
  
  # The same number of distances for each random set of spots for each
  #  gene as for the actual spots where it is expressed. n random 
  #  distance distributions are saved (as the number of counts in each
  #  bin) on row 2 to n+1 in count_matrix. In eukl_for_rand the actual
  #  distances are saved.
  
  eukl_for_rand <- matrix(0, nrow = n, ncol = length(eukl))
  
  for (k in 1:n){
    random_spots <- SpotProb(testdata,length(values))
    eukl_for_rand[k,] <- as.numeric(distrib(random_spots))
    counts_matrix[(k+1),] <- hist(eukl_for_rand[k,], breaks, 
                                  plot=FALSE)$counts
  }
  
  # Interpreting element-wise mean as taking the mean of the counts in 
  #  each bin
  counts_matrix[(n+2),] <- colMeans(counts_matrix[(2:(n+1)),])
  
  # Element-wise difference for the distance distribution of the GOI
  #  and the mean distribution, as well as for the random distributions
  #  and the mean.
  diff_matr <- matrix(0, (n+1), ncol(counts_matrix))
  diff_matr[1,] <- counts_matrix[1,]-counts_matrix[n+2,]
  for (m in 1:n){
    diff_matr[m+1,] <- counts_matrix[m+1,]-counts_matrix[n+2,]
  }
  
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


end_time <- Sys.time()