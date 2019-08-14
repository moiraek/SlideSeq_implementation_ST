#!/usr/bin/Rscript
#-----------------------------------------------------------------------
#
# Moira Ek
# 2019 07 18
# Direct implementation of the SlideSeq DGE method. The general idea 
#  is to, for each gene, first calculate the distribution of distances
#  between the spots in which it is expressed. In the paper by 
#  Rodriques et al. (Science 363, 1463–1467, 2019. 
#  doi 10.1126/science.aaw1219) 1000 random samples of the same number
#  of spots are then obtained, by sampling spots using probabilities
#  proportional to the total number of transcripts captured in each
#  spot, without replacement. For each such set of random spots, 
#  the distribution of distances between the spots is calculated, and 
#  the element-wise mean of the 1000 distributions is taken to obtain 
#  a mean distribution. The difference between the real distribution 
#  and the mean distribution as well as between each of the random 
#  samples and the mean distribution is calculated. The L1 norm of 
#  these difference vectors is then used as the test statistic, to 
#  determine whether the gene in question has a non-random spatial
#  distribution.
# 1. Initially implemented without any adaptations to ST, working with
#     raw count data, and using spot coordinates for the distances.
# 2. Optimization of the code, including not sampling anew for every
#     gene. The genes that are expressed in the same number of spots
#     are grouped together, and the same set of 1000 random sets
#     of spots are used for all of the genes within such a group.
# 
# 
#-----------------------------------------------------------------------

setwd("/home/moiraek/summerp19/SlideSeq_etc/Till_git/")

data <- as.data.frame(t(read.table("Rep1_MOB_count_matrix-1.tsv", check.names=FALSE)))

# Remove genes with expression in less than 5 spots.
testdata <- data[rowSums(data!=0)>4,]

#--------------------------------------------------------------------
# Create a distance matrix
#--------------------------------------------------------------------

start_time <- Sys.time()

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
# Calcuate the probability of sampling the spots. (Proportional to 
#  the total number of transcripts in each spot.) Indices are
#  calculated since the actual output from the sampling will be 
#  the index a specific spot has in the distance matrix. 
#--------------------------------------------------------------------

tot_spot <- colSums(testdata)
tot <- sum(tot_spot)
P <- vector(mode="numeric", length=ncol(testdata))
for (val in 1:length(P)){
  P[val] <- tot_spot[val]/tot
}
indices <- 1:length(P)


#--------------------------------------------------------------------
# Grouping of the genes according to the number of spots in which 
#  they are expressed, and thereby the number of spots that are to 
#  be sampled. 
#--------------------------------------------------------------------

# Calculate the maximum distance between included spots, in order
#  to define appropriate breakpoints for hist(). 0.5 due to the 
#  actual distances in the ST array.
maxdist <- max(euk)
breaks <- seq(0,maxdist+0.5, by=0.5)


# n = the number of random samples to draw for each gene.
# p = dataframe for the p-values.
# non_zero = the number of spots with non-zero expression of each gene.
n <- 1000
p <- vector(mode="numeric", length=nrow(testdata))
p <- as.data.frame(p, row.names=rownames(testdata), 
                   col.names="p_value")

non_zero <- matrix(nrow=nrow(testdata),ncol=1)
rownames(non_zero)<-rownames(testdata)
for (i in 1:nrow(non_zero)){
  non_zero[i,1] <- ncol(testdata[i, which(testdata[i,]!=0)])
}

# Group the genes expressed in the same number of spots on the rows of
#  equal_no_spots, each row is given the number of spots as its name.
#  unique --> each group of genes is only included once.
equal_no_spots <- matrix(NA,nrow=length(non_zero),
                         ncol=length(non_zero))
rownames(equal_no_spots)<-rownames(testdata)
for (j in 1:(nrow(non_zero))){
  equal<-rownames(non_zero)[which(non_zero==non_zero[j,1])]
  equal_no_spots[j,1:length(equal)] <- equal
}
rownames(equal_no_spots)<-as.character(non_zero[,1])

equal_no_spots<-unique(equal_no_spots)


#--------------------------------------------------------------------
# Go through the groups one by one. For each group the same number 
#  of indices as spots where these genes are expressed are first 
#  sampled. For each random sample, the counts in each bin, defined 
#  above, are saved as the rows of counts_matrix_rand, and this is
#  repeated n times (here, n=1000), until the matrix has been filled.
#  Column-wise means of this matrix gives the mean distribution,
#  mean_counts. The differences between each individual random 
#  distribution and the mean distribution are saved in diff_matr_rand.
#  The L1 norms of these differences are saved in L1_norms_rand.
#  The genes are then taken one by one. First, the indices of the
#  non-zero elements are obtained, i.e. the spots where the gene in
#  question is expressed. This is used to obtain the distances, which
#  are saved in eukl. The difference between the distribution of 
#  these distances (in the same bins as the random distributions)
#  and the mean distance are calculated, the L1 norm is calculated 
#  (L1_norm_real), and the p-value is calculated as 
#  p=(no random samples with L1>L1(true sample))/(total no random 
#  samples). If the numerator is equal to 0, it can however only
#  be said that p<1/n. These cases are currently saved as p=1/(10n).
#  In the supplementary material to the paper by Rodriques et al. 
#  (Science 363, 1463–1467, 2019. doi 10.1126/science.aaw1219)
#  it is stated that p=(no random samples with L1<L1(true sample))/
#  (total no random samples), but this contradicts their 
#  argumentation, as well as the example in figure S10.
#--------------------------------------------------------------------

for (i in 1:nrow(equal_no_spots)){
  # The counts for the genes in question
  values <- testdata[na.omit(equal_no_spots[i,]),]
  
  # The number of distances required
  len <- 0
  for (l in 1:(as.numeric(rownames(equal_no_spots)[i])-1)){
    len <- len + l
  }
  
  # The random distributions
  eukl_for_rand <- vector(mode="numeric", length=len)
  
  n_random <- as.numeric(rownames(equal_no_spots)[i])
  counts_matrix_rand <- matrix(0, nrow = n, ncol = (length(breaks)-1))
  for (k in 1:n){
    random_indices <- sort(sample(indices, size=n_random, 
                                  replace=FALSE, prob=P))
    interm_eukl <- euk[random_indices,random_indices]
    eukl_for_rand <- interm_eukl[lower.tri(interm_eukl, diag=FALSE)]
    counts_matrix_rand[k,] <- hist(eukl_for_rand, breaks, 
                                   plot=FALSE)$counts
  }
  
  mean_counts <- colMeans(counts_matrix_rand)
  
  # Element-wise difference between the distance distributions of 
  #  the random samples and the mean distribution.
  diff_matr_rand <- matrix(0, nrow=n, ncol=(length(breaks)-1))
  diff_matr_rand <- t(apply(counts_matrix_rand,1,'-',mean_counts))
  
  # The absolute values of these differences are summed over each row,
  #  yielding the L1 norms for the distances between each random 
  #  distribution and the mean distribution.
  abs_diff_rand <- abs(diff_matr_rand)
  L1_norms_rand <- rowSums(abs_diff_rand)
  
  
  # The true distributions.
  eukl <- vector(mode="numeric", length=len)
  vals <- vector(mode="numeric", 
                 length=as.numeric(rownames(equal_no_spots)[i]) )
 
  for (j in 1:nrow(values)){
    vals <- which(values[j,]!=0)
    eukli <- euk[vals,vals]
    eukl <- eukli[lower.tri(eukli, diag=FALSE)]
    
    eukdistr <- hist(eukl, breaks, plot=FALSE)
    
    diff_real <- eukdistr$counts - mean_counts
    abs_diff_real <- abs(diff_real)
    L1_norm_real <- sum(abs_diff_real)
    
    over_L1 <- L1_norms_rand[which(L1_norms_rand>L1_norm_real)]
    
    # Calculation of p values.
    if (length(over_L1)!=0){
      p[rownames(values)[j],1] <- length(over_L1)/n
    } else {
      p[rownames(values)[j],1] <- 1/(10*n)
    }
  }
}

diff_expr <- rownames(testdata[which(p<0.005),])
print(Sys.time() - start_time)


# --------------------------------------------------------------------
# Quick plotting, as a test
# --------------------------------------------------------------------
library(ggplot2)

gene <- "2010300C02Rik"
col <- as.numeric(as.vector(testdata[which(rownames(testdata)==gene),]))
# Create a colour gradient
rbPal <- colorRampPalette(c('yellow','red'))
color_vector <- rbPal(10)[as.numeric(cut(col,breaks = 10))]

xcoords <- as.numeric(sapply(strsplit(colnames(data), "x"), "[[", 1))
ycoords <- as.numeric(sapply(strsplit(colnames(data), "x"), "[[", 2))

# Plot the spot array with colours according to the gradient.
plot(x=xcoords, y=ycoords, col=alpha(color_vector, 1), lwd=1, asp=1,
     ylab="", xlab="", main=paste(gene), pch=19, cex.main=1.5, 
     xaxt="n", yaxt="n", bty="n", col.main="black")
