#!/usr/bin/Rscript
#-------------------------------------------------------------------------
# Moira Ek
# 2019 08 08
#
# Attempt to adapt the SlideSeq_with_oversampling_limited.R method to
#  accept multiple samples as input.
# For each gene, the distances (both of the real distribution and
#  within the 1000 random samples) are calculated separately for
#  each sample, and these are then combined into one distance
#  distribution, from which the method continues as if working with
#  one sample.
#
# N.B. - Är det lämpligt att kombinera slumpproverna från olika prover
#  ett och ett? Hur hantera gener som filtreras bort i vissa prover,
#  ej i andra - for-loop, men bör någon typ av varning ges om resultaten
#  för en gen inte grundar sig på alla prover?
# Önskad input: vekror med paths till countmatriser för de olika 
#  proverna. Ut: en lista av SV gener och deras p-värden, och 
#  optionally även p-värden för alla gener(?)
# Problem: för att inte behöva spara undan alla fördelningar för alla
#  slumpfördelningar för alla samples, etc, kan inte generna grupperas,
#  varje gen måste behandlas för sig. Försök i nästa steg minska
#  påverkan av detta på tidsåtgången (mycket stor) genom att 
# parallellisera koden.
#-----------------------------------------------------------------------

setwd("/home/moiraek/summerp19/SlideSeq_etc/Till_git")

samples <- c("/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep1_MOB_count_matrix-1.tsv",
          "/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep2_MOB_count_matrix-1.tsv")

# Create a list of the data frames, one for each sample
data <- list()
data <- lapply(samples, function(x) as.data.frame(t(read.table(x, check.names=FALSE))))


# Remove genes with expression in less than 10 spots, and spots with
#  expression of less than 200 genes.
testdata <- lapply(data, function(x) x[(rowSums(x!=0)>9),])
testdata <- lapply(testdata, function(x) x[,colSums(x!=0)>199])

# Normalize the counts using SCTransform
library(Seurat)
se <- lapply(testdata, CreateSeuratObject)
se <- lapply(se, SCTransform)

testdata <- lapply(se, function(x) GetAssayData(x, slot="counts", assay="SCT"))
testdata <- lapply(testdata, as.matrix)

#--------------------------------------------------------------------
# Create a list of distance matrices, one for each sample
#--------------------------------------------------------------------

start_time <- Sys.time()

euk <- list()
for (sample in 1:length(testdata)){
  testdat <- testdata[[sample]]
  spotnames <- colnames(testdat)
  eukli <- matrix(0, ncol=ncol(testdat), nrow=ncol(testdat))
  rownames(eukli) <- spotnames
  colnames(eukli) <-spotnames
  xcoords <- as.numeric(sapply(strsplit(colnames(testdat), "x"), 
                               "[[", 1))
  ycoords <- as.numeric(sapply(strsplit(colnames(testdat), "x"), 
                               "[[", 2))
  
  for (j in 1:(ncol(testdat) - 1)){
    rj <- c(xcoords[j], ycoords[j])
    for (k in (j + 1):ncol(testdat)){
      rk <- c(xcoords[k], ycoords[k])
      dist <- sqrt(sum((rj-rk)^2))
      eukli[spotnames[j],spotnames[k]]<-dist
      eukli[spotnames[k],spotnames[j]]<-dist
    }
  }
  
  euk[[sample]] <- eukli
}

# Calculate the maximum distance between included spots, in order
#  to define appropriate breakpoints for hist(). 0.5 due to the 
#  actual distances in the ST array.
maxdist <- sapply(euk, max)
maxdist <- max(maxdist)
breaks <- seq(0,maxdist+0.5, by=0.5)

#--------------------------------------------------------------------
# Calcuate the probability of sampling the spots. (Proportional to 
#  the total number of transcripts in each spot.) Indices are
#  calculated since the actual output from the sampling will be 
#  the index a specific spot has in the distance matrix.
# This is done for each sample, creating a list of probabilities
#  and indices.
#--------------------------------------------------------------------

P <- list()
indices <- list()
tot_spot <- list()
tot <- list()
for (sample in 1:length(testdata)){
  testdat <- testdata[[sample]]
  tot_spot[[sample]] <- colSums(testdat)
  tot[[sample]] <- sum(tot_spot[[sample]])
  P[[sample]] <- vector(mode="numeric", length=ncol(testdat))
  for (val in 1:length(P[[sample]])){
    P[[sample]][val] <- tot_spot[[sample]][val]/tot[[sample]]
  }
  
  indices[[sample]] <- 1:ncol(testdat)
}

# Based on the number of transcripts in the spots in which it is 
#  expressed, an oversampling factor is calculated for each gene,
#  and these are scaled in such a way that 99% of the factors 
#  are less than or equal to 1 (since this is the range where
#  the method produces the best results).
# This is done for each gene in each sample.

factor_seur <- list()
for (sample in 1:length(testdata)){
  testdat <- testdata[[sample]]
  factor_seur[[sample]] <- vector(mode='numeric', length=nrow(testdat))
  for (i in 1:nrow(testdat)){
    factor_seur[[sample]][i] <- median(testdat[i,testdat[i,]!=0])
  }
  scale_factor <- quantile(factor_seur[[sample]], 0.99)
  
  factor_seur[[sample]] <- factor_seur[[sample]]/scale_factor
}

# n = the number of random samples to draw for the background 
#  distribution.
# p = dataframe for the p-values.
# non_zero = a list of, for each sample, the number of spots with 
#  non-zero expression of each gene.
n <- 1000

genes <- list()
for (sample in 1:length(testdata)){
  genes[[sample]] <- rownames(testdata[[sample]])
}

all_genes <- unique(unlist(genes))
p <- matrix(0, nrow=length(all_genes), ncol=2)
p <- as.data.frame(p, row.names=all_genes)
colnames(p) <- c("p", "no_samples")

non_zero <- list()
non_zero_with_oversampl <- list()
for (sample in 1:length(testdata)){
  testdat <- testdata[[sample]]
  non_zero[[sample]] <- matrix(nrow=nrow(testdat),ncol=1)
  rownames(non_zero[[sample]])<-rownames(testdat)
  for (i in 1:nrow(non_zero[[sample]])){
    non_zero[[sample]][i,1] <- length(testdat[i, which(testdat[i,]!=0)])
  }
  non_zero_with_oversampl[[sample]] <- floor(non_zero[[sample]]*(1+factor_seur[[sample]]))
}

#--------------------------------------------------------------------
# Go through the groups one by one. For each group the same number 
#  of indices as spots where these genes are expressed are first 
#  sampled. (In this case, with replacement.) For each random sample,
#  the counts in each bin, defined above, are saved as the rows of 
#  counts_matrix_rand, and this is repeated n times (here, n=1000), 
#  until the matrix has been filled.
#  Column-wise means of this matrix gives the mean distribution,
#  mean_counts. The differences between each individual random 
#  distribution and the mean distribution are saved in diff_matr_rand.
#  The L1 norms of these differences are saved in L1_norms_rand.
#  The genes are then taken one by one. First, the indices of the
#  non-zero elements are obtained, i.e. the spots where the gene in
#  question is expressed. Then, an additional number of spots is 
#  sampled, based on the conditional probabilities, in the form of
#  indices. These are then used to obtain the distances, which
#  are saved in eukl. The difference between the distribution of 
#  these distances (in the same bins as the random distributions)
#  and the mean distance are calculated, the L1 norm is calculated 
#  (L1_norm_real), and the p-value is calculated as 
#  p=(no random samples with L1>=L1(true sample))/(total no random 
#  samples). If the numerator is equal to 0, it can however only
#  be said that p<1/n. These cases are currently saved as p=1/(10n).
#  In the supplementary material to the paper by Rodriques et al. 
#  (Science 363, 1463–1467, 2019. doi 10.1126/science.aaw1219)
#  it is stated that p=(no random samples with L1<L1(true sample))/
#  (total no random samples), but this contradicts their 
#  argumentation, as well as the example in figure S10.
#--------------------------------------------------------------------

for (i in 1:length(all_genes)){
  
  eukdistr <- list()
  counts_matrix_rand <- list()
  included_samples <- length(testdata)
  for (sample in 1:length(testdata)){
    testdat <- testdata[[sample]]
    values <- vector(mode="numeric", length=ncol(testdat))
    if (identical(which(rownames(testdat)==all_genes[i]), integer(0))){
      # If the gene is not present in this sample, no distances should
      #  be generated, and the values vector remains a vector of zeroes
      included_samples <- included_samples - 1
      n_spots <- 0
      
      #The number of distances required
      len <- 0
    
    } else {
      # The counts for the gene in question, if it is included
      values <- testdat[all_genes[i],]
      
      # The number of spots required.
      n_spots <- as.numeric(non_zero_with_oversampl[[sample]][all_genes[i],1])
     
       # The number of distances required
      len <- 0
      for (l in 1:(n_spots-1)){
        len <- len + l
      }
      
    }
  
    # The true distribution
    eukl <- vector(mode="numeric", length=len)
    eukli <- matrix(0, ncol=n_spots, nrow=n_spots)
  
    # The conditional probabilities for the gene in question
    P_g <- sum(values)/tot[[sample]]
    P_g_s <- vector(mode='numeric', length=length(values))
    P_cond <- vector(mode='numeric', length=length(values))
    for (l in 1:length(values)){
      P_g_s[l] <- as.numeric(values[l])/tot_spot[[sample]][l]
      P_s <- P[[sample]][l]
    }
  
   if (P_g==0){
      # Vector of zeroes, as to not yield an error message
      P_cond <- vector(mode="numeric", length=length(values))
    } else{
      P_cond <- P_g_s*P_s/P_g
    }
    
    if (identical(which(rownames(testdat)==all_genes[i]), integer(0))){
      n_expressed <- 0
      n_oversampling <- 0
    } else {
      n_expressed <- as.numeric(non_zero[[sample]][all_genes[i],1])
      n_oversampling <- n_spots - n_expressed
    }

    expressed_indices <- which(values!=0)
    oversampling_indices <- sample(indices[[sample]], size=n_oversampling, 
                                 replace=TRUE, prob=P_cond)
    gene_indices <- sort(c(expressed_indices,oversampling_indices))
    eukli <- euk[[sample]][gene_indices,gene_indices]
    eukl <- eukli[lower.tri(eukli, diag=FALSE)]
    eukdistr[[sample]] <- hist(eukl, breaks, plot=FALSE)$counts
  

    # The random distributions for each sample, saved in a list
    eukl_for_rand <- vector(mode="numeric", length=len)
  
    counts_matrix_rand[[sample]] <- matrix(0, nrow = n, ncol = (length(breaks)-1))
    for (k in 1:n){
      random_indices_1 <- sample(indices[[sample]], size=n_expressed, 
                                 replace=FALSE, prob=P[[sample]])
      random_indices_2 <- sample(indices[[sample]], size=n_oversampling, 
                                 replace=TRUE, prob=P[[sample]])
      random_indices <- sort(c(random_indices_1,random_indices_2))
      interm_eukl <- euk[[sample]][random_indices,random_indices]
      eukl_for_rand <- interm_eukl[lower.tri(interm_eukl, diag=FALSE)]
      counts_matrix_rand[[sample]][k,] <- hist(eukl_for_rand, breaks, 
                                          plot=FALSE)$counts
    }
    # For each gene, the number of samples that were actually included
    #  in the calculations is given.
    p[all_genes[i],2]<- included_samples

  }
  
  
  #Combining the distributions for each sample, both true and random:
  # True:
  eukdistrib <- colSums(do.call(rbind, eukdistr))
  
  # Random:
  rand_distrib <- matrix(0, nrow = n, ncol = (length(breaks)-1))
  for (m in 1:n){
    # Extracting the m:th row from each matrix in the list, taking the
    #  sum of the number of elements in each bin
    rand_distrib[m,] <- colSums(do.call(rbind, lapply(counts_matrix_rand, `[`,m,)))
  }
  
  mean_counts <- colMeans(rand_distrib)
  
  # Element-wise difference between the distance distributions of 
  #  the random samples and the mean distribution, as well as
  #  between the "true" distribution and the mean random distribution
  diff_matr_rand <- matrix(0, nrow=n, ncol=(length(breaks)-1))
  diff_matr_rand <- t(apply(rand_distrib, 1, '-', mean_counts))
  diff_real <- eukdistrib - mean_counts
  
  # The absolute values of these differences are summed over each row,
  #  yielding the L1 norms for the distances between each random 
  #  distribution/"true" distribution and the mean distribution.
  abs_diff_rand <- abs(diff_matr_rand)
  L1_norms_rand <- rowSums(abs_diff_rand)
  
  abs_diff_real <- abs(diff_real)
  L1_norm_real <- sum(abs_diff_real)
    
  over_L1 <- L1_norms_rand[which(L1_norms_rand>=L1_norm_real)]
    
  # Calculation of p value
  if (length(over_L1)!=0){
    p[all_genes[i],1] <- length(over_L1)/n
  } else {
    p[all_genes[i],1] <- 1/(10*n)
  }
  cat(c(i, " out of ", length(all_genes), " finished\n"))
}

# In the paper by Rodriques et al. p-values of less than 0.005 were
#  deemed significant.
diff_expr <- all_genes[which(p<0.005),]

# Calculating adjusted p values in order to correct for multiple
#  hypothesis testing, using the Benjamini-Hochberg procedure.
#  N.B. - this requires that the tests are independent or positively
#  dependent, not sure if this can really be said to be the case?
#  Otherwise, the Bejamini-Yekutieli procedure can be used.
# p-values saved as 1/(10n) are first converted to 1/n, as it is 
#  better to be slightly too conservative in this case.
# Accepting a FDR of 5%, all genes with an adjusted p value of less 
#  than 0.05 are said to be spatially variable.

p2 <- p[,1]
for (val in which(p[,1]==0.0001)){
  p2[val,1] <- 0.001
}

p_adj <- p.adjust(p2[,1], method="BH")
p_adj <- as.data.frame(p_adj, row.names=rownames(p2))
colnames(p_adj) <- "Adjusted p"
p_adj$no_samples <- p[,2]

diff <- which(p_adj[,1]<0.05)
p_adj_diff <- as.data.frame(p_adj[diff,], row.names=rownames(p_adj)[diff])
diff_expr_MHT <- rownames(p_adj_diff)

print(Sys.time() - start_time)

# --------------------------------------------------------------------
# Quick plotting, as a test
# --------------------------------------------------------------------
library(ggplot2)

gene <- "Champ1"

dat <- data[which(rownames(data)==gene),]
dat[1,which(dat>quantile(dat,0.99)[1,1])]<-quantile(dat,0.99)[1,1]

col <- as.numeric(as.vector(dat))
# Create a colour gradient
rbPal <- colorRampPalette(c('yellow','red'))
color_vector <- rbPal(10)[as.numeric(cut(col,breaks = 10))]

xcoords <- as.numeric(sapply(strsplit(colnames(data), "x"), "[[", 1))
ycoords <- as.numeric(sapply(strsplit(colnames(data), "x"), "[[", 2))

# Plot the spot array with colours according to the gradient.
plot(x=xcoords, y=ycoords, col=alpha(color_vector, 1), lwd=1, asp=1,
     ylab="", xlab="", main=paste(gene), pch=19, cex.main=1.5,
     xaxt="n", yaxt="n", bty="n", col.main="black")

par(mfrow=c(5,4))
par(mar=c(1,1,1,1))
for (gene in diff_expr_MHT[21:40]){
  dat <- data[which(rownames(data)==gene),]
  dat[1,which(dat>quantile(dat,0.99)[1,1])]<-quantile(dat,0.99)[1,1]
  
  col <- as.numeric(as.vector(dat))
  # Create a colour gradient
  rbPal <- colorRampPalette(c('yellow','red'))
  color_vector <- rbPal(10)[as.numeric(cut(col,breaks = 10))]
  
  xcoords <- as.numeric(sapply(strsplit(colnames(data), "x"), "[[", 1))
  ycoords <- as.numeric(sapply(strsplit(colnames(data), "x"), "[[", 2))
  
  # Plot the spot array with colours according to the gradient.
  #dev.new()
  plot(x=xcoords, y=ycoords, col=alpha(color_vector, 1), lwd=1, asp=1,
       ylab="", xlab="", main=paste(gene, round(p_adj_diff[gene,1], digits=4)), pch=19, cex.main=1,
       xaxt="n", yaxt="n", bty="n", col.main="black")
}
