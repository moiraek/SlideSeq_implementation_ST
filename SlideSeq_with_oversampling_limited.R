#!/usr/bin/Rscript
#-------------------------------------------------------------------------
# Moira Ek
# 2019 08 01
#
# Adaptation of the SlideSeq DGE method described by by Rodriques et 
#  al. (Science 363, 1463–1467, 2019. doi 10.1126/science.aaw1219) to
#  ST data. For each gene, each spot with expression is included once, 
#  and in addition, a number of spots (depending on the expression level
#  of the gene in question, but in most cases smaller than the number
#  of spots in which expression is observed) is sampled from the total
#  set of spots. For gene i, the probability of sampling spot j is the 
#  conditional probability P(spot(j)|gene(i)), given by Baye's theorem.
#  The obtained total set of spots is used to calculate a distribution 
#  of pairwise distances. As a background, the same number of spots is
#  sampled from the total set of spots, with probabilities given 
#  by the total number of transcripts in each spot divided by the 
#  total number of transcripts in all spots (As of now: with 
#  replacement. Should we first without replacement sample the same
#  number of spots as where the gene is expressed, and then the 
#  additional number with replacement? Would be more correct, but
#  take much more time as the genes will barely be allowed to be 
#  grouped anymore... How big an impact would it have?). 
#  n (1000 in the paper by Rodriques et al.) such sets are obtained,
#  their distance distributions are calculated, and the mean of these
#  distributions is taken. The distribution obtained for the gene 
#  under study is compared to this, and L1 norms are obtained and 
#  compared as in the original method, see e.g. SlideSeq_direct_with_MHT.R 
#  for a description.


# N.B. - Bakgrunden samplas nu helt och hållet med återläggning för att
#  kunna gruppera generna (--> beroende på antalet spots kan t.ex. ca
#  800 bakgrundsfördelningsset dras istället för drygt 13 000, vid ca 
#  600 spots). Efter normaliseringen skiljer sig inte spotsen så mycket 
#  åt i totalt antal transkript, så det är möjligt att detta inte har
#  någon nämnbar effekt. Då ex har tillgång till severn: testa båda 
#  varianter och jämför, för att se om detta är motiverat.
#
#-----------------------------------------------------------------------

setwd("/home/moiraek/summerp19/SlideSeq_etc/Till_git")

data <- as.data.frame(t(read.table("Rep1_MOB_count_matrix-1.tsv", 
                                   check.names=FALSE)))
#data <- as.data.frame(read.table("hippocampus_wt_rep1.tsv", 
#                                   check.names=FALSE))

# Remove genes with expression in less than 10 spots, and spots with
#  expression of less than 200 genes.
testdata <- data[rowSums(data!=0)>9,]
testdata <- testdata[,colSums(testdata!=0)>199]

# Normalize the counts using SCTransform
library(Seurat)
se <- CreateSeuratObject(testdata)
se <- SCTransform(se)
testdat <- GetAssayData(se, slot="counts", assay="SCT")
testdat <- as.matrix(testdat)
testdata <- testdat

#--------------------------------------------------------------------
# Create a distance matrix
#--------------------------------------------------------------------

start_time <- Sys.time()

spotnames <- colnames(testdata)
euk <- matrix(0, ncol=ncol(testdata), nrow=ncol(testdata))
rownames(euk) <- spotnames
colnames(euk) <-spotnames

xcoords <- as.numeric(sapply(strsplit(colnames(testdata), "x"), 
                             "[[", 1))
ycoords <- as.numeric(sapply(strsplit(colnames(testdata), "x"), 
                             "[[", 2))

for (j in 1:(ncol(testdata) - 1)){
  rj <- c(xcoords[j], ycoords[j])
  for (k in (j + 1):ncol(testdata)){
    rk <- c(xcoords[k], ycoords[k])
    dist <- sqrt(sum((rj-rk)^2))
    euk[spotnames[j],spotnames[k]]<-dist
    euk[spotnames[k],spotnames[j]]<-dist
  }
}

# Calculate the maximum distance between included spots, in order
#  to define appropriate breakpoints for hist(). 0.5 due to the 
#  actual distances in the ST array.
maxdist <- max(euk)
breaks <- seq(0,maxdist+1, by=1)

# Capping of distances, so that single spots far away won't have too
#  great an impact
#cap <- maxdist*0.9
#euk[euk>cap]<-cap

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

indices <- 1:ncol(testdata)


# Based on the number of transcripts in the spots in which it is 
#  expressed, an oversampling factor is calculated for each gene,
#  and these are scaled in such a way that 99% of the factors 
#  are less than or equal to 1 (since this is the range where
#  the method produces the best results).
factor_seur <- vector(mode='numeric', length=nrow(testdata))
for (i in 1:nrow(testdata)){
  factor_seur[i] <- median(testdata[i,testdata[i,]!=0])
}
scale_factor <- quantile(factor_seur, 0.99)
factor_seur <- factor_seur/scale_factor


# n = the number of random samples to draw for the background 
#  distribution.
# p = dataframe for the p-values.
# non_zero = the number of spots with non-zero expression of each gene.
n <- 1000
p <- vector(mode="numeric", length=nrow(testdata))
p <- as.data.frame(p, row.names=rownames(testdata), 
                   col.names="p")

non_zero <- matrix(nrow=nrow(testdata),ncol=1)
rownames(non_zero)<-rownames(testdata)
for (i in 1:nrow(non_zero)){
  non_zero[i,1] <- length(testdata[i, which(testdata[i,]!=0)])
}

# N.B.- Is this appropriate? or ceil or floor?
non_zero_with_oversampl <- round(non_zero*(1+factor_seur))

# Group the genes expressed in the same number of spots on the rows of
#  equal_no_spots, each row is given the number of spots as its name.
#  unique --> each group of genes is only included once.
equal_no_spots <- matrix(NA,nrow=length(non_zero),
                         ncol=length(non_zero))
rownames(equal_no_spots)<-rownames(testdata)
for (j in 1:(nrow(non_zero_with_oversampl))){
  equal<-rownames(non_zero_with_oversampl)[which(non_zero_with_oversampl
                                                 ==non_zero_with_oversampl[j,1])]
  equal_no_spots[j,1:length(equal)] <- equal
}
rownames(equal_no_spots)<-as.character(non_zero_with_oversampl[,1])

equal_no_spots<-unique(equal_no_spots)


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

for (i in 1:nrow(equal_no_spots)){
  # The counts for the genes in question
  if (length(na.omit(equal_no_spots[i,]))==1){
    values <- t(as.matrix(testdata[na.omit(equal_no_spots[i,]),]))
    rownames(values) <- as.character(na.omit(equal_no_spots[i,]))
  } else {
    values <- testdata[na.omit(equal_no_spots[i,]),]
  }
  
  # The number of spots required.
  n_spots <- as.numeric(rownames(equal_no_spots)[i])

  # The number of distances required
  len <- 0
  for (l in 1:(n_spots-1)){
    len <- len + l
  }
  
  # The random distributions
  eukl_for_rand <- vector(mode="numeric", length=len)
  
  counts_matrix_rand <- matrix(0, nrow = n, ncol = (length(breaks)-1))
  for (k in 1:n){
    random_indices <- sort(sample(indices, size=n_spots, 
                                  replace=TRUE, prob=P))
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
  eukli <- matrix(0, ncol=n_spots, nrow=n_spots)
  
  
  for (j in 1:nrow(values)){
    # The conditional probabilities for the gene in question
    vals <- values[j,]
    P_g <- sum(vals)/tot
    P_g_s <- vector(mode='numeric', length=length(vals))
    P_cond <- vector(mode='numeric', length=length(vals))
    P_s <- vector(mode='numeric', length=length(vals))
    for (l in 1:length(vals)){
      P_g_s[l] <- as.numeric(vals[l])/tot_spot[l]
      P_s[l] <- P[l]
    }
    P_cond <- P_g_s*P_s/P_g
    
    n_expressed <- as.numeric(non_zero[rownames(values)[j],1])
    n_oversampling <- n_spots - n_expressed
    
    
    expressed_indices <- which(vals!=0)
    oversampling_indices <- sample(indices, size=n_oversampling, 
                                   replace=TRUE, prob=P_cond)
    gene_indices <- sort(c(expressed_indices,oversampling_indices))
    eukli <- euk[gene_indices,gene_indices]
    eukl <- eukli[lower.tri(eukli, diag=FALSE)]
    eukdistr <- hist(eukl, breaks, plot=FALSE)
    
    diff_real <- eukdistr$counts - mean_counts
    abs_diff_real <- abs(diff_real)
    L1_norm_real <- sum(abs_diff_real)
    
    over_L1 <- L1_norms_rand[which(L1_norms_rand>=L1_norm_real)]
    
    # Calculation of p values.
    if (length(over_L1)!=0){
      p[rownames(values)[j],1] <- length(over_L1)/n
    } else {
      p[rownames(values)[j],1] <- 1/(10*n)
    }
  }
  cat(c(i, " out of ", nrow(equal_no_spots), " finished\n"))
}

# In the paper by Rodriques et al. p-values of less than 0.005 were
#  deemed significant.
diff_expr <- rownames(testdata[which(p<0.005),])

#Alternatively, the Benjamini-Hochberg procedure can be used to
# control the FDR. However, according to Benjamini and Yekutieli
# (Annals of Statistics, 29(4), 1165-1188, 2001) this is only 
# applicable if the tests are independent or are positively
# dependent. Not entirely sure of the case here... If this is not
# the case, the Benjamini-Yekutieli procedure can be used, by
# including a coefficient, but due to the limited range of observable
# p-values here, this doesn't work.

# p2 <- p[order(p$p), , drop = FALSE]
# count <- nrow(p2)
# 
# q <- 0.05
# test <- matrix(NA, nrow=nrow(p2),ncol=1)
# rownames(test)<- rownames(p2)
# 
# 
# # For the cases of which only p<0.001 (saved as p=0.0001) can be said
# #  the procedure is only meaningful if these cases aren't the only
# #  ones that pass the test. Since it is however certain that they
# #  have p-values of less than 0.001, the below will never lead to any
# #  conclusions of null hypothesis rejection if this cannot certainly
# #  be supported.
# for (val in which(p2[,1]==0.0001)){
#   p2[val,1]<-0.001
# }
# 
# # B-Y procedure
# #c <- 0
# #for (j in 1:count){
# #  c <- c + 1/j
# #}
# #for (i in 1:count){
# #  test[i,1] <- (p2[i,1]<=((i/(count*c))*q))
# #}
# 
# for (i in 1:count){
#   test[i,1] <- (p2[i,1]<=((i/count)*q))
# }
# 
# suppressWarnings(k <- max(which(test[,1]==TRUE)))
# if (k==-Inf){
#   print("No null hypotheses can be rejected with certainty")
# } else{
#   diff_expr_MHT <- rownames(p2)[1:k]
# }


p2 <- p
for (val in which(p[,1]==0.0001)){
  p2[val,1]<-0.001
}

p_adj <- p.adjust(p2[,1], method="BH")
p_adj <- as.data.frame(p_adj, row.names=rownames(p2), 
                       col.names="adjusted p")
diff <- which(p_adj<0.05)
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
for (gene in diff_expr_MHT[1:20]){
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
