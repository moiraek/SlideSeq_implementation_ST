#!/usr/bin/Rscript
#-------------------------------------------------------------------------
# Moira Ek
# 2019 08 06
#
# Adaptation of the SlideSeq_with_oversampling_limited.R method to 
#  accept a vector of gene names as input (in addition to a count 
#  matrix - if making this into a function: add a parameter for
#  whether the data is normalized or not, with value T/F), and give
#  a p-value for each gene as output, indicating whether or not it is
#  spatially variable.
#  
# N.B. - look at how to extend this to accept count matrices from 
#  multiple sections as input.
# N.B. - Fix so that a single gene also works

# Testar här anpassningen att köra generna var för sig, och ta
#  bakgrunden på samma vis som de faktiska proverna - först ett visst
#  antal utan återläggning, sedan ytterligare ett par med.
#
#-----------------------------------------------------------------------

# To have something to test with, and in case the data hasn't been 
#  normalized yet.
setwd("/home/moiraek/summerp19/SlideSeq_etc/Till_git")

samples <- c("/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep1_MOB_count_matrix-1.tsv",
             "/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep2_MOB_count_matrix-1.tsv",
             "/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep3_MOB_count_matrix-1.tsv",
             "/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep4_MOB_count_matrix-1.tsv",
             "/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep5_MOB_count_matrix-1.tsv",
             "/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep6_MOB_count_matrix-1.tsv",
             "/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep7_MOB_count_matrix-1.tsv",
             "/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep8_MOB_count_matrix-1.tsv",
             "/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep9_MOB_count_matrix-1.tsv",
             "/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep10_MOB_count_matrix-1.tsv",
             "/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep11_MOB_count_matrix-1.tsv",
             "/home/moiraek/summerp19/SlideSeq_etc/Till_git/Rep12_MOB_count_matrix-1.tsv")

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

# All of the included data should be used for calculating the 
#  probabilities of sampling certain spots
testdata_all <- testdata

# Looking only at specified genes, gene_vector has to be given as input
rm(.Random.seed, envir=globalenv())
gene_vector <- rownames(testdata_all[[1]])[round(runif(5, min=1, max=length(rownames(testdata_all[[1]]))))]

genes <- list()
for (sample in 1:length(testdata)){
  genes[[sample]] <- gene_vector[gene_vector %in% rownames(testdata[[sample]])]
  testdata[[sample]] <- testdata[[sample]][genes[[sample]],]
}

par(mfrow=c(length(gene_vector),(length(testdata)+1)))
par(mar=c(1,1,1,1))
library(ggplot2)

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
#  to define appropriate breakpoints for hist(). 1 due to the 
#  actual distances in the ST array.
maxdist <- sapply(euk, max)
maxdist <- max(maxdist)
breaks <- seq(0,maxdist+1, by=1)

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
  tot_spot[[sample]] <- colSums(testdata_all[[sample]])
  tot[[sample]] <- sum(tot_spot[[sample]])
  P[[sample]] <- vector(mode="numeric", length=ncol(testdata[[sample]]))
  for (val in 1:length(P[[sample]])){
    P[[sample]][val] <- tot_spot[[sample]][val]/tot[[sample]]
  }
  
  indices[[sample]] <- 1:ncol(testdata[[sample]])
}


# Based on the number of transcripts in the spots in which it is 
#  expressed, an oversampling factor is calculated for each gene,
#  and these are scaled in such a way that 99% of the factors 
#  are less than or equal to 1 (since this is the range where
#  the method produces the best results).
# This is done for each gene in each sample.

factor_seur <- list()
for (sample in 1:length(testdata)){
  testdat <- testdata_all[[sample]]
  factor_seur[[sample]] <- vector(mode='numeric', length=nrow(testdat))
  factor_seur[[sample]] <- as.data.frame(factor_seur[[sample]], 
                                         row.names=rownames(testdata_all[[sample]]))
  for (i in 1:nrow(testdat)){
    factor_seur[[sample]][i,1] <- median(testdat[i,testdat[i,]!=0])
  }
  scale_factor <- quantile(factor_seur[[sample]][,1], 0.99)
  
  factor_seur[[sample]][,1] <- factor_seur[[sample]][,1]/scale_factor
  factor_seur[[sample]] <- factor_seur[[sample]][genes[[sample]],]
}

# n = the number of random samples to draw for the background 
#  distribution.
# p = dataframe for the p-values.
# non_zero = a list of, for each sample, the number of spots with 
#  non-zero expression of each gene.
n <- 1000

# Alternativt bara gene_vector, vad sker i så fall om den innehåller en
#  gen som inte finns med i någon av proverna? - se mer på detta
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

gg <- list()
for (i in 1:length(all_genes)){
  
  eukdistr <- list()
  counts_matrix_rand <- list()
  values <- list()
  
  # For each gene, keep count of the number of samples that were
  #  included in the analysis.
  included_samples <- length(testdata)
  for (sample in 1:length(testdata)){
    testdat <- testdata[[sample]]
    values[[sample]] <- vector(mode="numeric", length=ncol(testdat))
    if (identical(which(rownames(testdat)==all_genes[i]), integer(0))){
      # If the gene is not present in this sample, no distances should
      #  be generated, and the values vector remains a vector of zeroes
      included_samples <- included_samples - 1
      
    } else {
      # The counts for the gene in question, if it is included
      values[[sample]] <- testdat[all_genes[i],]
      
      # The number of spots required.
      n_spots <- as.numeric(non_zero_with_oversampl[[sample]][all_genes[i],1])
      
      # The number of distances required
      len <- 0
      for (l in 1:(n_spots-1)){
        len <- len + l
      }
      
      # The true distribution
      eukl <- vector(mode="numeric", length=len)
      eukli <- matrix(0, ncol=n_spots, nrow=n_spots)
      
      # The conditional probabilities for the gene in question
      P_g <- sum(values[[sample]])/tot[[sample]]
      P_g_s <- vector(mode='numeric', length=length(values[[sample]]))
      P_cond <- vector(mode='numeric', length=length(values[[sample]]))
      P_s <- vector(mode='numeric', length=length(values[[sample]]))
      for (l in 1:length(values[[sample]])){
        P_g_s[l] <- as.numeric(values[[sample]][l])/tot_spot[[sample]][l]
        P_s[l] <- P[[sample]][l]
      }
      P_cond <- P_g_s*P_s/P_g
      
      
      n_expressed <- as.numeric(non_zero[[sample]][all_genes[i],1])
      n_oversampling <- n_spots - n_expressed
      
      expressed_indices <- which(values[[sample]]!=0)
      oversampling_indices <- sample(indices[[sample]], size=n_oversampling, 
                                     replace=TRUE, prob=P_cond)
      gene_indices <- sort(c(expressed_indices, oversampling_indices))
      eukli <- euk[[sample]][gene_indices, gene_indices]
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
        interm_eukl <- euk[[sample]][random_indices, random_indices]
        eukl_for_rand <- interm_eukl[lower.tri(interm_eukl, diag=FALSE)]
        counts_matrix_rand[[sample]][k,] <- hist(eukl_for_rand, breaks, 
                                                 plot=FALSE)$counts
        
      }
      
    }
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
  
  # Plot the expression of each gene on the array, as well as
  #  the difference between the gene's distribution and the mean 
  #  distribution, along with the sample number:
  
  gg[[i]] <- list()
  for (sample in 1:length(testdata)){
    if (identical(which(genes[[sample]]==all_genes[i]), integer(0))){
      df <- data.frame()
      gg1 <- ggplot(df) + 
        geom_point() + 
        xlim(0, 10) + ylim(0, 100)+
        theme_classic()+
        ggtitle(paste(all_genes[i], sample), "\n Sample not included \n in analysis")+
        theme(plot.title = element_text(hjust = 0.5))+
        theme(legend.position="none", axis.line=element_blank(),
              axis.text.x=element_blank(), axis.text.y=element_blank(),
              axis.ticks=element_blank(), axis.title.x=element_blank(),
              axis.title.y=element_blank())
      gg[[i]][[sample]] <- gg1
    
    } else{
      xcoords <- as.numeric(sapply(strsplit(colnames(data[[sample]]), 
                                            "x"), "[[", 1))
      ycoords <- as.numeric(sapply(strsplit(colnames(data[[sample]]), 
                                            "x"), "[[", 2))
      df <- as.data.frame(cbind(xcoords, ycoords))
      df$tot_counts <- as.numeric(as.vector(data[[sample]][all_genes[i],]))
      gg1 <- ggplot(df, aes(x=xcoords, y=ycoords))+
        geom_point(data=df, aes(x=xcoords, y=ycoords, color=tot_counts), size=4)+
        scale_color_gradient(low="yellow", high="red")+
        coord_fixed()+
        theme_classic()+
        ggtitle(paste(all_genes[i], sample))+
        theme(plot.title = element_text(hjust = 0.5))+
        theme(legend.position="none", axis.line=element_blank(),
              axis.text.x=element_blank(), axis.text.y=element_blank(),
              axis.ticks=element_blank(), axis.title.x=element_blank(),
              axis.title.y=element_blank())
      gg[[i]][[sample]] <- gg1
    }
  }
  
  df <- as.data.frame(diff_real)
  df$distance <- breaks[2:(length(breaks))]-0.5
  gg2 <- ggplot(df, aes(x=bins, y=diff_real))+
    geom_bar(stat="identity")+
    theme_bw()+
    theme(axis.title.y=element_blank())
  gg[[i]][[length(testdata)+1]]
  
  cat(c(i, " out of ", length(all_genes), " finished\n"))
}


# In the paper by Rodriques et al. p-values of less than 0.005 were
#  deemed significant.
diff_expr <- genes[which(p[,1]<0.005)]

# Alternatively: Calculating adjusted p values in order to correct for 
#  multiple hypothesis testing, using the Benjamini-Hochberg procedure.
#  N.B. - this requires that the tests are independent or positively
#  dependent, not sure if this can really be said to be the case?
#  Otherwise, the Bejamini-Yekutieli procedure can be used.
# p-values saved as 1/(10n) are first converted to 1/n, as it is 
#  better to be slightly too conservative in this case.
# Accepting a FDR of 5%, all genes with an adjusted p value of less 
#  than 0.05 are said to be spatially variable.

for (val in which(p[,1]==0.0001)){
  p[val,1] <- 0.001
}

p_adj <- p.adjust(p[,1], method="BH")
p_adj <- as.data.frame(p_adj, row.names=rownames(p))
colnames(p_adj) <- "p_adj"
p_adj$no_samples <- p[,2]

# As output, the adjusted p values for each gene and the number of
#  samples used in the calculation should be given, sorted in order
#  of ascending adjusted p value.
p_adj <- p_adj[order(p_adj$p_adj),]
print(p_adj)


# diff <- which(p_adj[,1]<0.05)
# p_adj_diff <- as.data.frame(p_adj[diff,], row.names=rownames(p_adj)[diff])
# diff_expr_MHT <- rownames(p_adj_diff)

print(Sys.time() - start_time)


# Useful for later
# plot(x=xcoords, y=ycoords, col=alpha(color_vector, 1), lwd=1, asp=1,
#      ylab="", xlab="", main=paste(gene, round(p_adj_diff[gene,1], digits=4)), pch=19, cex.main=1,
#      xaxt="n", yaxt="n", bty="n", col.main="black")