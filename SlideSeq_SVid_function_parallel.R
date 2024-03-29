#!/usr/bin/Rscript
#-------------------------------------------------------------------------
# Moira Ek
# 2019 08 13
#
# Adaptation of the SlideSeq method to ST data, by taking the
#  expression level into account as described in 
#  SlideSeq_function_genes_as_input.R. Here, a function taking
#  a vector of paths to multiple count matrices as input, and giving a
#  list of spatially variable genes as output is described.
#
# For each gene, the distances (both of the real distribution and
#  within the 1000 random samples) are calculated separately for
#  each sample, and these are then combined into one distance
#  distribution, from which the method continues as if working with
#  one sample. New in this version: try if taking the mean of the
#  distributions from the different samples gives fewer SV genes. If
#  the large number when using multiple samples is due to the increased
#  no. distances this should be an improvement. Also new: 
#  parallelization of the code.
#
# Input:
#  Samples = a vector of paths to the count matrices for all samples that
#   are to be included in the analysis.
#  Filter_and_normalize = TRUE or FALSE, if TRUE, Samples will be 
#   filtered to remove spots with expression of less than 200 genes, 
#   and genes with expression in less than 10 spots, followed by 
#   normalization by SCTransform (Seurat package). If FALSE, it is
#   assumed that Samples has already been filtered and normalized, and
#   the count matrices found in Samples are used as they are.
#  Genes_as_rows = TRUE or FALSE. If TRUE, the count matrices specified
#   by Samples are assumed to have genes as row names and spot 
#   coordinates as column names. If FALSE, the reverse is assumed, and
#   the matrices will be transposed.
#  FDR_limit = The accepted false discovery rate, adjusted p values
#   smaller than this value will be regarded as significant.
#  NoCores = If a number, the calculations will be divided over that
#    number of cores. If FALSE, the number of available cores will be
#    obtained and used.

# Output:
#  p_adj = A data frame consisting of the names and (Benjamini-
#   Hochberg-) adjusted p values for the genes that were identified as
#   spatially variable, in order of ascending p-value. The data frame
#   also includes a column for the number of samples that were included
#   in the analysis leading to the identification of this particular 
#   gene as spatially variable.
#
#  
# N.B. - Fix so that a single gene also works
# Problem: in order to not have to save all of the distributions for 
#  each random sample, sample, etc, the genes can't be grouped. Thus,
#  each gene must be treated individually. Try to decrease the impact
#  of this on the running time by parallelizing the code.
# N.B. - This gives a large number of "spatially variable" genes. 
#  A similar problem is descibed in the SlideSeq paper - look at how
#  they solved it.
#
#-----------------------------------------------------------------------

library(parallel)
library(Seurat)

SlideSeq_SVid <- function(Samples, Filter_and_normalize, 
                          Genes_as_rows, FDR_limit, NoCores){
  
  #---------------------------------------------------------------------
  #  Data organization
  #---------------------------------------------------------------------
  
  # Create a list of the count matrices, one for each element in Samples
  data <- list()
  if (Genes_as_rows == TRUE){
    data <- lapply(Samples, function(x) as.data.frame(read.table(x, check.names=FALSE)))
  } else if (Genes_as_rows == FALSE) {
    data <- lapply(Samples, function(x) as.data.frame(t(read.table(x, check.names=FALSE))))
  } else {
    stop("Genes_as_rows must be either TRUE or FALSE")
  }
  
  # Filtering and normalization, if desired
  if (Filter_and_normalize == TRUE){
    # Remove genes with expression in less than 10 spots, and spots with
    #  expression of less than 200 genes.
    testdata <- lapply(data, function(x) x[(rowSums(x!=0)>9),])
    testdata <- lapply(testdata, function(x) x[,colSums(x!=0)>199])
    
    # Normalize the counts using SCTransform
    se <- lapply(testdata, CreateSeuratObject)
    se <- lapply(se, SCTransform)
    
    testdata <- lapply(se, function(x) GetAssayData(x, slot="counts", assay="SCT"))
    testdata <- lapply(testdata, as.matrix)
  } else if (Filter_and_normalize == FALSE){
    testdata <- data
  } else {
    stop("Filter_and_normalize must be either TRUE or FALSE")
  }
  
  #--------------------------------------------------------------------
  # Distance calculations
  #--------------------------------------------------------------------
  
  # Create a list of distance matrices, one for each sample
  
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
  
  # Calculate the maximum distance between included spots, across all
  #  samples, in order to define appropriate breakpoints for hist(). 
  #  1 due to the actual distances in the ST array.
  maxdist <- sapply(euk, max)
  maxdist <- max(maxdist)
  breaks <- seq(0,maxdist+1, by=1)
  
  #--------------------------------------------------------------------
  # Sampling probabilities
  #--------------------------------------------------------------------
  
  # Calcuate the probability of sampling the spots as part of the 
  #  background. (Proportional to the total number of transcripts at 
  #  each spot.) Indices are calculated since the actual output from the
  #  sampling will be the index a specific spot has in the distance 
  #  matrix. 
  # This is done for each sample, creating two lists, of probabilities
  #  and indices, respectively.
  
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
  
  # The number of random samples to draw for the background 
  #  distribution.
  n <- 1000
  
  # All genes included in the analysis for each sample
  genes <- list()
  for (sample in 1:length(testdata)){
    genes[[sample]] <- rownames(testdata[[sample]])
  }
  
  # All genes for which a p value can be calculated (i.e. the genes 
  #  that are found in at least one sample after filtering) 
  all_genes <- unique(unlist(genes))
  
  # # A dataframe for the p-values.
  # p <- matrix(0, nrow=length(all_genes), ncol=2)
  # p <- as.data.frame(p, row.names=all_genes)
  # colnames(p) <- c("p", "no_samples")
  
  # A list of, for each sample, the number of spots with non-zero 
  #  expression of each gene (non_zero). non_zero_with_oversampl is
  #  a list of, for each sample, the total number of spots that will be
  #  sampled.
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
  # Go through the genes one by one. For each gene, each sample is
  #  first treated individually: the true distance distribution for the 
  #  gene in the sample in question is obtained by including all spots 
  #  with expression of the gene, and an additional number of spots
  #  (as decided by factor_seur), sampled based on the conditional 
  #  probabilities for the gene in question. The distances between all
  #  of these spots are extracted from the distance matrix for the 
  #  sample, and the distances are binned. For each sample, a background 
  #  is then obtained, by first sampling the same number of spots as in 
  #  which the gene in question is expressed, without replacement and 
  #  based on the total number of transcripts at each spot. An additional
  #  number of spots (the same as for the true distribution) is then 
  #  sampled with replacement. The distance distribution is saved after 
  #  binning, and n-1 more such randomdistributions are sampled, forming
  #  the background distribution.
  #  Then, the distance distributions for each sample (still looking at 
  #  one gene) are summed (the true distributions, and each of the n
  #  background distributions). All samples are now treated as one, and
  #  the mean background distribution is obtained. The difference between
  #  the true distribution and this mean distribution is calculated, as
  #  well as between each random distribution and the mean distribution.
  #  The L1 norms of these distances are obtained (by taking the sum of
  #  the absolute values), and the p-value is calculated as 
  #  p=(no random samples with L1>=L1(true sample))/(total no random 
  #  samples). If the numerator is equal to 0, it can however only
  #  be said that p<1/n. These cases are currently saved as p=1/(10n).
  #  In the supplementary material to the paper by Rodriques et al. 
  #  (Science 363, 1463–1467, 2019. doi 10.1126/science.aaw1219)
  #  it is stated that p=(no random samples with L1<L1(true sample))/
  #  (total no random samples), but this contradicts their 
  #  argumentation, as well as the example in figure S10.
  #--------------------------------------------------------------------
  
  if (NoCores==FALSE){
    # Identify the number of available cores
    ncores <- detectCores()
  } else if (is.numeric(NoCores)==TRUE){
    ncores <- NoCores
  } else{
    stop("NoCores must be either a number of cores, or equal to FALSE")
  }
  
  # The samples need to be treated together, but all genes don't have to
  #  be analyzed in the same core.
  # Divide the genes, ceiling as to not miss any.
  ngenes <- ceiling(length(all_genes)/ncores)
  
  # Group the genes
  genes_core <- list()
  for (k in 1:ncores){
    genes_core[[k]] <- all_genes[(1+(k-1)*ngenes):(k*ngenes)]
    genes_core[[k]] <- na.omit(unique(genes_core[[k]]))
  }
  
  
  until_unadj_p <- function(genes_to_analyze) {
    p <- matrix(0, nrow=length(genes_to_analyze), ncol=2)
    p <- as.data.frame(p, row.names=genes_to_analyze)
    colnames(p) <- c("p", "no_samples")
    
    #Go through the genes one by one
    for (i in 1:length(genes_to_analyze)){
      
      #Prepare lists
      eukdistr <- list()
      counts_matrix_rand <- list()
      
      # For each gene, keep count of the number of samples that were
      #  included in the analysis.
      included_samples <- length(testdata)
      for (sample in 1:length(testdata)){
        testdat <- testdata[[sample]]
        values <- vector(mode="numeric", length=ncol(testdat))
        if (identical(which(rownames(testdat)==genes_to_analyze[i]), integer(0))){
          # If the gene is not present in this sample, no distances should
          #  be generated, and the values vector remains a vector of zeroes
          # One sample not included in the analysis for this gene
          included_samples <- included_samples - 1
        } else {
          # The counts for the gene in question, if it is included
          values <- testdat[genes_to_analyze[i],]
          
          # The number of spots required.
          n_spots <- as.numeric(non_zero_with_oversampl[[sample]][genes_to_analyze[i],1])
          
          # The number of distances required
          len <- 0
          for (l in 1:(n_spots-1)){
            len <- len + l
          }
          
          # The true distribution
          eukl <- vector(mode="numeric", length=len)
          eukli <- matrix(0, ncol=n_spots, nrow=n_spots)
          
          # The conditional probabilities for the gene in question
          P_g <- sum(values)/tot[[sample]]
          P_g_s <- vector(mode='numeric', length=length(values))
          P_cond <- vector(mode='numeric', length=length(values))
          P_s <- vector(mode='numeric', length=length(values))
          for (l in 1:length(values)){
            P_g_s[l] <- as.numeric(values[l])/tot_spot[[sample]][l]
            P_s[l] <- P[[sample]][l]
          }
          P_cond <- P_g_s*P_s/P_g
          
          # The number of spots with expression of the gene in question
          n_expressed <- as.numeric(non_zero[[sample]][genes_to_analyze[i],1])
          
          # The number of additional spots to sample
          n_oversampling <- n_spots - n_expressed
          
          # Saving the indices of the spots with expression, as well as
          #  of the additional sampled spots. The distances are then
          #  extracted from the correct distance matrix
          expressed_indices <- which(values!=0)
          oversampling_indices <- sample(indices[[sample]], size=n_oversampling, 
                                         replace=TRUE, prob=P_cond)
          gene_indices <- sort(c(expressed_indices, oversampling_indices))
          eukli <- euk[[sample]][gene_indices,gene_indices]
          eukl <- eukli[lower.tri(eukli, diag=FALSE)]
          
          # Binning of these distances
          eukdistr[[sample]] <- hist(eukl, breaks, plot=FALSE)$counts
          
          
          # The random distributions for each sample, saved in a list
          eukl_for_rand <- vector(mode="numeric", length=len)
          
          counts_matrix_rand[[sample]] <- matrix(0, nrow = n, 
                                                 ncol = (length(breaks)-1))
          
          # The same number of indices as spots with expression is first
          #  sampled without replacement, then an additional number is
          #  sampled with replacement. Distances are extracted and binned,
          #  and this is repeated n times in total, to obtain the background
          #  distributions.
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
          
        }
        
        # The number of samples included in the analysis for the gene 
        #  in question is saved
        p[genes_to_analyze[i],2]<- included_samples
      }
      
      
      #Combining the distributions for each sample, both true and random:
      # True:
      eukdistrib <- colSums(do.call(rbind, eukdistr))
      
      # Random:
      rand_distrib <- matrix(0, nrow = n, ncol = (length(breaks)-1))
      for (m in 1:n){
        # Extracting the m:th row from each matrix in the list, taking the
        #  mean of the number of elements in each bin
        rand_distrib[m,] <- colSums(do.call(rbind, lapply(counts_matrix_rand, `[`,m,)))
      }
      
      # The mean background distribution
      mean_counts <- colSums(rand_distrib)
      
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
      
      # The number of L1 norms of random samples that are greater than 
      #  (or equal to, as to not create any strange effects, e.g.  in 
      #  the extreme case where all L1 norms are 0) the
      #  L1 norm for the true distribution (minus the mean background)
      over_L1 <- L1_norms_rand[which(L1_norms_rand>=L1_norm_real)]
      
      # Calculation of p value
      if (length(over_L1)!=0){
        p[genes_to_analyze[i],1] <- length(over_L1)/n
      } else {
        p[genes_to_analyze[i],1] <- 1/(10*n)
      }
      
      # Progress tracking
      #cat(c(i, " out of ", length(all_genes), " finished\n"))
    }
    return(p)
  }
  
  RNGkind("L'Ecuyer-CMRG")
  p_unadj <- mclapply(genes_core, until_unadj_p, mc.cores=ncores)
  
  p <- do.call("rbind", p_unadj)
  
  #----------------------------------------------------------------------
  # p values
  #----------------------------------------------------------------------
  # In the paper by Rodriques et al. p values of less than 0.005 were
  #  deemed significant.
  #diff_expr <- genes[which(p[,1]<0.005)]
  
  # Alternatively: Calculating adjusted p values in order to correct for 
  #  multiple hypothesis testing, using the Benjamini-Hochberg procedure.
  #  N.B. - this requires that the tests are independent or positively
  #  dependent, not sure if this can really be said to be the case?
  #  Otherwise, the Bejamini-Yekutieli procedure can be used.
  # p-values saved as 1/(10n) are first converted to 1/n, as it is 
  #  better to be slightly too conservative in this case.
  
  for (val in which(p[,1]==0.0001)){
    p[val,1] <- 0.001
  }
  
  # Benjamini-Hochberg MHT correction
  p_adj <- p.adjust(p[,1], method="BH")
  p_adj <- as.data.frame(p_adj, row.names=rownames(p))
  colnames(p_adj) <- "p_adj"
  # Including the number of samples included in the analysis for each 
  #  gene in this data frame
  p_adj$no_samples <- p[,2]
  
  # Sorting in order of ascending adjusted p value
  p_adj <- p_adj[order(p_adj$p_adj),]
  
  diff <- which(p_adj[,1]<FDR_limit)
  p_adj_diff <- as.data.frame(p_adj[diff,], row.names=rownames(p_adj)[diff])
  #diff_expr_MHT <- rownames(p_adj_diff)
  
  return(p_adj_diff)
  
}
