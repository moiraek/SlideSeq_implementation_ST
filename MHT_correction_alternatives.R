#--------------------------------------------------------------------
# Holm-Bonferroni method, to correct for multiple hypothesis testing
# (Not e.g. Benjamini-Hochberg, since the tests cannot reasonably be 
# assumed to be independent of each other?).
#--------------------------------------------------------------------


p2 <- p[order(p$p), , drop = FALSE]
count <- nrow(p2)

alpha <- 0.05
test <- matrix(NA, nrow=nrow(p2),ncol=1)
rownames(test)<- rownames(p2)
for (i in 1:count){
  test[i,1] <- (p2[i,1]>(alpha/(count+1-i)))
}

k <- min(which(test[,1]==TRUE))
if (k==1){
  print("No null hypothesis can be rejected")
} else{
  diff_expr_MHT <- p2[1:(k-1),]
}

# --> Too stringent under these conditions, where the smallest 
#  observable p-value is p<0.001 (for n=1000).


#----------------------------------------------------------------------
# Implementation of the Bejamini-Hochberg procedure. But can the tests
#  really be said to be independent --> is this actually applicable?

# NB - forts läs 1.3 i Benjamini&Yekutieli
#----------------------------------------------------------------------

p2 <- p[order(p$p), , drop = FALSE]
count <- nrow(p2)

q <- 0.05
test <- matrix(NA, nrow=nrow(p2),ncol=1)
rownames(test)<- rownames(p2)

# For the cases of which only p<0.001 (saved as p=0.0001) can be said
#  the procedure is only meaningful if these cases aren't the only
#  ones that pass the test. Since it is however certain that they
#  have p-values of less than 0.001, the below will never lead to any
#  conclusions of null hypothesis rejection if this cannot certainly
#  be supported.
for (val in which(p2[,1]==0.0001)){
  p2[val,1]<-0.001
}

for (i in 1:count){
  test[i,1] <- (p2[i,1]<=((i/count)*q))
}

suppressWarnings(k <- max(which(test[,1]==TRUE)))
if (k==-Inf){
  print("No null hypotheses can be rejected with certainty")
} else{
  diff_expr_MHT <- rownames(p2)[1:k]
}

#----------------------------------------------------------------------
# Implementation of the Bejamini-Yekutieli procedure, under arbitrary
#  dependence.
#----------------------------------------------------------------------

p2 <- p[order(p$p), , drop = FALSE]
count <- nrow(p2)

q <- 0.05
test <- matrix(NA, nrow=nrow(p2),ncol=1)
rownames(test)<- rownames(p2)


# For the cases of which only p<0.001 (saved as p=0.0001) can be said
#  the procedure is only meaningful if these cases aren't the only
#  ones that pass the test. Since it is however certain that they
#  have p-values of less than 0.001, the below will never lead to any
#  conclusions of null hypothesis rejection if this cannot certainly
#  be supported.
for (val in which(p2[,1]==0.0001)){
  p2[val,1]<-0.001
}

c <- 0
for (j in 1:count){
  c <- c + 1/j
}
 
for (i in 1:count){
  test[i,1] <- (p2[i,1]<=((i/(count*c))*q))
}

suppressWarnings(k <- max(which(test[,1]==TRUE)))
if (k==-Inf){
  print("No null hypotheses can be rejected with certainty")
} else{
  diff_expr_MHT <- rownames(p2)[1:k]
}
