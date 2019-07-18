#!/bin/Rscript
#-----------------------------------------------------------------------
#
# Moira Ek
# 2019 07 18
# Försök att implementera SlideSeq-DGE-metoden
# 1. Direkt, utan att ännu ha gjort några anpassningar till ST,
#    jobbar med råa countdata initialt. Använder spotkoordinater för
#    avstånd.
# 2. Försök att optimera koden för att försnabba processen.
# 3. Test: går det att inte behöva sampla på nytt för varje gen?
#     -testa att gruppera ihop de gener som har samma antal spots där
#     de uttrycks.
# 
#-----------------------------------------------------------------------

setwd("/home/moiraek/summerp19/SlideSeq_etc/Till_git")

data <- as.data.frame(t(read.table("Rep1_MOB_count_matrix-1.tsv", check.names=FALSE)))

#Ta bort gener med yttryck i färre än 5 spots
testdata <- data[rowSums(data!=0)>4,]

#--------------------------------------------------------------------
# Skapa en avståndsmatris
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
# Beräkna sannolikheter för framslupning av spots. Indices då den
#  faktiska framslumpningen är av det index en viss spot har i
#  avståndsmatrisen.
#--------------------------------------------------------------------

tot_spot <- colSums(testdata)
tot <- sum(tot_spot)
P <- vector(mode="numeric", length=ncol(testdata))
for (val in 1:length(P)){
  P[val] <- tot_spot[val]/tot
}
indices <- 1:length(P)


#--------------------------------------------------------------------
# Gruppering av generna efter antalet spots de uttrycks i --> antalet
#  spots som ska samplas.
#--------------------------------------------------------------------

# Beräkna det maximala avståndet mellan inkluderade spots, för att
#  sätta upp lämpliga breakpoints för hist. 0.5 p.g.a. de faktiska
#  avstånden i ST-arrayen.
maxdist <- max(euk)
breaks <- seq(0,maxdist+0.5, by=0.5)


# n = antal slumpfördelningar att dra för varje gen.
# p = dataframe för p-värdena
# non_zero = antalet spots med nollskilt uttryck för varje gen
n <- 1000
p <- vector(mode="numeric", length=nrow(testdata))
p <- as.data.frame(p, row.names=rownames(testdata), 
                   col.names="p-value")

non_zero <- matrix(nrow=nrow(testdata),ncol=1)
rownames(non_zero)<-rownames(testdata)
for (i in 1:nrow(non_zero)){
  non_zero[i,1] <- ncol(testdata[i, which(testdata[i,]!=0)])
}

# Grupperar generna som uttrycks i samma antal spots på raderna i
#  equal_no_spots, ger varje rad antalet spots som namn. unique -->
#  varje grupp av gener finns bara med en gång.
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
# Gå igenom grupperna en och en. För varje grupp samplas först
#  lika många index som spots där dessa gener uttrycks,
#  countsen i varje bin definierad enligt ovan sparas på raderna i 
#  counts_matrix_rand, och detta upprepas 1000 (n) gånger, tills
#  matrisen är fylld. Kolumnvisa medelvärden av denna matris ger 
#  medelfördelningen, mean_counts. diff_matr_rand sparar undan
#  differenserna mellan varje enskild slumpfördelning och medel-
#  fördelningen. L1-normerna av detta sparas i L1_norms_rand.
#  Sedan gås generna igenom en och en. Tar först ut indexen för
#  de nollskilda elementen, d.v.s. de spots där genen i fråga 
#  uttrycks. Detta används för att ta fram avstånden, som sparas i
#  eukl. Differensen mellan fördelningen av dessa avstånd (i samma
#  bins som för slumpfördelningarna) och medelfördelningen tas,
#  L1-normen beräknas (L1_norm_real), och p-värdet beräknas, som 
#  p=(antal slumpmässiga prover med L1 > L1(sanna provet))/
#  (totalt antal slumpmässiga prover). Om täljaren är 0 kan dock 
#  bara sägas att p<1/n. Dessa fall sparas i nuläget som p=1/(10n),
#  för hanteringens skull. I suppl. till artikeln säger de att 
#  p=(antal slumpm prover med L1<L1(sanna))/(totalt antal slumpmässiga),
#  men detta motsäger resonemanget de för, och exemplet i figur S10. 
#--------------------------------------------------------------------

for (i in 1:nrow(equal_no_spots)){
  # Countsen för generna i fråga
  values <- testdata[na.omit(equal_no_spots[i,]),]
  
  # Antalet avstånd som behövs
  len <- 0
  for (l in 1:(as.numeric(rownames(equal_no_spots)[i])-1)){
    len <- len + l
  }
  
  # Slumpfördelningarna
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
  
  # Elementvis differens mellan avståndsfördelningarna för
  #  slumpsamplen, och medelfördelningen.
  diff_matr_rand <- matrix(0, nrow=n, ncol=(length(breaks)-1))
  diff_matr_rand <- t(apply(counts_matrix_rand,1,'-',mean_counts))
  
  # Ta beloppet av dessa, och summera över varje rad --> L1-normen
  #  för avståndet mellan de slumpade fördelningarna och
  #  den genomsnittliga fördelningen. 
  abs_diff_rand <- abs(diff_matr_rand)
  L1_norms_rand <- rowSums(abs_diff_rand)
  
  
  # De verkliga fördelningarna
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
    
    # p-värdesberäkning
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
# Snabbplottning, test
# --------------------------------------------------------------------
library(ggplot2)

gene <- "2010300C02Rik"
col <- as.numeric(as.vector(testdata[which(rownames(testdata)==gene),]))
# Skapa färggradient
rbPal <- colorRampPalette(c('yellow','red'))
color_vector <- rbPal(10)[as.numeric(cut(col,breaks = 10))]

xcoords <- as.numeric(sapply(strsplit(colnames(data), "x"), "[[", 1))
ycoords <- as.numeric(sapply(strsplit(colnames(data), "x"), "[[", 2))

# Plotta arrayen med färg enligt gradienten
plot(x=xcoords, y=ycoords, col=alpha(color_vector, 1), lwd=1, asp=1,
     ylab="", xlab="", main=paste(gene), pch=19, cex.main=1.5, 
     xaxt="n", yaxt="n", bty="n", col.main="black")
