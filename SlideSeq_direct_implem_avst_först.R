#-----------------------------------------------------------------------
#
# Moira Ek
# 2019 07 12
# Försök att implementera SlideSeq-DGE-metoden
# 1. Direkt, utan att ännu ha gjort några anpassningar till ST,
#    jobbar med råa countdata initialt. Använder spotkoordinater för
#    avstånd.
# Försöker förbättra hanteringen, möjliggöra för alla gener. 
# 
#-----------------------------------------------------------------------

setwd("/home/moiraek/summerp19/SlideSeq_etc/")
library(STutility)
library(Seurat)

data <- as.data.frame(t(read.table("Rep1_MOB_count_matrix-1.tsv", check.names=FALSE)))

testdata <- data[1:10,]
#--------------------------------------------------------------------
# Funktioner:
#--------------------------------------------------------------------


# Skapa funktion för att ge en slumpmässig uppsättning spots, med
#  SL enligt mängden transkript de fångat. Utan återläggning.
#  STdata är en dataframe med spots som kolumner, gener som rader 
#  (rådatan), nr är antalet spots som önskas. Nr kan aldrig vara
#  större än det totala antalet spots i den tänkta tillämpningen
#  --> behöver ej lägga in en if-sats.

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


# # Funktion för att ge en fördelning av de Euklidiska avstånden mellan
# #  en given uppsättning spots, spot_coords är en vektor av spotnamn 
# #  på formen "x-koordinat x y-koordinat". För varje spot beräknas 
# #  det euklidiska avståndet till alla spots till höger om den (--> 
# #  avståndet för varje par av spots tas med en gång). 
# #  Output: avstånden. 
# 
# distrib <- function(spot_coords){
#   len <- 0
#   for (val in 1:(length(spot_coords)-1)){
#     len <- len + val
#   }
#   euk <- vector(mode="numeric", length=len)
#   coords <- GetCoords(spot_coords, delim="x")
#   xcoords <- coords$x
#   ycoords <- coords$y
#   count <- 1
#   for (j in 1:(length(spot_coords)-1)){
#     rj <- c(xcoords[j], ycoords[j])
#     for (k in (j+1):length(spot_coords)){
#       rk <- c(xcoords[k], ycoords[k])
#       euk[count] <- sqrt(sum((rj-rk)^2)) 
#       count <- count+1
#     }
#   }
#   return(euk)
# }

#--------------------------------------------------------------------
# Skapa en avståndsmatris
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
# Försök
#--------------------------------------------------------------------

# Beräkna det maximala avståndet mellan inkluderade spots, för att
#  sätta upp lämpliga breakpoints för hist. Lät tidigare den verkliga
#  fördelningen sätta bingränserna (gav antal bins för den, tog sedan
#  $breaks som input till hist för de slumpade spotuppsättningarna), 
#  detta ger dock risken att alla avståndsvärden i de slumpade inte
#  täcks in, tex om en gen uttrycks i en begränsad region. (De
#  spotsen dras ju från hela arrayen, efter mängden fångade 
#  transkript.)
# NB - se mer på vad som är ett lämpligt binantal.
maxdist <- max(euk)
breaks <- seq(0,maxdist, by=(maxdist/200))


# Ta för varje gen ut spotsen med nollskilt värde ("values"), beräkna
#  avstånden. eukl innehåller avstånden mellan spotsen med uttryck
#  för varje gen, följt av nollor. Skapa först en matris för countsen
#  av (antal bins) kolumner, n+2 rader ty en fördelning för
#  de faktiska avstånden, n slumpmässiga spotfördelningar, och en
#  medel av dessa. n = antal slumpmässiga
#  spotuppsättningar som ska dras, slideseq använder 1000.
# p - för att spara undan p-värdet för varje gen.
n <- 1000
p <- vector(mode="numeric", length=nrow(testdata))

start_time <- Sys.time()
for (i in 1:nrow(testdata)){
  values <- testdata[i,which(testdata[i,]!=0)]
  spots <- colnames(values)
  
  len <- 0
  for (val in 1:(length(spots)-1)){
    len <- len + val
  }
 
  eukl <- vector(mode="numeric", length=len)
  count <- 1
  for (l in 1:(length(spots)-1)){
    for (m in (l+1):length(spots)){
      eukl[count]<-euk[spots[l],spots[m]]
      count <- count + 1
    }
  }
  
  eukdistr <- hist(eukl, breaks, plot=FALSE)
  counts_matrix <- matrix(0, nrow = (n + 2), ncol = (length(breaks)-1))
  counts_matrix[1,] <- eukdistr$counts
  
  # Samma antal avstånd för varje slumpmässig uppsättning spots för
  #  varje gen som för de faktiska spotsen där den uttrycks, n
  #  slumpmässiga avståndsfördelningar (i form av antalet counts i 
  #  varje bin) sparas på rad 2 till n+1 i count_matrix.
  # I eukl_for_rand sparas de faktiska avstånden. 
  
  eukl_for_rand <- matrix(0, nrow = n, ncol = length(eukl))
  
  for (k in 1:n){
    random_spots <- SpotProb(testdata,length(values))
    count <- 1
    for (l in 1:(length(random_spots)-1)){
      for (m in (l+1):length(random_spots)){
        eukl_for_rand[k, count]<-euk[spots[l],spots[m]]
        count <- count + 1
      }
    }
    counts_matrix[(k+1),] <- hist(eukl_for_rand[k,], breaks, 
                                  plot=FALSE)$counts
  }
  
  # Tolkar elementwise mean som att ta medel av countsen per bin  
  counts_matrix[(n+2),] <- colMeans(counts_matrix[(2:(n+1)),])
  
  # Elementvis differens för avståndsfördelningen för GOI och 
  #  medelfördelningen, samt slumpfördelningarna och medel.
  diff_matr <- matrix(0, (n+1), ncol(counts_matrix))
  diff_matr[1,] <- counts_matrix[1,]-counts_matrix[n+2,]
  for (m in 1:n){
    diff_matr[m+1,] <- counts_matrix[m+1,]-counts_matrix[n+2,]
  }
  
  # Ta beloppet av dessa, och summera över varje rad --> L1-normen
  #  för avståndet mellan den verkliga fördelningen (element 1) 
  #  respektive de slumpade fördelningarna (element 2-(n+1)) och
  #  den genomsnittliga slumpfördelningen. n/2 bins, tex
  abs_diff <- abs(diff_matr)
  L1_norms <- rowSums(abs_diff)
  #hist(L1_norms[2:(n+1)], n/2, xlim=c(min(L1_norms),max(L1_norms)))
  #abline(v=L1_norms[1], col="red")
  
  # p-värdet är som minst p<1/n, där n är antalet slumpmässiga 
  #  spotuppsättningar. Generellt är p=(slumpmässiga prover med
  #  L1>L1(sanna provet))/(totalt antal slumpmässiga prover). Om
  #  täljaren är 0 kan dock bara sägas att p<1/n. I suppl. till 
  #  artikeln säger de att p=(slumpm prover med L1<L1(sanna))/(totalt
  #  antal slumpmässiga), men detta motsäger resonemanget de för,
  #  och exemplet i figur S10. För p<1/n: sparar i nuläget värdet
  #  som p=1/(10n), för hanteringens skull. Kanske fixa detta på  
  #  något snyggare sätt?
  
  over_L1 <- L1_norms[which(L1_norms>L1_norms[1])]
  if (length(over_L1)!=0){
    p[i] <- length(over_L1)/(length(L1_norms)-1)
  } else {
    p[i] <- 1/(10*n)
  }
  
  
}


# Lägg till som metadata till assayen i seuratobjektet
p <- as.data.frame(p, row.names=rownames(testdata), col.names="p-value")
#se[["RNA"]] <- AddMetaData(se[["RNA"]], p)

#se@assays$RNA[[]]

diff_expr <- rownames(testdata[which(p<0.005),])


end_time <- Sys.time()


# --------------------------------------------------------------------
# Snabbplottning, test
# --------------------------------------------------------------------
library(ggplot2)

gene <- "Chrm1"
col <- as.numeric(as.vector(testdata[which(rownames(testdata)==gene),]))
# Skapa färggradient
rbPal <- colorRampPalette(c('yellow','red'))
color_vector <- rbPal(10)[as.numeric(cut(col,breaks = 10))]

coords <- GetCoords(colnames(data), delim="x")
xcoords <- coords$x
ycoords <- coords$y

# Plotta arrayen med färg enligt gradienten
plot(x=xcoords, y=ycoords, col=alpha(color_vector, 1), lwd=1, asp=1,
     ylab="", xlab="", main=paste(gene), pch=19, cex.main=1.5, 
     xaxt="n", yaxt="n", bty="n", col.main="black")
