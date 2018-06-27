library("TDA")

#Read in centroids of PA vtds
data <- read.table("vtds_pa.txt", sep = " ")

plot(data, xlab="", ylab="")

#compute filtration and homology
Diag <- alphaComplexDiag(X = data, maxdimension = 1,
                 library = c("GUDHI", "Dionysus"), printProgress = FALSE, location = TRUE)

#plot stuff
par(mfrow=c(1,3))
plot(data, xlab="", ylab="")
plot(Diag[["diagram"]])
plot(Diag[["diagram"]], barcode = TRUE)


DiagAlphaShape <- Diag
X <- data


#plot representative feature
par(mfrow = c(1, 2))
plot(DiagAlphaShape[["diagram"]])
plot(X[, 1:2], col = 2, main = "Representative loop of alpha shape filtration")
one <- which(DiagAlphaShape[["diagram"]][, 1] == 1)
one <- one[which.max(
  DiagAlphaShape[["diagram"]][one, 3] - DiagAlphaShape[["diagram"]][one, 2])]
 for (i in seq(along = one)) {
   for (j in seq_len(dim(DiagAlphaShape[["cycleLocation"]][[one[i]]])[1])) {
     lines(
       DiagAlphaShape[["cycleLocation"]][[one[i]]][j, , 1:2], pch = 19,
       cex = 1, col = i)
     }
   }
par(mfrow = c(1, 1))