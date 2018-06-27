library("TDA")

doTDA <-function(fname, sname){
  #Read in centroids of vtds
  X <- read.table(fname, sep = " ")
  
  #compute filtration and homology
  DiagAlphaShape <- alphaComplexDiag(X = X, maxdimension = 1,
                                     library = c("GUDHI", "Dionysus"), printProgress = TRUE, location = TRUE)
  
  #plot stuff
  png(file=sprintf("vtd_topology_%s.png", state),width=1200,height=400)
  par(mfrow=c(1,3))
  par(oma=c(0,0,2,0))
  
  plot(X[, 1:2],  main = "Representative Loops", col = "black", cex = 0.5, xlab = "Lat", ylab = "Lon")
  H1 <- which(DiagAlphaShape[["diagram"]][, 1] == 1)
  lifespans = DiagAlphaShape[["diagram"]][H1, 3] - DiagAlphaShape[["diagram"]][H1, 2]
  
  n <- length(unique(lifespans))
  
  first <- H1[which.max(lifespans)]
  second <- H1[which(lifespans == sort(unique(lifespans),partial=n-1)[n-1])]
  third <- H1[which(lifespans == sort(unique(lifespans),partial=n-2)[n-2])]
  
  plotGenerator <- function(gen, c){
    for (i in seq(along = gen)){
      for (j in seq_len(dim(DiagAlphaShape[["cycleLocation"]][[gen[i]]])[1])) {
        lines(
          DiagAlphaShape[["cycleLocation"]][[gen[i]]][j, , 1:2], pch = 19,
          cex = 1, col = c)
      }
    }
  }
  
  plotGenerator(first, "red")
  plotGenerator(second, "green")
  plotGenerator(third, "blue")
  
  plot(DiagAlphaShape[["diagram"]], main = "Persistance")
  plot(DiagAlphaShape[["diagram"]], barcode = TRUE, main = "Barcode")
  
  title("PA VTD Topology", outer = TRUE)
  dev.off()
  
  #save objects
  saveRDS(DiagAlphaShape, file = sprintf("%sDiag", state), ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
}