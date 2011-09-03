Flim <- function(singleton.counts,
                 pairwise.counts,
                 document.count,
                 beta.1 = 0.1,
                 beta.2 = 0.1,
                 num.iterations = 15) {
  N <- length(singleton.counts)

  if (max(pairwise.counts[,1]) > N ||
      min(pairwise.counts[,1]) < 0 ||
      max(pairwise.counts[,2]) > N ||
      min(pairwise.counts[,2]) < 0) {
    stop("Atrocity!  Pairwise count indices must be between 1 and N.");
  }

  flim.obj <- new(.module$Flim, N, beta.1, beta.2)
  

  flim.obj$loadCorpus(singleton.counts, 
                      pairwise.counts[,1],
                      pairwise.counts[,2],
                      pairwise.counts[,3],
                      document.count)
  for (ii in 1:num.iterations) {
    cat("Iteration ")
    print(ii)
    print(system.time({
      flim.obj$estimateExpectations()
      cat("Total change: ")
      print(flim.obj$optimizeAll())
    }))
  }
  return(flim.obj)
}
