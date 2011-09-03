# This function takes an lda document structure and returns a sparse matrix
# giving counts of how many documents in which each pair of words co-occur.
count.pairs <- function(documents) {
  ## documents is uniquified per document already, so we can ignore the counts.
  w <- lapply(documents, function(x) x[1,])
  ## infer the size of the matrix from the  documents
  V <- max(unlist(w)) + 1L
  M <- Matrix(0, V, V)

  ## create a list giving the outer product indices for each document.
  w.pairs <- lapply(w, function(ww) {
    cbind(rep(ww + 1L, length(ww)), 
          rep(ww + 1L, each=length(ww)))
  })

  ## cross tabulate to get counts.
  M <- xtabs(~ X1 + X2,
             data.frame(do.call(rbind, w.pairs)), 
             sparse=T)
}
