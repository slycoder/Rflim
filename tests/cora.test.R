require(lda)
require(Rflim)
require(Matrix)

data(cora.documents)
data(cora.vocab)

counts <- count.pairs(cora.documents)
singleton.counts <- diag(counts)
counts <- as(counts, 'dgTMatrix')
pairwise.counts <- subset(data.frame(
  i = counts@i + 1L,
  j = counts@j + 1L,
  x = counts@x), x > 0 & i < j)

save(pairwise.counts, singleton.counts, file="counts.Rdata")

flim.instance <- Flim(singleton.counts,
                      pairwise.counts,
                      length(cora.documents))

lambda <- flim.instance$getLambda()
save(lambda, file="lambda.Rdata")
