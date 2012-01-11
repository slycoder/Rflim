require(lda)
require(Rflim)
require(Matrix)

data(cora.documents)
data(cora.vocab)

## documents <- head(cora.documents)
documents <- cora.documents

counts <- count.pairs(documents)
singleton.counts <- diag(counts)
counts <- as(counts, 'dgTMatrix')
pairwise.counts <- subset(data.frame(
  i = counts@i + 1L,
  j = counts@j + 1L,
  x = counts@x), x > 0 & i < j)

save(pairwise.counts, singleton.counts, file="counts.Rdata")

flim.instance <- flim(singleton.counts,
                      pairwise.counts,
                      length(documents),
                      0.0, 0.0, 15)

lambda <- flim.instance$getLambda()
save(lambda, file="lambda.Rdata")
