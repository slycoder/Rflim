NAMESPACE <- environment()
.module <- new("Module")

.onLoad <- function(libname, pkgname) {
  require(methods)
  unlockBinding(".module", NAMESPACE)
  assign(".module", Rcpp:::Module("Rflim"), NAMESPACE)
  lockBinding(".module", NAMESPACE)
}

