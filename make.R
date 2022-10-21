args = commandArgs(trailingOnly = TRUE)

cat("*** Compiling ... ***\n")
Rcpp::compileAttributes()

cat("*** Loading ... ***\n")
devtools::document()

if (length(args) > 0 && nchar(args) >= 2 && substring(args, 1, 1) == "-") {
  if (grepl("t", args)) {
    cat("*** Testing ... ***\n")
    devtools::test()
  }
  
  if (grepl("p", args)) {
    cat("*** Pkgdown ... ***\n")
    pkgdown::build_site()
  }
}
