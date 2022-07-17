# Complications arising with in-place operations between R and Rcpp. If you want
# to do something in place in Rcpp, then you can do it, but: 
# 1. Beware of implicit copying: int-to-float conversion for NumericVector makes 
#    a copy! So it won't be modified by reference.
# 2. R never lets you modify by reference in a usual function environment, so
#    you need to let .Call() handle it. Somehow with .Call(), or a single-line  
#    R function that calls a .Call() function, you can do it without extra
#    tricks. But once you have a multi-line function, you can't do it without
#    some extra tricks ... 
# Therefore the design (just set things equal to one line .Call() functions) in 
# dot_functions.R is purposeful ... 
