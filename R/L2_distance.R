L2_distance <- function(a, b, df = 0){
  # a and b are k*n matrices, representing n k-dimensional points
  # D(i,j) is the L2 distance between a(i) and b(j)
  k = dim(a)[1] 
  n = dim(a)[2] 
  
  aa = colSums(a^2)
  bb = colSums(b^2)
  ab = t(a) %*% b
  D = matrix(aa, n, n) + matrix(bb, n, n, byrow = TRUE) - 2 * ab
  D = sqrt(D)
  return(D)
}