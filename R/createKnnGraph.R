createKnnGraph <- function(N1){
  # N1 is a k*n matrix which shows the k-nearest neighbors of n points
  library('Matrix')
  n1 = dim(N1)[1]
  W1 = Matrix(0, n1, n1)
  for (i in 1:n1){
    W1[i, N1[, i]] = 1
  }
}