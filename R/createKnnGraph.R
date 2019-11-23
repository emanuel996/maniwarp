createKnnGraph <- function(N1){
  # N1 is a n*k matrix which shows the k-nearest neighbors of n points
  library('Matrix')
  n1 = dim(N1)[1]
  W1 = Matrix(0, n1, n1)
  for (i in 1:n1){
    W1[i, N1[i, ]] = 1
  }
  return(W1)
  # example
  # N1 = [2 1 4 3]^T   (graph: 1-2 3-4 and n=4,K=1) 
  # W1 = [0 1 0 0
  #       1 0 0 0
  #       0 0 0 1
  #       0 0 1 0]
}