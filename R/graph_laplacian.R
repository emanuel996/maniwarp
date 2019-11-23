graph_laplacian <- function(W){
  # W is an adjacency matrix
  # L = diag(sum(W)) - W
  library('Matrix')
  library('crqa')
  
  N = dim(W)[1]
  D = sqrt(colSums(W))
  D[D==0] = 1
  D = Diagonal(x = 1/D)
  L = Diagonal(N) - D %*% W %*% D
  # example 
  # W1=[   0    1    0    0
  #        1    0    0    0
  #        0    0    0    1
  #        0    0    1    0]
  # W1=[   1    -1    0    0
  #        -1    1    0    0
  #        0    0   1    -1
  #        0    0    -1    1]
  return(L)
}