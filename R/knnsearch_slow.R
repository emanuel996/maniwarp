knnsearch_slow <- function(X, k = 10){
  # slower than cpp
  aa = rowSums(X^2)
  n = length(aa)
  D = matrix(aa, n, n) + matrix(aa, n, n, byrow = TRUE) - 2 * tcrossprod(X)
  index = matrix(0, n, k)
  for (i in 1:n){
    output = sort(D[i, ], index.return = TRUE)
    index[i, ] = output$ix[2 : (k+1)]
  }
  G = matrix(0, n, n)
  for (i in 1:n){
    for (j in 1:k){
      if (! is.na(D[i, index[i,j]])){
        G[i, index[i, j]] = D[i, index[i, j]]
      }
    }
  }
  # for (i in 1:n){
  #   G[i, index[i, ]] = D[i, index[i, ]]
  # }
  return(list(index = index, matrix = G))
}



