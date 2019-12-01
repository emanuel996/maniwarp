# manifold alignment for non-linear case

manifold_nonlinear <- function(X1, X2, W1, W2, W12, mu = 1, max_dim = 20, epsilon = 1e-8){
  library(Rcpp)
  #source('R/aliDif.R')
  #source('R/my_components.R')
  #source('R/createKnnGraph.R')
  #source('R/graph_laplacian.R')
  #source('R/L2_distance.R')
  #source('R/laplacian_eigen.R')
  #source('R/rowBdSlow.R')
  #sourceCpp('src/rowBd.cpp')
  #sourceCpp('src/knnsearch.cpp')
  # Feature-level Manifold Projections. Two domains.
  # X1: M1*P1 matrix, M1 examples in a P1 dimensional space.
  # X2: M2*P1 matrix
  # W1: M1*M1 matrix. weight matrix for each domain.
  # W2: M2*M2 matrix
  # W12: M1*M2 sparse matrix modeling the correspondence of X1 and X2.
  # mu: used to balance matching corresponding pairs and preserving manifold topology.
  # max_dim: max dimensionality of the new space. (default: 200)
  # epsilon: precision. (default: 1e-8)
  # get sizes for convenience later
  M1 = size(X1)[1]
  M2 = size(X2)[1]
  P1 = size(X1)[2]
  P2 = size(X2)[2]
  # Create weight matrix
  mu = mu * (sum(W1) + sum(W2))/(2 * sum(W12))
  W = rbind(cbind(W1, mu * W12), cbind(mu * t(W12), W2))
  L = graph_laplacian(W)
  rm(W1, W2, W12)
  # Eigen decomposition
  output = eigen(L, only.values = FALSE)
  vecs = output$vectors
  vals = output$values
  output2 = sort(vals, decreasing = FALSE,index.return = TRUE)
  vals = output2$x
  index = output2$ix[1 : min(max_dim * 2, dim(L)[1])]
  vecs = t( t(vecs)/sqrt(colSums(vecs^2)) )
  # filter out eigenvalues that are ~= 0
  for (i in 1:length(vals)){
    if (vals[i] > epsilon){
      break
    }
  }
  # Compute mappings
  start = i
  m = length(index) - start +1
  #
  g1 = vecs[1:M1, index[start:(start+m-1)] ]
  g2 = vecs[(M1+1):(M1+M2), index[start:(start+m-1)] ]
  return(list(map1 = g1, map2 = g2))
}
