laplacian_eigen <- function(X, no_dims = 2, k = 10, sigma = 1){
  # LAPLACIAN_EIGEN Performs non-linear dimensionality reduction using Laplacian Eigenmaps
  #
  # [mappedX, mapping] = laplacian_eigen(X, no_dims, k, sigma, eig_impl)
  #
  # Performs non-linear dimensionality reduction using Laplacian Eigenmaps.
  # The data is in matrix X, in which the rows are the observations(n) and the
  # columns the dimensions(d). The variable dim indicates the preferred amount
  # of dimensions to retain (default = 2). The variable k is the number of
  # neighbours in the graph (default = 12).
  # The reduced data is returned in the matrix mappedX.
  #
  source('R/L2_distance.R')
  source('R/my_components.R')
  source('R/knnsearch_slow.R')
  source('R/createKnnGraph.R')
  library('FNN')
  library('geigen')
  
  # construct neighborhood graph
  #G = get.knn(X, k)
  #G = G^2
  #G = G/pmax(G)
  G = knnsearch_slow(X, k)$matrix
  G = as.matrix(G)
  G = G/max(G)
  # Only embed largest connected component of the neighborhood graph
  blocks = my_components(G)
  count = matrix(0, 1, max(blocks))
  for (i in 1:max(blocks)){
    count[i] = length(which(blocks == i))
  }
  block_no = which.max(count)
  conn_comp = which(blocks == block_no)
  G = G[conn_comp, conn_comp]
  
  # Compute weights (W = G)
  # Compute Gaussian kernel (heat kernel-based weights)
  G[G != 0] = exp(-G[G != 0] / (2 * sigma ^ 2))
  
  # Construct diagonal weight matrix
  D = diag(rowSums(G));
  
  # Compute Laplacian
  L = D - G
  L = matrix(L, ncol = dim(L)[1])
  D = matrix(D, ncol = dim(D)[1])
  L[is.na(L)] = 0
  G[is.na(G)] = 0
  L[is.infinite(L)] = 0
  G[is.infinite(G)] = 0
  #L[which.nan(L)] = 0
  #D(is.nan(D)) = 0
  #L(is.infinite(L)) = 0
  #D(is.infinite(D)) = 0
  
  # Construct eigenmaps (solve Ly = lambda*Dy)
  output = geigen(L, D, only.values = FALSE)
  
  # only need bottom (no_dims + 1) eigenvectors
  lambda = output$values
  vectors = output$vectors
  # magnitude = colSums(vectors)
  # vec_index = sort(magnitude, decreasing = TRUE, index.return = TRUE)$ix
  # lambda = lambda[vec_index[1 : (no_dims + 1)]]
  # vectors = vectors[, vec_index[1 : (no_dims + 1)]]

  # Sort eigenvectors in ascending order
  output2 = sort(lambda, decreasing = FALSE, index.return = TRUE)

  # Final embedding
  lambda = lambda[output2$ix[2 : (no_dims + 1)]]
  mappedX = vectors[, output2$ix[2 : (no_dims + 1)]]
  #########
  return(list(K = G, vec = mappedX, val = lambda, X = X, sigma = sigma, k = k, conn_comp = conn_comp))
}