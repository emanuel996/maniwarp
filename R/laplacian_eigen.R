laplacian_eigen <- function(X, no_dims = 2, k = 12, sigma = 1){
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
  source('L2_distance.R')
  source('components.R')
  library('FNN')
  library('geigen')
  
  # construct neighborhood graph
  G = get.knn(X, k)
  G = G^2
  G = G/max(G)
  # Only embed largest connected component of the neighborhood graph
  blocks = components(G)
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
  L = D - G;
  L(is.nan(L)) = 0; D(is.nan(D)) = 0;
  L(is.infinite(L)) = 0; D(is.infinite(D)) = 0;
  
  # Construct eigenmaps (solve Ly = lambda*Dy)
  output = geigen(L, D, only.values = FALSE)
  
  # only need bottom (no_dims + 1) eigenvectors
  lambda = output$values
  vectors = output$vectors
  magnitude = colSums(vectors)
  vec_index = sort(magnitude, decreasing = TRUE, index.return = TRUE)$ix
  lambda = lambda[vec_index[1 : (no_dims + 1)]]
  vectors = vectors[, vec_index[1 : (no_dims + 1)]]

  # Sort eigenvectors in ascending order
  output2 = sort(lambda, decreasing = FALSE, index.return = TRUE)

  # Final embedding
  lambda = lambda[output2$ix[2 : (no_dims + 1)]]
  mappedX = vectors[, output2$ix[2 : (no_dims + 1)]]
  #########
  return(list(K = G, vec = mappedX, val = lambda, X = X, sigma = sigma, k = k, conn_comp = conn_comp))
}