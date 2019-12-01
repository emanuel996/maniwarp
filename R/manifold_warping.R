manifold_warping <- function(X1, X2, mode = 'linear', target_dim = NULL, k = 8, mu = 0.3, thresh = 0.01, max_its = 100){
  library(Rcpp)
  library(pracma)
  library(dtw)
  #source('R/aliDif.R')
  #source('R/my_components.R')
  #source('R/createKnnGraph.R')
  #source('R/graph_laplacian.R')
  #source('R/L2_distance.R')
  #source('R/laplacian_eigen.R')
  #source('R/knnsearch_slow.R')
  #source('R/rowBdSlow.R')
  #source('R/manifold_linear.R')
  #source('R/manifold_nonlinear.R')
  #source('R/my_dtw.R')
  #sourceCpp('src/rowBd.cpp')
  #sourceCpp('src/knnsearch.cpp')
  # Manifold alignment + DTW
  #
  # Input
  #   X1,X2     - sequences,  each Xi is di x ni.
  #   mode      - one of {'linear','nonlinear','embed'}
  #
  # Output
  #   P        -  new alignment path
  #   Y1,Y2    -  new sequences
  #   V1,V2    -  new mappings (only if mode = 'linear')

  # parameters
  dim = min(dim(X1)[2], dim(X2)[2])
  epsilon = 1e-8
  if (is.null(target_dim)){
    target_dim = dim - 1
  }
  # initial alignment
  W1 = createKnnGraph(knnsearch_slow(X1, k)$index)
  W2 = createKnnGraph(knnsearch_slow(X2, k)$index)

  if (strcmp(mode, 'embed')){
    Y1 = laplacian_eigen(X1, target_dim, k)$vec
    Y2 = laplacian_eigen(X2, target_dim, k)$vec
    # check compatibility
    if (dim(Y1)[2] != dim(W1)[2] || dim(Y2)[2] != dim(W2)[2]){
      print('Embedding didn\'t preserve all instances.')
      return(list(Y1 = Y1, Y2 = Y2))
    }
  }

  is_linear = strcmp(mode, 'linear')

  W12 = matrix(0, dim(X1)[1], dim(X2)[1])
  W12[1,1] = 1                         # Align the 2 first points
  W12[dim(X1)[1], dim(X2)[1]] = 1      # Align the 2 last points
  #[W12,P0] = my_dtw(X1,X2)

  # coordinate-descent search
  nIt = 0
  tmp1 = matrix(c(1:dim(X1)[1]), ncol = 1)
  tmp2 = matrix(1, dim(X1)[1], 1)
  P0 = cbind(tmp1, tmp2)
  while (nIt < max_its){
    nIt = nIt + 1
    # manifold alignment
    if (is_linear){
      output1 = manifold_linear(X1, X2, W1, W2, W12, mu, dim+4, epsilon)
      V1 = output1$map1[, 1:target_dim]
      V2 = output1$map2[, 1:target_dim]
      Y1 = X1 %*% V1
      Y2 = X2 %*% V2
    }else{
      output2 = manifold_nonlinear(X1, X2, W1, W2, W12, mu, dim+4, epsilon)
      Y1 = output2$map1[, 1:target_dim]
      Y2 = output2$map2[, 1:target_dim]
    }
    # temporal warping
    output3 = my_dtw(Y1, Y2)
    W12 = output3$matrix
    P = output3$path
    # stop condition
    dif = aliDif(P, P0)$dif
    if (dif <= thresh){
      break
    }
    P0 = P
  }
  ############
  return(list(P = P, Y1 = Y1, Y2 = Y2))

}
