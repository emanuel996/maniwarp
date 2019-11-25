manifold_warping <- function(X1, X2, mode, target_dim = NULL, k = 10, mu = 1, thresh = 0.01, max_its = 100){
  library(Rcpp)
  library(pracma)
  library(dtw)
  source('R/aliDif.R')
  source('R/components.R')
  source('R/createKnnGraph.R')
  source('R/graph_laplacian.R')
  source('R/L2_distance.R')
  source('R/laplacian_eigen.R')
  source('R/rowBdSlow.R')
  sourceCpp('src/rowBd.cpp')
  sourceCpp('src/knnsearch.cpp')
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
  dim = min(dim(X1)[1], dim(X2)[1])
  epsilon = 1e-6
  if (is.null(target_dim)){
    target_dim = dim - 1
  }
  # initial alignment
  W1 = createKnnGraph(knnsearch_same(X1, k))
  W2 = createKnnGraph(knnsearch_same(X2, k))
  
  if (strcmp(mode, 'embed')){
    X1 = t(laplacian_eigen(t(X1), target_dim, k))
    X2 = t(laplacian_eigen(t(X2), target_dim, k))
    # check compatibility
    if (dim(X1)[2] != dim(W1)[2] || dim(X2)[2] != dim(W2)[2]){
      stop('embedding didn''t preserve all instances.')
    }
  }
  
  is_linear = strcmp(mode, 'linear')
  
  W12 = Matrix(0, dim(X1)[2], dim(X2)[2])
  W12[1,1] = 1;                         # Align the 2 first points
  W12[dim(X1)[2], dim(X2)[2]] = 1;      # Align the 2 last points
  #[W12,P0] = my_dtw(X1,X2);
  
  # coordinate-descent search
  nIt = 0
  tmp1 = matrix(c(1:dim(X1)[2]), ncol = 1)
  tmp2 = matrix(1, dim(X1)[2], 1)
  P0 = cbind(tmp1, tmp2)
  while (nIt < max_its){
    nIt = nIt + 1
    # manifold alignment
    if (is_linear){
      output1 = manifold_linear(X1, X2, W1, W2, W12, mu, dim+4, epsilon)
      V1 = output1$V1[, 1:target_dim]
      V2 = output1$V1[, 1:target_dim]
      Y1 = t(V1) %*% X1
      Y2 = t(V2) %*% X2
    }else{
      output2 = manifold_nonlinear(W1, W2, W12, mu, dim+4, epsilon)
      Y1 = output2$Y1[, 1:target_dim]
      Y2 = output2$Y2[, 1:target_dim]
    }
    # temporal warping
    output = my_dtw(Y1, Y2)
    W12 = output$W12
    P = output$P
    # stop condition
    dif = aliDif(P, P0)
    if (dif <= thersh){
      break
    }
    P0 = P
  }
  
  if (! is_linear){
    V1 = NULL
    V2 = NULL
  }
  # return outputs
  return(list(P = P, Y1 = Y1, Y2 = Y2, V1 = V1, V2 = V2))
}
