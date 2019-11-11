# manifold alignment for linear case
resource("sparse.R")
resource("laplacian.R")

manifold_linear <- function(X1, X2, W1, W2, W12, mu, max_dim, epsilon){
  # Feature-level Manifold Projections. Two domains.
  # X1: P1*M1 matrix, M1 examples in a P1 dimensional space.
  # X2: P2*M2 matrix
  # W1: M1*M1 matrix. weight matrix for each domain.
  # W2: M2*M2 matrix
  # W12: M1*M2 sparse matrix modeling the correspondence of X1 and X2.
  # mu: used to balance matching corresponding pairs and preserving manifold topology.
  # max_dim: max dimensionality of the new space. (default: 200)
  # epsilon: precision. (default: 1e-8)
  p1 = size(X1)[1]
  p2 = size(X2)[1]
  M1 = size(X1)[2]
  M2 = size(X2)[2]
  # get sizes for convenience later
  mu = mu * (sum(W1) + sum(W2))/(2 * sum(W12))
  W = sparse(W1, mu * W12, mu * W12, W2)
  L = laplacian(W)
  # Create weight matrix
  Z = sparse(X1, cbind(P1, M2), cbind(P2, M1), X2)
  svd_X = svd(tcrossprod(X))
  Fplus = pinv(svd_X$u %*% svd_X$v)
  TT = Fplus %*% Z %*% L %*% t(Z) %*% Fplus
  # prepare for decomposition
  vecs = egis( (TT + t(TT))/2)$vec
  vals = egis( (TT + t(TT))/2)$val
  vecs = t(Fplus) %*% vecs
  for (i in 1:nrow(vecs)){
    vecs[, i] = vecs[, i]/norm(vect[, i], 2)
  }
  # Eigen decomposition
  for (i in 1:length(vals)){
    if (vals[i]) > epsilon){
      break
    }
  }
  start = i
  # filter out eigenvalues that are ~= 0
  m = min(max_dim, nrow(vecs) - start +1)
  map1 = vecs[1:P1, start:start+m-1]
  map2 = vecs[P1+1:P1+P2, start:start+m-1]
  # Compute mappings
  return(list(map1, map2))
}