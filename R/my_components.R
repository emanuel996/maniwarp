my_components <-function(A){
  # COMPONENTS Finds connected components in a graph defined by a adjacency matrix
  
  # The function outputs an n-vector of integers 1:k in blocks, meaning that
  # A has k components. The vector blocks labels the vertices of A according 
  # to component.
  # If the adjacency matrix A is undirected (i.e. symmetric), the blocks are 
  # its connected components. If the adjacency matrix A is directed (i.e. 
  # unsymmetric), the blocks are its strongly connected components.
  library('igraph')
  # Adjacency matrix must be square
  n = dim(A)[1]
  m = dim(A)[2]
  if (n != m){
    stop('Adjacency matrix must be square')
  }
  g = graph.adjacency(A)
  # we may use plot(g) for visulization
  # plot(g)
  block_ind = components(g)$csize
  return(block_ind)
}