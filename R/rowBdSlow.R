rowBdSlow <- function(P){
  # Bound the correspondence in rows.
  # Example
  #   input   -  P = [1 1 2 3 4 4 4; ...
  #                   1 2 3 4 4 5 6]^T;
  # call    -  B = rowBd(C)
  # output  -  B = [1 3 4 4; ...
  #                 2 3 4 6]^T;
  #
  # Input
  #   P       -  correspondence matrix, n0 x 2
  #
  # Output
  #   B       -  boundary of each row, n1 x 2
  n0 = dim(P)[1]
  n1 = P[n0, 1]
  B = matrix(0, n1, 2)
  
  head = 1
  while (head <= n0){
    i = P[head, 1]
    
    tail = head
    while (tail <= n0 && P[tail, 1] == i){
      tail = tail +1
    }
    B[i, ] = P[c(head, tail - 1), 2]
    head = tail
  }
  return(B)
  # example2
  # P = [1 2 2 2 2 2 3 4 4
  #      1 1 2 3 4 5 5 5 6]
  # B = [1 1 5 5
  #      1 5 5 6]

}