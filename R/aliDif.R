aliDif <- function(P1, P2){
  # Evaluate the difference between two alignment.
  #
  # Example
  #   input   -  P2 = [1 1 2 3 4 4 4; ...
  #                    1 2 3 4 4 5 6]';
  #           -  P2 = [1 2 2 2 2 2 3 4 4; ...
  #                    1 1 2 3 4 5 5 5 6]';
  #   call    -  [dif, nDif] = aliDif(P1, P2)
  #   output  -  dif = .3, nDif = 8
  #
  # Input
  #   ali1    -  1st alignment
  #   ali2    -  2nd alignment
  #
  # Output
  #   dif     -  difference rate
  #   nDif    -  difference number between two alignment
  
  n = dim(P1)[1]
  n1 = P1[n, 1]
  n2 = P1[n, 2]
  
  # Note: this sometimes segfaults. Try running with rowBdSlow to debug.
  B1 = rowBdSlow(P1)
  B2 = rowBdSlow(P2) 
  #B1 = rowBd(P1)
  #B2 = rowBd(P2) 
  
  nDif = 0;
  for (i in 2:n1){
    gap0 = abs(B1[i-1, 2] - B2[i-1, 2])
    gap = abs(B1[i, 1] - B2[i, 1])
    nDif = nDif + gap0 + 0.5 * (gap0 != gap)
  }
  nAll = (n1 -1) * (n2 -1)
  dif = nDif/nAll
  # return
  #   dif     -  difference rate
  #   nDif    -  difference number between two alignment
  return(list(dif = dif, nDif = nDif))
}