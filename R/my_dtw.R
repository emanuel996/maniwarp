my_dtw <- function(Y1, Y2){
  library('dtw')
  # Dynamic Time Warping
  # find dtw path and update W12
  newY1 = dtw(Y1, Y2)$index1
  newY2 = dtw(Y1, Y2)$index2
  P = cbind(newY1, newY2)
  W12 = matrix(0, dim(Y1)[1], dim(Y2)[1])
  for (i in 1:length(newY1)){
    W12[newY1[i], newY2[i]] = 1
  }
  #############
  return(list(path = P, matrix = W12))
}