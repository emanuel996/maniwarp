dataset2 <- function(){
  # sine roll data
  library('plotly')

  num_points = 100
  t = seq(from = 0, to = 10, by = 10/num_points)
  #y1 = cbind(x, sin(x^2), rep(0, num_points))
  #y2 = cbind(x*cos(x), sin(x^2), x*sin(x))
  y1 = data.frame(x = t, y = sin(t^2), z = rep(0, length(t)))
  y2 = data.frame(x = t*cos(t), y = sin(t^2), z = t*sin(t))
  # plot the original figure
  p1 <- plot_ly(y1, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines', name = 'on a plane')
  p2 <- plot_ly(y2, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines', name = 'on a swiss roll')
  p <- subplot(p1, p2)
  # return data set
  return(list(X1 = as.matrix(y1), X2 = as.matrix(y2), p = p) )
}