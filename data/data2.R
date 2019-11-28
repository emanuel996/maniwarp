# sine roll data
library('plotly')

num_points = 100
t = linspace(0, 10, num_points)
#y1 = cbind(x, sin(x^2), rep(0, num_points))
#y2 = cbind(x*cos(x), sin(x^2), x*sin(x))
y1 = data.frame(x = t, y = sin(t^2), z = rep(0, num_points))
y2 = data.frame(x = t*cos(t), y = sin(t^2), z = t*sin(t))

plot_ly(y1, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines')
plot_ly(y2, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines')
