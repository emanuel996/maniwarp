# Manifold_Warping

Knowledge transfer is computationally challenging, due in part to the curse of dimensionality. Recent work on manifold learning has shown that data collected in real-world settings often have high-dimensional representations but lie on low-dimensional manifolds.  
This package is designed to align two sequentially ordered high-dimensional data sets by combining traditional manifold alignment and dynamic time warping algorithms. In each iteration, it firstly aligns two data sets by graph Laplacian and then uses a dynamic time warping method to pair them. Finally, it updates the lost matrix. One can choose linear or non-linear Laplacian method by parameter mode (\'linear\' or \'nonlinear\'). In order to compare with embedding-only algorithm, it also can show embedding by choosing mode \'embed\'. The idea and theoretical formulation are from https://people.cs.umass.edu/~ccarey/pubs/ManifoldWarping.pdf.

## Installation:

You can use following codes to install the package. 
```{r}
## You may need following codes to install dependent packages.
library(devtools)
install_github("emanuel996/maniwarp")
```
After this, we can use this package.
```{r}
library(maniwarp)
```
## Main functions:

|Code| Function |
|:-|:-|
|manifold_warping|Performing manifold alignment and warping|
|manifold_linear|Performing linear manifold alignment|
|manifold_nonlinear|Performing non-linear manifold alignment|
|laplacian_eigen|Performing manifold embedding|
|my_dtw|Performing dynamic time warping|

## Toy example:
Here we use the pre-designed synthetic data set (sine functions) to show our main function \"manifold_warping\". 

```{r}
X1 = dataset2()$X1
X2 = dataset2()$X2
```

It will return a list with 3 matrices. The first matrix is the warping path. Here is a part of the warping path.

```
            [1] [2] [3] [4] [5] [6] [7] [8] [9] [10] [11] [12] [13] [14] [15] [16]
newY1    1    2    3    4    5    6    6    6    6     6     6     6     6     7     7     7
newY2    1    1    1    1    1    2    3    4    5     6     7     8     9    10    11    12
```

The second and third matrices are the projections from the original data to the low-dimensional representation. Instead of printing the matrices, we can visualize the matrices.


For more details, please read [vignettes](https://github.com/emanuel996/maniwarp/blob/master/vignettes/manifoldwarping_pdf.pdf).




