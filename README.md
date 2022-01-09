# hdfmat

* **Version:** 0.2-3
* **License:** [BSD 2-Clause](https://opensource.org/licenses/BSD-2-Clause)
* **Author:** Drew Schmidt


This package provides some out-of-core matrix operations via HDF5. The main purpose of the package is to provide the architecture to calculate shrink correlation matrices when the number of columns of a data matrix is very large. In this case the resulting correlation matrix can easily exceed the amount of RAM available.

Note that I/O will be the main bottleneck since this is out of core. You could do this much faster using a distributed resource like a cluster, using say the [pbdDMAT package](https://github.com/RBigData/pbdDMAT). But if all you have available is a small compute resource like a workstation, then this can get the job done.



## Installation

You will need a system installation of HDF5 to build the package. The package probably does not build on Windows at the moment. If you know how to set `src/Makevars.win` to link against HDF5, please let me know.

The package also depends several R packages, including the light-weight [R6 package](https://cran.r-project.org/package=R6), and links with the [float](https://cran.r-project.org/package=float) and [fmlh](https://hpcran.org/packages/fmlh/index.html) packages.

You can install the stable version from [the HPCRAN](https://hpcran.org) using the usual `install.packages()`:

```r
install.packages("hdfmat", repos=c("https://hpcran.org", "https://cran.rstudio.com"))
```

The development version of the package is maintained on GitHub.

```r
remotes::install_github("wrathematics/hdfmat")
```



## Basic IO

The main package functions operate out-of-core, that is, on data stored on disk (out of memory) in an HDF5 file. In R, this file is managed as an hdfmat object. There are a few ways of getting the data into the correct format:

1. Build your dataset using other tools (e.g. [hdf5r](https://cran.r-project.org/package=hdf5r)) and "inherit" it as an hdfmat with `hdfmat_open()`.
2. Create your new dataset using the `hdfmat()` constructor, and using the `fill()` method to add blocks of rows to the file.
3. In the specialized case of having a matrix that fits in memory whose crossproduct you need computed out-of-core, you can use `crossprod_ooc()` or `tcrossprod_ooc()`. This makes sense only in the case where you need to calculate
    * `crossprod()` when $m < n$ and $n$ is very large, or
    * `tcrossprod()` when $m > n$ and $m$ is very large, or



## Examples

Here's a a basic example. We'll start by creating a simple matrix:

```r
m = 3
n = 10
x = matrix(1:(m*n), m, n)
x
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
## [1,]    1    4    7   10   13   16   19   22   25    28
## [2,]    2    5    8   11   14   17   20   23   26    29
## [3,]    3    6    9   12   15   18   21   24   27    30
str(x)
## int [1:3, 1:10] 1 2 3 4 5 6 7 8 9 10 ...
```

With our dataset created, we will compute its crossproduct out-of-core. Since we are in the case where we have more columns than rows, this makes sense --- we just need to imagine that the magnitudes of these numbers are far larger.

```r
library(hdfmat)

file = tempfile()
dataset = "mydata"

h = crossprod_ooc(x, file, name=dataset)
h
## # An hdfmat object
##   * Location: /tmp/Rtmp14Voup/file4ade4f50060a
##   * Dimension: 10x10
##   * Type: double
```

Notice that although the input was of type `int` (seen in the call to `str()` above), the storage on disk is type `double`. That's because the default for an hdfmat object is type double. In this case, the input matrix is cast to `double` before constructing the crossproduct out-of-core.

Since the matrix is symmetric, we can compute its eigenvalues. Note that the `eigen()` method does not verify symmetry. The implementation uses the Lanczos method to approximate the first `k` eigenvalues. Here we ask for the first 3.

```r
set.seed(1234)
h$eigen(k=3)
## [1] 9450.285857    4.714143    0.000000
```

Of course, we have the real dataset still in memory, and since it's small we can compare this approximation to the actual values:

```r
eigen(crossprod(x))$values[1:3]
## [1] 9.450286e+03 4.714143e+00 5.237617e-14
```

If we want, we can close the file which will make the hdfmat object handle unusable, although the file still exists and contains data.

```r
h$close()
h
## # An invalid hdfmat object - perhaps it is closed?

system(paste("h5dump", file))
## HDF5 "/tmp/Rtmp14Voup/file4ade4f50060a" {
## GROUP "/" {
##    DATASET "mydata" {
##       DATATYPE  H5T_IEEE_F64LE
##       DATASPACE  SIMPLE { ( 10, 10 ) / ( 10, 10 ) }
##       DATA {
##       (0,0): 14, 32, 50, 68, 86, 104, 122, 140, 158, 176,
##       (1,0): 32, 77, 122, 167, 212, 257, 302, 347, 392, 437,
##       (2,0): 50, 122, 194, 266, 338, 410, 482, 554, 626, 698,
##       (3,0): 68, 167, 266, 365, 464, 563, 662, 761, 860, 959,
##       (4,0): 86, 212, 338, 464, 590, 716, 842, 968, 1094, 1220,
##       (5,0): 104, 257, 410, 563, 716, 869, 1022, 1175, 1328, 1481,
##       (6,0): 122, 302, 482, 662, 842, 1022, 1202, 1382, 1562, 1742,
##       (7,0): 140, 347, 554, 761, 968, 1175, 1382, 1589, 1796, 2003,
##       (8,0): 158, 392, 626, 860, 1094, 1328, 1562, 1796, 2030, 2264,
##       (9,0): 176, 437, 698, 959, 1220, 1481, 1742, 2003, 2264, 2525
##       }
##    }
## }
## }
```

Float types are also supported, either as user inputs via the float package, and/or in hdfmat objects. However, we can't use `crossprod_ooc()` without manually casting the input via `float::fl(x)`. This is easy enough, but instead we will show off the other way to do this:

```r
file = tempfile()
dataset = "mydata"

h = hdfmat(file, name=dataset, n, n, type="float")
h
## # An hdfmat object
##   * Location: /tmp/Rtmp14Voup/file4aded98d9e1
##   * Dimension: 10x10
##   * Type: float
```

Now we can call the `crossprod()` method and the cast will be done for us. The disadvantage here is that we have to manually size the hdfmat object correctly. Incorrectly doing so will cause an error.

```r
h$fill_crossprod(x)

set.seed(1234)
h$eigen(k=3)
## # A float32 vector: 3
## [1] 9.4503e+03 4.7147e+00 3.9768e-04
```
