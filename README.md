# hdfmat

* **Version:** 0.1-0
* **License:** [BSD 2-Clause](https://opensource.org/licenses/BSD-2-Clause)
* **Author:** Drew Schmidt


This package provides some out-of-core matrix operations via HDF5. The main purpose of the package is to provide the architecture to calculate shrink correlation matrices when the number of columns of a data matrix is very large. In this case the resulting correlation matrix can easily exceed the amount of RAM available.

Note that I/O will be the main bottleneck since this is out of core. You could do this much faster using a distributed resource like a cluster, using say the [pbdDMAT package](https://github.com/RBigData/pbdDMAT). But if all you have available is a small compute resource like a workstation, then this can get the job done.


## Installation

The package depends on the light-weight [R6 package](https://cran.r-project.org/web/packages/R6/index.html) and links with the [float](https://cran.r-project.org/web/packages/float/index.html) and [fmlr](https://hpcran.org/packages/fmlr/index.html) packages:

```r
repos = c("https://hpcran.org", "https://cran.rstudio.com")
install.packages("R6", repos=repos)
install.packages("float", repos=repos)
install.packages("fmlr", repos=repos)
```

The development version of the package is maintained on GitHub:

```r
remotes::install_github("wrathematics/hdfmat")
```

You will need a system installation of HDF5 to build the package. The package probably does not build on Windows at the moment. If you know how to set `src/Makevars.win` to link against HDF5, please let me know.


## Package Use

Things are fairly ad hoc at the moment. Currently only double precision data is supported, and there are no checks in place ("yolo mode"). Eventually I will make this more robust and safe, probably.

Here's a a basic example:

```r
type = "double"

nr = 3
nc = 10
x = matrix(1:(nr*nc), nr, nc)
storage.mode(x) = type
x
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
## [1,]    1    4    7   10   13   16   19   22   25    28
## [2,]    2    5    8   11   14   17   20   23   26    29
## [3,]    3    6    9   12   15   18   21   24   27    30
str(x)
## num [1:3, 1:10] 1 2 3 4 5 6 7 8 9 10 ...
```

With our dataset created, we will compute its crossproduct out-of-core. Since we are in the case where we have more columns than rows, this makes sense --- we just need to imagine that the magnitudes of these numbers are far larger.

```r
library(hdfmat)

file = tempfile()
dataset = "mydata"

h = crossprod_ooc(x, file, name=dataset)
h
## ## A 10x10 hdfmat object /tmp/RtmpvpfPca/file3b43ef31938bd2
```

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
