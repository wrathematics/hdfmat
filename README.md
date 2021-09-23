# hdfmat

* **Version:** 0.1-0
* **License:** [BSD 2-Clause](https://opensource.org/licenses/BSD-2-Clause)
* **Author:** Drew Schmidt


This package provides some out-of-core matrix operations via HDF5. The main purpose of the package is to provide the architecture to calculate shrink correlation matrices when the number of columns of a data matrix is very large. In this case the resulting correlation matrix can easily exceed the amount of RAM available.

Note that I/O will be the main bottleneck since this is out of core. You could do this much faster using a distributed resource like a cluster, using say the [pbdDMAT package](https://github.com/RBigData/pbdDMAT). But if all you have available is a small compute resource like a workstation, then this can get the job done.


## Installation

The package depends on the light-weight [R6 package](https://cran.r-project.org/web/packages/R6/index.html) and links with the [fmlr package](https://hpcran.org/packages/fmlr/index.html):

```r
install.packages("R6")
install.packages("fmlr", repos="https://hpcran.org")
```

The development version of the package is maintained on GitHub:

```r
remotes::install_github("wrathematics/hdfmat")
```

You will need a system installation of HDF5 to build the package. The package probably does not build on Windows at the moment. If you know how to set `src/Makevars.win` to link against HDF5, please let me know.


## Package Use

Things are fairly ad hoc at the moment. Currently only double precision data is supported, and there are no checks in place ("yolo mode"). Eventually I will make this more robust and safe, probably.

You can look at the tests folder to see how the package can be used. I will write more when I implement some more important things like SVD.
