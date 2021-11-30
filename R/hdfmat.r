#' hdfmat class
#' 
#' Storage and methods for matrix data in an HDF5 file.
#' 
#' @details
#' Data is held in an external pointer.
#' 
#' @useDynLib hdfmat R_hdfmat_cp
#' @useDynLib hdfmat R_hdfmat_eigen_sym
#' @useDynLib hdfmat R_hdfmat_fill
#' @useDynLib hdfmat R_hdfmat_fill_diag
#' @useDynLib hdfmat R_hdfmat_fill_linspace
#' @useDynLib hdfmat R_hdfmat_fill_rnorm
#' @useDynLib hdfmat R_hdfmat_fill_runif
#' @useDynLib hdfmat R_hdfmat_fill_val
#' @useDynLib hdfmat R_hdfmat_finalize
#' @useDynLib hdfmat R_hdfmat_inherit
#' @useDynLib hdfmat R_hdfmat_init
#' @useDynLib hdfmat R_hdfmat_open
#' @useDynLib hdfmat R_hdfmat_read
#' @useDynLib hdfmat R_hdfmat_scale
#' @useDynLib hdfmat R_hdfmat_svd
#' @useDynLib hdfmat R_hdfmat_tcp
#' 
#' @rdname hdfmat-class
#' @name hdfmat-class
hdfmatR6 = R6::R6Class("cpumat",
  public = list(
    #' @details
    #' Class initializer.
    #' @param open Are you working with an existing matrix? If so, \code{nrows},
    #' \code{ncols}, and \code{type} are ignored (and intuited).
    #' @param file File to store data in.
    #' @param name Dataset name on disk.
    #' @param nrows,ncols The dimension of the matrix.
    #' @param type Storage type for the matrix. Should be one of 'int', 'float',
    #' or 'double'.
    #' @param compression The compression level, an integer from 0 (no compression)
    #' to 9 (highest compression). Run-time performance degrades with increased
    #' compression levels.
    initialize = function(open, file, name, nrows, ncols, type, compression)
    {
      file = normalizePath(file, winslash="/", mustWork=FALSE)
      
      if (isTRUE(open))
        private$inherit(file=file, name=name)
      else
        private$create(file=file, name=name, nrows=nrows, ncols=ncols, type=type, compression=compression)
      
      invisible(self)
    },
    
    
    #' @details
    #' Closes the file.
    close = function()
    {
      private$finalize()
    },
    
    
    #' @details
    #' Print some basic info about an hdfmat object.
    print = function()
    {
      if (is.null(private$fp))
        cat(paste0("# An invalid hdfmat object - perhaps it is closed?\n"))
      else
        cat(paste0("# An hdfmat object\n", 
          "  * Location: ", private$file, "\n",
          "  * Dimension: ", private$nrows, "x", private$ncols, "\n",
          "  * Type: ", type_int2str(private$type), "\n",
          "\n"))
    },
    
    
    #' @details
    #' Fill an hdfmat matrix (or part of it) with the input matrix. The number
    #' of rows of the input must be less than or equal to the number of rows
    #' allocated for the hdfmat, and the number of columns must be equal.
    #' @param x The input matrix. Must be double or float.
    #' @param row_offset The number of rows to skip. Useful when building the
    #' hdfmat in chunks.
    fill = function(x, row_offset=0)
    {
      if (!is.matrix(x) && !is.float(x))
        x = as.matrix(x)
      
      if (!is.numeric(row_offset) || is.na(row_offset) || length(row_offset) != 1 || row_offset < 0)
        stop("row_offset should be an int")
      else if (private$nrows < nrow(x) || private$ncols != ncol(x))
        stop("dimensions of x and hdfmat are non-conformable")
      else if (nrow(x) + row_offset > private$nrows)
        stop("x and given row_offset incompatible with hdfmat")
      
      row_offset = floor(as.double(row_offset))
      
      if (private$type == TYPE_DOUBLE)
      {
        if (float::is.float(x))
          x = float::dbl(x)
        else if (typeof(x) != "double")
          storage.mode(x) = "double"
      }
      else # if (private$type == TYPE_FLOAT)
      {
        if (!float::is.float(x))
          x = float::fl(x)@Data
      }
      
      x = t(x)
      dim(x) = rev(dim(x))
      .Call(R_hdfmat_fill, private$ds, x, row_offset, private$type)
      invisible(self)
    },
    
    
    #' @details
    #' Read an hdfmat-stored matrix into memory.
    #' @param row_start,row_stop The first/last row (1-based) to read. If
    #' missing, the values 1 and total number of rows will be used,
    #' respectively.
    read = function(row_start, row_stop, col_start, col_stop)
    {
      if (missing(row_start))
        row_start = 1
      if (missing(row_stop))
        row_stop = private$nrows
      if (missing(col_start))
        col_start = 1
      if (missing(col_stop))
        col_stop = private$ncols
      
      if (!is.numeric(row_start) || !is.numeric(row_stop) || length(row_start) != 1 || length(row_stop) != 1)
        stop("'row_start' and 'row_stop' must be single numbers")
      if (!is.numeric(col_start) || !is.numeric(col_stop) || length(col_start) != 1 || length(col_stop) != 1)
        stop("'col_start' and 'col_stop' must be single numbers")
      
      if (row_stop < row_start || row_start < 1 || row_stop < 1 || row_start > private$nrows)
        stop("must have 1 <= row_start <= row_stop <= nrows")
      if (col_stop < col_start || col_start < 1 || col_stop < 1 || col_start > private$ncols)
        stop("must have 1 <= col_start <= col_stop <= ncols")
      
      row_start = as.double(row_start) - 1.0
      row_stop = as.double(row_stop) - 1.0
      col_start = as.double(col_start) - 1.0
      col_stop = as.double(col_stop) - 1.0
      
      ret = .Call(R_hdfmat_read, row_start, row_stop, col_start, col_stop, private$nrows, private$ncols, private$ds, private$type)
      if (private$type == TYPE_FLOAT)
        ret = float::float32(ret)
      
      t(ret)
    },
    
    
    #' @details
    #' Scale (multiply) all values of an hdfmat-stored matrix by the input
    #' scalar.
    #' @param v Scalar. Fundamental type can be double, float, or int.
    scale = function(v)
    {
      v = as.double(v)
      .Call(R_hdfmat_scale, private$nrows, private$ncols, private$ds, v, private$type)
      invisible(self)
    },
    
    
    #' @details
    #' Set every element to the given scalar value.
    #' @param v Scalar. Fundamental type can be double, float, or int.
    fill_val = function(v)
    {
      v = as.double(v)
      .Call(R_hdfmat_fill_val, private$nrows, private$ncols, private$ds, v, private$type)
      invisible(self)
    },
    
    
    #' @details
    #' Fill the matrix (column-wise) with linearly-spaced values.
    #' @param start,stop Beginning/end of the linear spacing.
    fill_linspace = function(start, stop)
    {
      if (start == stop)
        self$fill_val(start)
      else
      {
        start = as.double(start)
        stop = as.double(stop)
        .Call(R_hdfmat_fill_linspace, private$nrows, private$ncols, private$ds, start, stop, private$type)
      }
      
      invisible(self)
    },
    
    
    #' @details
    #' Fill the matrix with random uniform values.
    #' @param min,max Minimum/maximum values for the generator.
    fill_runif = function(min=0, max=1)
    {
      if (min == max)
        self$fill_val(min)
      else if (min < max)
      {
        min = as.double(min)
        max = as.double(max)
        .Call(R_hdfmat_fill_runif, private$nrows, private$ncols, private$ds, min, max, private$type)
      }
      else
        stop("need min <= max")
      
      invisible(self)
    },
    
    
    #' @details
    #' Fill the matrix with random normal values.
    #' @param mean,sd Mean/standard deviation values for the generator.
    fill_rnorm = function(mean=0, sd=1)
    {
      if (sd == 0)
        self$fill_val(mean)
      else if (sd > 0)
      {
        mean = as.double(mean)
        sd = as.double(sd)
        .Call(R_hdfmat_fill_rnorm, private$nrows, private$ncols, private$ds, mean, sd, private$type)
      }
      else
        stop("need sd >= 0")
      
      invisible(self)
    },
    
    
    #' @details
    #' Set the diagonal of an hdfmat-stored matrix to the input vector. Values
    #' will be recycled as necessary.
    #' @param v A vector. Fundamental type can be double, float, or int.
    fill_diag = function(v)
    {
      v = as.double(v)
      .Call(R_hdfmat_fill_diag, private$nrows, private$ncols, private$ds, v, private$type)
      invisible(self)
    },
    
    
    #' @details
    #' Calculate the crossproduct of an input matrix with result stored in an
    #' hdfmat. Useful when the number of columns of the input is very large.
    #' @param x Input matrix. Fundamental type can be double, float, or int.
    fill_crossprod = function(x)
    {
      n = ncol(x)
      if (n != private$nrows || n != private$ncols)
        stop(paste0("hdfmat dimension ", private$nrows, "x", private$ncols, " different from crossprod of input ", n, "x", n))
      
      if (private$type == TYPE_DOUBLE)
      {
        if (float::is.float(x))
          x = float::dbl(x)
        else if (typeof(x) != "double")
          storage.mode(x) = "double"
      }
      else # if (private$type == TYPE_FLOAT)
      {
        if (!float::is.float(x))
          x = float::fl(x)@Data
        else
          x = x@Data
      }
      
      .Call(R_hdfmat_cp, x, private$ds, private$type)
      invisible(self)
    },
    
    
    #' @details
    #' Calculate the transposed crossproduct of an input matrix with result
    #' stored in an hdfmat. Useful when the number of columns of the input is
    #' very large.
    #' @param x Input matrix. Fundamental type can be double, float, or int.
    fill_tcrossprod = function(x)
    {
      m = nrow(x)
      if (m != private$nrows || m != private$ncols)
        stop(paste0("hdfmat dimension ", private$nrows, "x", private$ncols, " different from crossprod of input ", m, "x", m))
      
      if (private$type == TYPE_DOUBLE)
      {
        if (float::is.float(x))
          x = float::dbl(x)
        else if (typeof(x) != "double")
          storage.mode(x) = "double"
      }
      else # if (private$type == TYPE_FLOAT)
      {
        if (!float::is.float(x))
          x = float::fl(x)@Data
        else
          x = x@Data
      }
      
      .Call(R_hdfmat_tcp, x, private$ds, private$type)
      invisible(self)
    },
    
    
    #' @details
    #' Compute approximations to the eigenvalues of a square symmetric
    #' hdfmat-stored matrix using the Lanczos method. The matrix is not checked
    #' for symmetry.
    #' @param k The number of Lanczos iterations.
    eigen = function(k=3)
    {
      if (private$nrows != private$ncols)
        stop("matrix is non-square")
      
      k = as.integer(k)
      n = as.double(private$nrows)
      values = .Call(R_hdfmat_eigen_sym, k, n, private$ds, private$type)
      if (private$type == TYPE_FLOAT)
        values = float::float32(values)
      
      values
    },
    
    
    #' @details
    #' Compute approximations to the singular values of a rectangular
    #' hdfmat-stored matrix using the Lanczos method.
    #' @param k The number of Lanczos iterations.
    svd = function(k=3)
    {
      k = as.integer(k)
      values = .Call(R_hdfmat_svd, k, private$nrows, private$ncols, private$ds, private$type)
      if (private$type == TYPE_FLOAT)
        values = float::float32(values)
      
      values
    }
  ),
  
  
  
  private = list(
    open = function(file, name, mode)
    {
      private$file = file
      private$name = name
      
      private$fp = .Call(R_hdfmat_open, file, mode)
    },
    
    
    inherit = function(file, name)
    {
      private$file = file
      private$name = name
      
      private$open(file=file, name=name, mode=FILE_MODE_RW)
      
      ret = .Call(R_hdfmat_inherit, private$fp, name)
      private$ds = ret[[1]]
      private$nrows = ret[[2]][1]
      private$ncols = ret[[2]][2]
      private$type = ret[[3]]
    },
    
    
    create = function(file, name, nrows, ncols, type, compression)
    {
      type = match.arg(tolower(type), c("double", "float"))
      type = type_str2int(type)
      
      compression = as.integer(compression)
      if (!(compression %in% 0L:9L))
        stop("'compression' must be an integer from 0 to 9")
      
      nrows = as.double(nrows)
      ncols = as.double(ncols)
      
      private$nrows = nrows
      private$ncols = ncols
      private$type = type
      
      private$open(file=file, name=name, mode=FILE_MODE_CR)
      private$ds = .Call(R_hdfmat_init, private$fp, name, nrows, ncols, type, compression)
    },
    
    
    finalize = function()
    {
      if (is.null(private$fp))
        return(invisible(self))
      
      .Call(R_hdfmat_finalize, private$fp, private$ds)
      
      private$ds = NULL
      private$fp = NULL
      invisible(gc())
      
      invisible(self)
    },
    
    file = "",
    name = "",
    nrows = 0,
    ncols = 0,
    type = 0L,
    fp = NULL,
    ds = NULL
  )
)



#' hdfmat
#' 
#' Constructor for hdfmat objects.
#' 
#' @param file File to store data in.
#' @param name Dataset name on disk.
#' @param nrows,ncols The dimension of the matrix.
#' @param type Storage type for the matrix. Should be one of 'int', 'float', or
#' 'double'.
#' @param compression The compression level, an integer from 0 (no compression)
#' to 9 (highest compression). Run-time performance degrades with increased
#' compression levels.
#' 
#' @return An hdfmat class object.
#' 
#' @export
hdfmat = function(file, name, nrows, ncols, type="double", compression=0L)
{
  hdfmatR6$new(open=FALSE, file=file, name=name, nrows=nrows, ncols=ncols, type=type, compression=compression)
}



#' hdfmat_open
#' 
#' Constructor for hdfmat objects.
#' 
#' @param file File to store data in.
#' @param name Dataset name on disk.
#' @return An hdfmat class object.
#' 
#' @export
hdfmat_open = function(file, name)
{
  hdfmatR6$new(open=TRUE, file=file, name=name)
}
