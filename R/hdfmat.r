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
#' @useDynLib hdfmat R_hdfmat_finalize
#' @useDynLib hdfmat R_hdfmat_init
#' @useDynLib hdfmat R_hdfmat_open
#' @useDynLib hdfmat R_hdfmat_read
#' @useDynLib hdfmat R_hdfmat_scale
#' @useDynLib hdfmat R_hdfmat_set_diag
#' 
#' @rdname hdfmat-class
#' @name hdfmat-class
hdfmatR6 = R6::R6Class("cpumat",
  public = list(
    #' @details
    #' Class initializer.
    #' @param file File to store data in.
    #' @param name Dataset name on disk.
    #' @param nrows,ncols The dimension of the matrix.
    #' @param type Storage type for the matrix. Should be one of 'int', 'float', or 'double'.
    initialize = function(file, name, nrows, ncols, type="double")
    {
      type = match.arg(tolower(type), c("double", "float"))
      type = type_str2int(type)
      
      private$file = file
      private$name = name
      private$nrows = as.double(nrows)
      private$ncols = as.double(ncols)
      private$type = type
      
      private$fp = .Call(R_hdfmat_open, file)
      private$ds = .Call(R_hdfmat_init, private$fp, name, nrows, ncols, type)
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
    #' Fill an hdfmat file with the input matrix. Really only meant for testing.
    #' @param x The input matrix. Must be double or float.
    fill = function(x)
    {
      if (!is.matrix(x) && !is.float(x))
        x = as.matrix(x)
      
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
      .Call(R_hdfmat_fill, private$ds, x, private$type)
      invisible(self)
    },
    
    
    #' @details
    #' Read an hdfmat-stored matrix into memory. Really only meant for testing.
    read = function()
    {
      ret = .Call(R_hdfmat_read, private$nrows, private$ncols, private$ds, private$type)
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
    #' Set the diagonal of an hdfmat-stored matrix to the input scalar.
    #' @param v Scalar. Fundamental type can be double, float, or int.
    set_diag = function(v)
    {
      v = as.double(v)
      .Call(R_hdfmat_set_diag, private$nrows, private$ncols, private$ds, v, private$type)
      invisible(self)
    },
    
    
    #' @details
    #' Calculate the crossproduct of an input matrix with result stored in an
    #' hdfmat. Useful when the number of columns of the input is very large.
    #' @param x Input matrix. Fundamental type can be double, float, or int.
    crossprod = function(x)
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
      }
      
      .Call(R_hdfmat_cp, x, private$ds, private$type)
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
    }
  ),
  
  
  private = list(
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
    type = "",
    fp = NULL,
    ds = NULL
  )
)



#' hdfmat
#' 
#' Constructor for cpumat objects.
#' 
#' @details
#' Data is held in an external pointer.
#' 
#' @param nrows,ncols The dimensions of the matrix.
#' @param type Storage type for the matrix. Should be one of 'int', 'float', or 'double'.
#' @return A cpumat class object.
#' 
#' @export
hdfmat = function(file, name, nrows, ncols, type="double")
{
  hdfmatR6$new(file=file, name=name, nrows=nrows, ncols=ncols, type=type)
}
