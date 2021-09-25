#' hdfmat class
#' 
#' Storage and methods for CPU matrix data.
#' 
#' @details
#' Data is held in an external pointer.
#' 
#' @useDynLib hdfmat R_hdfmat_cp
#' @useDynLib hdfmat R_hdfmat_eigen_sym
#' @useDynLib hdfmat R_hdfmat_fill
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
      private$file = file
      private$name = name
      private$nrows = as.double(nrows)
      private$ncols = as.double(ncols)
      private$type = type
      
      private$fp = .Call(R_hdfmat_open, file)
      private$ds = .Call(R_hdfmat_init, private$fp, name, nrows, ncols, type)
    },
    
    
    
    print = function()
    {
      cat(paste0("## A ",
        private$nrows, "x", private$ncols, " hdfmat object ", 
        private$file, "\n"))
    },
    
    
    
    fill = function(x)
    {
      x = t(x)
      dim(x) = rev(dim(x))
      .Call(R_hdfmat_fill, private$ds, x)
      invisible(self)
    },
    
    
    
    read = function()
    {
      ret = .Call(R_hdfmat_read, private$nrows, private$ncols, private$ds)
      t(ret)
    },
    
    
    
    scale = function(v)
    {
      v = as.double(v)
      .Call(R_hdfmat_scale, private$nrows, private$ncols, private$ds, v)
      invisible(self)
    },
    
    
    
    set_diag = function(v)
    {
      v = as.double(v)
      .Call(R_hdfmat_set_diag, private$nrows, private$ncols, private$ds, v)
      invisible(self)
    },
    
    
    
    crossprod = function(x)
    {
      storage.mode(x) = "double"
      .Call(R_hdfmat_cp, x, private$ds)
      invisible(self)
    },
    
    
    
    eigen = function(k=3)
    {
      if (private$nrows != private$ncols)
        stop("matrix is non-square")
      
      k = as.integer(k)
      n = as.double(private$nrows)
      .Call(R_hdfmat_eigen_sym, k, n, private$ds)
    }
  ),
  
  
  private = list(
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
