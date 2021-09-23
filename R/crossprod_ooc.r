#' crossprod_ooc
#' 
#' Out-of-core crossproduct. Useful when the number of columns is very large.
#' 
#' @param x
#' The input matrix. Should be in double precision.
#' @param file
#' Name of the file to use for the out-of-core storage.
#' @param name
#' The dataset name within the HDF5 file.
#' 
#' @return Returns an hdfmat object.
#' 
#' @export
crossprod_ooc = function(x, file, name="crossprod")
{
  if (!is.matrix(x) || !is.double(x))
    stop("'x' must be a numeric matrix")
  
  n = as.double(ncol(x))
  
  h = hdfmat:::hdfmat(file, name, n, n, "double")
  h$crossprod(x)
  
  h
}
