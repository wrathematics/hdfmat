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
#' @rdname crossprod_ooc
#' @export
crossprod_ooc = function(x, file, name="crossprod")
{
  if (!is.matrix(x) && !float::is.float(x))
    x = as.matrix(x)
  
  if (float::is.float(x))
    x = x@Data
  else if (!is.double(x))
    storage.mode(x) = "double"
  
  n = as.double(ncol(x))
  
  h = hdfmat(file, name, n, n, "double")
  h$fill_crossprod(x)
  
  h
}



#' @rdname crossprod_ooc
#' @export
tcrossprod_ooc = function(x, file, name="tcrossprod")
{
  if (!is.matrix(x) && !float::is.float(x))
    x = as.matrix(x)
  
  if (float::is.float(x))
    x = x@Data
  else if (!is.double(x))
    storage.mode(x) = "double"
  
  m = as.double(nrow(x))
  
  h = hdfmat(file, name, m, m, "double")
  h$fill_tcrossprod(x)
  
  h
}
