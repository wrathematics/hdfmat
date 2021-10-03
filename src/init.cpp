#include <cstdlib>

#include "hdfmat.h"
#include "extptr.h"
#include "types.h"


extern "C" SEXP R_hdfmat_open(SEXP filename)
{
  SEXP ret;
  
  // H5::Exception::dontPrint();
  H5::H5File *file = new H5::H5File(CHARPT(filename, 0), H5F_ACC_TRUNC);
  
  newRptr(file, ret, hdf_object_finalizer<H5::H5File>);
  UNPROTECT(1);
  return ret;
}



extern "C" SEXP R_hdfmat_init(SEXP fp, SEXP name, SEXP nrows, SEXP ncols, SEXP type)
{
  SEXP ret;
  
  // H5::Exception::dontPrint();
  H5::H5File *file = (H5::H5File*) getRptr(fp);
  
  hsize_t dim[2];
  dim[0] = DBL(nrows);
  dim[1] = DBL(ncols);
  H5::DataSpace data_space(2, dim);
  
  H5::DataSet *dataset = new H5::DataSet;
  if (INT(type) == TYPE_DOUBLE)
  {
    H5::DataType datatype(H5::PredType::IEEE_F64LE);
    *dataset = file->createDataSet(CHARPT(name, 0), datatype, data_space);
  }
  else // if (INT(type) == TYPE_FLOAT)
  {
    H5::DataType datatype(H5::PredType::IEEE_F32LE);
    *dataset = file->createDataSet(CHARPT(name, 0), datatype, data_space);
  }
  
  newRptr(dataset, ret, hdf_object_finalizer<H5::DataSet>);
  UNPROTECT(1);
  return ret;
}



extern "C" SEXP R_hdfmat_finalize(SEXP fp, SEXP ds)
{
  H5::H5File *file = (H5::H5File*) getRptr(fp);
  H5::DataSet *dataset = (H5::DataSet*) getRptr(ds);
  
  dataset->close();
  file->close();
  
  return R_NilValue;
}
