#include <cstdlib>

#include "hdfmat.h"
#include "extptr.h"
#include "types.h"


extern "C" SEXP R_hdfmat_open(SEXP filename, SEXP fm)
{
  SEXP ret;
  
  // H5::Exception::dontPrint();
  auto mode = INT(fm) == FILE_MODE_CR ? H5F_ACC_TRUNC : H5F_ACC_RDWR;
  H5::H5File *file = new H5::H5File(CHARPT(filename, 0), mode);
  
  newRptr(file, ret, hdf_object_finalizer<H5::H5File>);
  UNPROTECT(1);
  return ret;
}



static inline H5::DSetCreatPropList get_plist(const hsize_t *dim, const int compression)
{
  const hsize_t max_contig_rows = 1;
  
  hsize_t dim_chunk[2];
  dim_chunk[0] = dim[0] > max_contig_rows ? max_contig_rows : dim[0];
  dim_chunk[1] = dim[1];
  
  H5::DSetCreatPropList plist;
  plist.setChunk(2, dim_chunk);
  plist.setDeflate(compression);
  
  return plist;
}

extern "C" SEXP R_hdfmat_init(SEXP fp, SEXP name, SEXP nrows, SEXP ncols, SEXP type, SEXP compression)
{
  SEXP ret;
  
  // H5::Exception::dontPrint();
  H5::H5File *file = (H5::H5File*) getRptr(fp);
  
  hsize_t dim[2];
  dim[0] = DBL(nrows);
  dim[1] = DBL(ncols);
  H5::DataSpace data_space(2, dim);
  
  int cp = INT(compression);
  
  H5::DataSet *dataset = new H5::DataSet;
  if (INT(type) == TYPE_DOUBLE)
  {
    H5::DataType datatype(H5::PredType::IEEE_F64LE);
    if (cp > 0)
    {
      auto plist = get_plist(dim, INT(compression));
      *dataset = file->createDataSet(CHARPT(name, 0), datatype, data_space, plist);
    }
    else
      *dataset = file->createDataSet(CHARPT(name, 0), datatype, data_space);
  }
  else // if (INT(type) == TYPE_FLOAT)
  {
    H5::DataType datatype(H5::PredType::IEEE_F32LE);
    if (cp > 0)
    {
      auto plist = get_plist(dim, INT(compression));
      *dataset = file->createDataSet(CHARPT(name, 0), datatype, data_space, plist);
    }
    else
      *dataset = file->createDataSet(CHARPT(name, 0), datatype, data_space);
  }
  
  newRptr(dataset, ret, hdf_object_finalizer<H5::DataSet>);
  UNPROTECT(1);
  return ret;
}



extern "C" SEXP R_hdfmat_inherit(SEXP fp, SEXP name)
{
  SEXP ds, Rdims, type, ret;
  
  // H5::Exception::dontPrint();
  H5::H5File *file = (H5::H5File*) getRptr(fp);
  
  H5::DataSet *dataset = new H5::DataSet;
  *dataset = file->openDataSet(CHARPT(name, 0));
  
  auto type_class = dataset->getTypeClass();
  if (type_class != H5T_FLOAT)
    error("only float types are supported");
  auto h5flt_type = dataset->getFloatType();
  size_t sz = h5flt_type.getSize();
  
  PROTECT(type = allocVector(INTSXP, 1));
  if (sz == 8)
    INT(type) = TYPE_DOUBLE;
  else if (sz == 4) 
    INT(type) = TYPE_FLOAT;
  else
    error("Unsupported float size");
  
  auto dataspace = dataset->getSpace();
  int ndims = dataspace.getSimpleExtentNdims();
  if (ndims != 2)
    error("invalid number of dimensions in hdf5 file");
  
  hsize_t dims[2];
  dataspace.getSimpleExtentDims(dims, NULL);
  
  newRptr(dataset, ds, hdf_object_finalizer<H5::DataSet>);
  
  PROTECT(Rdims = allocVector(REALSXP, 2));
  for (int i=0; i<ndims; i++)
    REAL(Rdims)[i] = (double) dims[i];
  
  PROTECT(ret = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(ret, 0, ds);
  SET_VECTOR_ELT(ret, 1, Rdims);
  SET_VECTOR_ELT(ret, 2, type);
  
  UNPROTECT(4);
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
