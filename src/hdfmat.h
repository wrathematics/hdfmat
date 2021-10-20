#ifndef HDFMAT_H
#define HDFMAT_H
#pragma once


#include <H5Cpp.h>
#include <H5DataSpace.h>

#include <R.h>
#include <Rinternals.h>

#include <float/float32.h>

#define CHARPT(x,i) ((char*)CHAR(STRING_ELT(x,i)))
#define INT(x) (INTEGER(x)[0])
#define DBL(x) (REAL(x)[0])

#define TRY_CATCH(expr) try { expr; } catch(const std::exception& e) { error(e.what()); }


#endif
