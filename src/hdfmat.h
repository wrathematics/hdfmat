#ifndef HDFMAT_H
#define HDFMAT_H
#pragma once


#include <H5Cpp.h>
#include <H5DataSpace.h>

#include <R.h>
#include <Rinternals.h>

#define CHARPT(x,i) ((char*)CHAR(STRING_ELT(x,i)))
#define INT(x) (INTEGER(x)[0])
#define DBL(x) (REAL(x)[0])


#endif
