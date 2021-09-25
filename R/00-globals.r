TYPE_DOUBLE = 1L
TYPE_FLOAT = 2L

TYPES_STR = c("double", "float")
TYPES_INT = 1:length(TYPES_STR)
names(TYPES_INT) = TYPES_STR

type_int2str = function(type) TYPES_STR[type]
type_str2int = function(type) TYPES_INT[[type]]

type_robj2int = function(robj)
{
  if (is.double(robj))
    TYPE_DOUBLE
  else if (float::is.float(robj))
    TYPE_FLOAT
  else
    stop("bad fundamental type")
}

type_robj2str = function(robj) type_int2str(type_robj2int(robj))
