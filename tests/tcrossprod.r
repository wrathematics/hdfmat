library(hdfmat)

f = tempfile()
n = "mydata"
type = "double"

nr = 10
nc = 3
x = matrix(1:(nr*nc), nr, nc)
storage.mode(x) = type

h = tcrossprod_ooc(x, f, name=n)

test = h$read()
truth = tcrossprod(x)

stopifnot(all.equal(test, truth))

h$close()
unlink(f)
