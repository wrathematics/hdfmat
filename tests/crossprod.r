library(hdfmat)

f = tempfile()
n = "mydata"
type = "double"

nr = 3
nc = 10
x = matrix(1:(nr*nc), nr, nc)
storage.mode(x) = type

h = crossprod_ooc(x, f, name=n)

test = h$read()
truth = crossprod(x)

stopifnot(all.equal(test, truth))

unlink(f)
