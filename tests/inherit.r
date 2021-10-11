library(hdfmat)
suppressMessages(library(float))

f = tempfile()
n = "mydata"
type = "float"

nr = 10
nc = 5
x = matrix(1:(nr*nc), nr, nc)

h = hdfmat::hdfmat(f, n, nr, nc, type)
h$fill(x)
h$close()

k = hdfmat::hdfmat_open(f, n)

test = k$read()
truth = float::fl(x)

stopifnot(all.equal(test, truth))
unlink(f)
