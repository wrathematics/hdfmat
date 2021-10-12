library(hdfmat)
suppressMessages(library(float))

f = tempfile()
n = "mydata"

nr = 5
nc = 3
h = hdfmat::hdfmat(f, n, nr, nc, "float")

h$fill_linspace(1, nr*nc)

test = h$read()
truth = float::fl(matrix(1:(nr*nc), nr, nc))

stopifnot(all.equal(test, truth))

h$close()
unlink(f)
