library(hdfmat)
suppressMessages(library(float))

f = tempfile()
n = "mydata"

nr = 5
nc = 3
h = hdfmat::hdfmat(f, n, nr, nc, "float")

v = 13
h$fill_val(v)

test = h$read()
truth = float::fl(matrix(v, nr, nc))

stopifnot(all.equal(test, truth))

h$close()
unlink(f)
