library(hdfmat)

f = tempfile()
n = "mydata"
type = "float"

nr = 3
nc = 10
x = matrix(1:(nr*nc), nr, nc)

h = hdfmat:::hdfmat(f, n, nr, nc, type)
h$fill(x)
alpha = 10
h$scale(alpha)

test = h$read()
truth = float::fl(alpha * x)

stopifnot(all.equal(test, truth))

unlink(f)
