f = tempfile()
n = "mydata"
type = "double"

nr = 3
nc = 5
x = matrix(1:(nr*nc), nr, nc)
storage.mode(x) = type

h = hdfmat:::hdfmat(f, n, nr, nc, type)
h$fill(x)

alpha = 2.0
h$scale(alpha)

test = h$read()
truth = alpha*x

stopifnot(all.equal(test, truth))

unlink(f)
