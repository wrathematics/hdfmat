f = tempfile()
n = "mydata"
type = "double"

nr = 3
nc = 5
x = matrix(1:(nr*nc), nr, nc)
storage.mode(x) = type

h = hdfmat:::hdfmat(f, n, nr, nc, type)
h$fill(x)

v = pi
h$set_diag(v)

test = h$read()
truth = x
diag(truth) = v

stopifnot(all.equal(test, truth))

unlink(f)