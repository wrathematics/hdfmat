f = tempfile()
n = "mydata"
type = "double"

nr = 3
nc = 10
x = matrix(1:(nr*nc), nr, nc)
storage.mode(x) = type

h = hdfmat:::hdfmat(f, n, nc, nc, "double")
h$crossprod(x)

test = h$read()
truth = crossprod(x)

stopifnot(all.equal(test, truth))

unlink(f)
