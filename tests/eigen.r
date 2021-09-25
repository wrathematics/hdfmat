library(hdfmat)
set.seed(1234)

f = tempfile()
n = "mydata"
type = "double"

nr = 3
nc = 10
x = matrix(1:(nr*nc), nr, nc)
storage.mode(x) = type

h = crossprod_ooc(x, f, name=n)
test = h$eigen(k=3)

truth = c(9450.2858568857428, 4.714143114255676, 5.0978833483225391e-13)

tol = 1e-6
stopifnot(all.equal(test, truth, tol))

unlink(f)
