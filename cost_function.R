set.seed(1)
X <- matrix(rnorm(1000), ncol=10) # some random data
Y <- sample(0:1, 100, replace=TRUE)

# Implement Sigmoid function
sigmoid <- function(z) {
  g <- 1/(1+exp(-z))
  return(g)
}

cost.glm <- function(theta,X) {
  m <- nrow(X)
  g <- sigmoid(X%*%theta)
  (1/m)*sum((-Y*log(g)) - ((1-Y)*log(1-g)))
}

X1 <- cbind(1, X)
out<-optim(par=rep(0,ncol(X1)), fn = cost.glm, method='CG',
      X=X1, control=list(trace=TRUE))

glm(Y~X, family=binomial)$coefficients
out
