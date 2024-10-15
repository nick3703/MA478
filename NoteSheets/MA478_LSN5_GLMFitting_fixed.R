## Implement IRLS for GLM Fitting

# A family object in R is coded as a list of useful components
pfamily <- poisson()
names(pfamily)
pfamily$variance
pfamily$linkinv

# Iteratively Reweighted Least Squares Algorithm (based on Faraway pg. 155)

myglmfit <- function(y, X, family, tol=1e-8, maxit=50){
  
  # Set initial values
  mu <- y
  oldbeta <- rep(0,ncol(X)) ## setup 'oldbeta' and 'beta' so they're different
  beta <- rep(1,ncol(X))
  it <- 1 ## number of iterations
  
  while (it < maxit && max(abs((1-beta/oldbeta))) > tol) {
    
    oldbeta <- beta                           ## set oldbeta to incumbent
    eta <- family$linkfun(mu)                 ## calc. linear predictor
    mm <- family$mu.eta(eta)                  ## calc. d(mu)/d(eta)
    adjdev <- eta + (y-mu)/mm                 ## form adjusted dependent variable
    W <- c(1/((1/mm)^2*family$variance(mu)))      ## form the weights
    beta <- lm.wfit(X,adjdev,W)$coefficients  ## conduct WLS to re-estimate beta
    mu <- family$linkinv(X %*% beta)          ## compute the new mu
    it <- it+1                                ## update the iteration
  }
  beta
}

# Example: this data set gives the number of warp breaks per loom, where a loom 
# corresponds to a fixed length of yarn.
X <- model.matrix(~wool*tension,
                  data=warpbreaks) # wool, tension and interaction variables
y <- warpbreaks$breaks # the number of breaks (count response)

# Fit the GLM using the IRLS Algorithm
myglmfit(y, X, poisson())
coef(glm(breaks~wool*tension,data=warpbreaks,family=poisson))
