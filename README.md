# The Q Test 

This is the github repo for the uniform inference for nonlinear endogenous effects with high-dimensional covariates using double bias correction by Fan et al. (2023).

## Required packages 

You can install the required package with 

``` r
install.pacakges(c("splines","splines2","glmnet","MASS","igraph","matrixStats","Jmisc"))
```
In addition, users need to install the "Rmosek" package manually. The instructions for installation are available at  [Installation of MOSEK Rmosek package.](https://docs.mosek.com/latest/rmosek/install-interface.html)


"Estimate_Functions.R" and "Inference_Functions.R" contain necessary R functions for Lasso estimators and projection vectors for bias correction.  

"example.R" provides an example that is shown below.  


## Example

This is a basic example that shows you how to use the Q test. 

```{r example}
rm(list = ls()) 
source("Estimate_Functions.R")
source("Inference_Functions.R")  


### Setup: 
## Y = g(D) + X\theta + u,
## D = \sum_{\ell} \psi_{\ell}(Z_\ell) + X\phi + v

## Target: Uniform inference for the marginal effect function g'(.)



set.seed(2024)
n = 1000
px = 150
pz = 1
sx = 10 # sparsity index for \theta and \phi 


# Define coefficients and functions 
phi =  c( rep(c(1,-1),sx/2),rep(0,px-sx))
theta = c(rep(1,sx),rep(0,px-sx))

psi.fun <- function(Z, normalization = FALSE, base = NULL){
  if (normalization){
    if (is.null(base)){
      base = Z
    }
    Z = normalize(Z, base = base)
  }
  return( 4 * ( 2 * Z - 1)^2)
}

g.fun <- function(D, normalization = FALSE, base = NULL){
  if (normalization){
    if (is.null(base)){
      base = D 
    }
    D = normalize(D, base = base)
  }
  return( 0.05 * (D-1)^3 )
}
g.deriv <- function(D,normalization = FALSE, base = NULL){ 
  range = 1
  if (normalization){
    if (is.null(base)){
      base = D 
    }
    D = normalize(D, base = base)
    range = max(base) - min(base)
  }
  return( 0.15 * (D-1)^2 ) 
}

## control function: E(u|v) = q(v)
q.fun = function(v, normalization = FALSE, base = NULL){
  if (normalization){
    if (is.null(base)){
      base = v
    }
    v = normalize(v, base = base)
  }
  return( v^2 - 1 ) 
}

### Generate data

Unif <- matrix(runif(n*(px+pz+1)), nrow = n)
UU = Unif[,px+pz+1]

X = ((Unif[,1:px] + 0.3 * UU) / 1.3 - 0.5 ) 
Z =  (Unif[,-c(1:px,px+pz+1)] + 0.3 * UU) / 1.3
v = sqrt(12) * (runif(n) - 0.5)

D = psi.fun(Z)  + X %*% phi + v
u = q.fun(v) + rnorm(n, sd = 1)
Y = g.fun(D) + X %*% theta + u


## grid points
d.seq = seq(-1,3,4/999)

## Uniform inference
out <- Inference.gfun(Y,D,Z,X,d.seq = d.seq, sig.level = 0.05)
```


Plot the results 
```{r}
## Plot the results
g.deriv.plot.dta = cbind(g.deriv(d.seq),out$gderiv.init,out$gderiv.est,
                         out$gderiv.est - out$H.quantile*out$gderiv.sd,
                         out$gderiv.est + out$H.quantile*out$gderiv.sd) 

pdf(file = "Example.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 6) # The height of the plot in inches
 
matplot(d.seq, g.deriv.plot.dta[,1:3], ylab = "g'(D)", xlab = "D",
        type = "l",
        col = c("black","red","blue","blue","blue"), 
        ylim = c(-0.5,1.2))
matlines(d.seq, g.deriv.plot.dta[, 1:3], type = "l", col = c("black","red","blue"),
         lwd = 2,
         lty = 1)
matlines(d.seq, g.deriv.plot.dta[, 4:5], type = "l", col = "blue",
         lwd = 2,
         lty = 2)
legend("bottomright",legend = c("True","Initial","Debiased","95% CB"), 
       col = c("black","red","blue","blue"), 
       lty = c(1,1,1,2,2), lwd = 2)

dev.off()
```


 ![Image load failed](https://raw.github.com/ZiweiMEI/repositpry/master/HDNPIV/Example.pdf) 