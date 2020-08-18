# derandomKnock
An R package to implement the derandomized knockoffs.

### Overview
`derandomKnock` is an R package to implement the variable selection procedure <em>Derandomized Knockoffs</em>, proposed in our
paper: [Derandomized Knockoffs](). Given the covariate maitrx and response vector, it automatically returns a set of selected variables with type-I 
error guarantee. 


### Installation
To install the package, run the following command in your R console:
```{r}
devtools::install_github("zhimeir/derandomKnock")
```

### Usage Example
We illustrate the usage of `derandomKnock` via a synthetic example. We first generate the data from a linear model.
```{r}
# Generate the data
n=100;p=50;s=10;
rho=0.5;
Sigma=toeplitz(rho^(1:p-1))
X=matrix(rnorm(n*p),n,p)%*%chol(Sigma)
beta=rep(0,p)
beta[1:s]=3.5/sqrt(n)
y=X%*%beta+rnorm(n)
```
Suppose we want a selection set with PFER controlled by 1:
```{r}
res <- derandomKnock(X,y,type = "pfer",v=1, 
                  #type of knockoffs generated
                  knockoff_method = "gaussian",
                  #details about the distribution of X
                  mu = rep(0,p),Sigma = Sigma)
res
```
Suppose we want a selection set with 1-FWER controlled by 0.1:
```{r}
res <- derandomKnock(X,y,type = "kfwer", k=1, alpha = 0.1, 
                    #type of knockoffs generated
                    knockoff_method = "gaussian",
                    #details about the distribution of X
                    mu = rep(0,p),Sigma = Sigma)
res
```
