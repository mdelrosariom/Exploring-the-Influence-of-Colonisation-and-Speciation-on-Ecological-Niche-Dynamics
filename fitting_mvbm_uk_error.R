# Multivariate BM process with known and unknown ME
mvbm <- function(tree, data, error=NULL, nuisance=TRUE,...){
  
  if(is.null(error)) error = matrix(0,ncol=ncol(data),nrow=nrow(data))
  # The model is matrix variate normal: X~MN(B, C, R) where B are the ancestral states,
  # C is the phylogenetic covariance matrix, and R is the (evolutionary) traits covariances
  # here I use the "pic" algorithm that uses Felsenstein pruning to speed up the computations
  descendent <- tree$edge[,2]
  extern <- (descendent <= Ntip(tree))
  height = max(nodeHeights(tree))
  
  fun <- function(mserr){
    
    tree_list <- lapply(1:ncol(data), function(i) {
      tree2 = tree
      tree2$edge.length[extern] <-  tree$edge.length[extern] + (error[,i] + mserr^2)/(var(data[,i]) / height) # trick here to scale the relative variance at the tips compared to the process variance. It's a rough approximation...
      tree2
      })
    
    ll <- mvLL(tree_list, data, method="pic", param=list(estim=TRUE))
    if(is.finite(ll$logl)) return(-ll$logl) else return(1e6) # remember that optim is minimizing the function, we should take the negative log-likelihood here
  }
  
  # should I estimate the nuisance parameter?
  if(nuisance){
    # Optimization of the function to find the "nuisance" parameter - I use 1/10 of the average trait variance as a starting value, but ideally this should be explored
    opt <- optim(par=0.1*mean(apply(data,2,var)), fn = fun, method = "L-BFGS-B")
  
    # retrieve the parameters
    mserr = opt$par
    tree_list <- lapply(1:ncol(data), function(i) {
      tree2 = tree
      tree2$edge.length[extern] <-  tree$edge.length[extern] + (error[,i] + mserr^2)/(var(data[,i]) / height)
      tree2
    })
  }else{
    mserr = 0
    tree_list <- lapply(1:ncol(data), function(i) {
      tree2 = tree
      tree2$edge.length[extern] <-  tree$edge.length[extern] + (error[,i] + mserr^2)/(var(data[,i]) / height)
      tree2
    })
  }
  
  # fit the model (with estimated nuisance parameter if needed)
  multiv_bm <- mvLL(tree_list, data, method="pic", param=list(estim=TRUE))
  
  results <- list(mserr_sq=mserr^2, sigma=multiv_bm$sigma, theta=multiv_bm$theta, logl=multiv_bm$logl, 
                  param=list(model="BM1", ntraits=ncol(data), nregimes=1, smean = TRUE, traits=colnames(data)))
  class(results) <- "mvmorph"
  return(results)
 
}



## Test example

# Simulated dataset
set.seed(14)
# Generating a random tree
tree<-pbtree(n=150)

# Simulate the traits
sigma<-matrix(c(0.1,0.05,0.05,0.1),2)
theta<-c(0,0)
data<-mvSIM(tree, param=list(sigma=sigma, theta=theta, model="BM1", nsim=1))

# fit the model
fit <- mvbm(tree, data)



# add some (isotropic) noise to the data
data2 = data + matrix(rnorm(2*Ntip(tree), sd=sqrt(0.05)), ncol=2)
mvbm(tree, data2)
mvbm(tree, data2, nuisance = FALSE)
mvgls(data2~1,tree=tree, model="BM",error=TRUE,method = "LL")

# we can compare the ancestral states or the sigma matrices
# fit1=mvbm(tree, data2)
# fit2=mvbm(tree, data2, nuisance = FALSE)
# sum((fit1$R-sigma)^2)
# sum((fit2$R-sigma)^2)

# add some (isotropic) noise to the data + known measurment error
sdval = rchisq(2*Ntip(tree), df=1)/sqrt(100) # known s.e.m.
data3 = data + matrix(rnorm(2*Ntip(tree), sd=sqrt(0.05 + sdval)), ncol=2) # I add the known (sampling variance/measurement error) and unknown variance to the noise
mvbm(tree, data3)
mvbm(tree, data3, nuisance = FALSE)
mvbm(tree, data3, nuisance = TRUE, error=matrix(sdval,ncol=2))
mvgls(data3~1,tree=tree, model="BM",error=TRUE,method = "LL")


#data4 = data + matrix(rnorm(2*Ntip(tree), sd=sqrt(sdval)), ncol=2) 
#mvBM(tree, data4, model="BM1", error=matrix(sdval,ncol=2))
#mvbm(tree, data4, nuisance = FALSE, error=matrix(sdval,ncol=2))
