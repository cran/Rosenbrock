dfull <- function(x,a,b,mu){
  n <- length(x)
  xi <- utils::head(x,-1)
  xi[1] <- xi[1] - mu
  xi1 <-x[-1]
  const <- pi^((n-1)/2) / prod(sqrt(a*b))

  return(
    const*exp( -sum( a * (xi1 - xi^2)^2 + b * (1 - xi)^2 ) )
  )
}




deven <- function(x,a,b,mu){
  n <- length(x)
  if(n %% 2 != 0){
    stop("Vector needs to be of even length")
  }else{
    x1 <- x[seq(1,n-1,2)]
    x2 <-x[seq(2,n,2)]
    const <- pi^(n/2)/prod(sqrt(a*b))

    return(
      const*exp(- sum(a * (x1 - mu)^2 - b* (x2 - x1^2)^2) )
    )
  }

}


reven <- function(n,a,b,mu){
    ns <- 2*length(a)
    mat <- matrix(NA,nrow = n,ncol = ns)

    x1 <- MASS::mvrnorm(n = n,mu = mu, Sigma = diag(1/(2*a)))
    x2 <- t(apply(X = x1,MARGIN = 1,function(t){MASS::mvrnorm(n = 1,mu = t^2, Sigma = diag(1/(2*b)))}))
    mat[,seq(1, ns, 2)] <- x1
    mat[,seq(2, ns, 2)] <- x2
    return(mat)
}



dhybrid <- function(xprime,x,a,b,mu){
  n2 <- length(x)
  n1 <- length(x[[1]])

  const <- pi^((n2*(n1-1)+1)/2)/sqrt(a)/prod(unlist(lapply(X = b,FUN = function(x){prod(sqrt(x))})))

  innerprod <- function(z,p){
    return(sum(p*(z[-1] - utils::head(z,-1)^2)^2))
  }


  return(
    exp(-(a*(xprime - mu)^2 + sum(mapply(z = x,p = b, FUN = innerprod))))
  )
}

rhybrid<-function(n,a,b,mu){
  blocks <- length(b)
  blocklengths <- unlist(lapply(X = b,FUN = length))
  mat <- matrix(NA,nrow = n,ncol = 1)

  mat[,1] <- stats::rnorm(n = n,mean = mu,sd = 1/(2*a))

  for(j in 1:blocks){
    for(i in 1:blocklengths[j]){
      if(i == 1){
        prev = mat[,1]
      }else{
        prev = mat[,ncol(mat)]
      }
      mat <- cbind(mat,MASS::mvrnorm(n = 1,mu = prev^2,diag(n)/(2*b[[j]][i])))
    }
  }

  return(mat)
}
