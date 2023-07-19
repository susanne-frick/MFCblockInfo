#likelihood for Thurstonian IRT, assuming local independence
posterior.tirt <- function (theta, responsesj, s.prior, m.prior, design.mat, comp.mat.a, int, uni, est="MAP") {
  prob <- pnorm((comp.mat.a%*%theta + as.vector(design.mat %*% int))/sqrt(as.vector(design.mat^2 %*% uni)))
  rprobs <- prob^responsesj*(1-prob)^(1-responsesj)
  rprobs[rprobs<=0] <- 1e-100 #number close to zero
  loglik <- sum(log(rprobs))
  if(est=="MAP") {
    post <- loglik + mvtnorm::dmvnorm(x=theta, mean=m.prior, sigma=s.prior, log = TRUE)
  } else {
    post <- loglik
  }
  return(post)
}

#MAP estimator (for local dependence)
# function to estimate MAP scores
est.MAP.old <-  function(FUN, responses, s.prior, ...) {

  nt <- ncol(s.prior) # number of traits
  J <- nrow(responses) #number of persons

  results <- matrix(nrow=J,ncol=nt,dimnames = list(rep(1:J), c("N", "E", "O", "A", "C")))
  ses <- matrix(nrow=J,ncol=nt,dimnames = list(rep(1:J), c("N", "E", "O", "A", "C")))
  errors <- NULL
  warns <- NULL
  messages <- NULL
  for (j in 1:J) {
    responsesj <- responses[j,]
    tryCatch({
      MAP <- optim(par=rep(0,length=nt),fn=FUN, responsesj=responsesj, s.prior=s.prior,
      #             method="BFGS", control=list(fnscale=-1,trace=10),hessian=T, ...)
      method="L-BFGS-B", lower=rep(-3,nt), upper=rep(3,nt), control=list(fnscale=-1),hessian=T, ...) #to create errors
      results[j,] <- MAP$par
      ses[j,] <- sqrt(diag(MASS::ginv(-MAP$hessian)))
      messages <- c(messages, MAP$message)
    }, error=function(e){
      errors <<- c(errors, conditionMessage(e))
      results[j,] <- NA
    }, warning=function(w)
      warns=c(warns, conditionMessage(w)))
  }
  return(list("traits"=results, "ses"=ses, "errors"=errors, "warns"=warns, "messages"=messages))
}

est.hessian.old <-  function(FUN, traits, responses, s.prior, ...) {
  
  nt <- ncol(s.prior) # number of traits
  J <- nrow(responses) #number of persons
  
  ses <- matrix(nrow=J,ncol=nt,dimnames = list(rep(1:J), c("N", "E", "O", "A", "C")))
  infos <- array(NA, dim = c(J, nt, nt))
  errors <- NULL
  warns <- NULL

  for (j in 1:J) {
    responsesj <- responses[j,]
    tryCatch({
      result <- optimHess(par=traits[j,],fn=FUN, responsesj=responsesj, s.prior=s.prior, ...)
      infos[j,,] <- -result
      ses[j,] <- sqrt(diag(MASS::ginv(-result)))
    }, error=function(e){
      errors <<- c(errors, conditionMessage(e))
      ses[j,] <- NA
    }, warning=function(w)
      warns=c(warns, conditionMessage(w)))
  }
  return(list("traits"=traits, "ses"=ses, "errors"=errors, "warns"=warns, "infos"=infos))
}

#true and empirical reliability
rel.emp <- function(scores, ses) {
  svar <- apply(scores, 2, var)
  svar/(svar + colMeans(ses^2, na.rm=T))
}

rel.true <- function(scores, traits) diag(cor(traits,scores))^2
