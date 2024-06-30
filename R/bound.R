bound <- function(target,n.cohort,cohortsize,cutoff,K,cutoff_e=5/24){

  lower <- upper <- 0
  lower.ex <- upper.ex <- 0
  for(i in 1:n.cohort){
    t <- i*cohortsize
    x <- 0:(i*cohortsize)
    y <- lapply(x,prbf01,n=(i*cohortsize),target=target)
    a <- 0

    for(j in 1:length(x)){
      if(y[[j]]<cutoff){ #e,exp(1)*(1/t*log(1+t))^(1/(2*t)),exp(1)*(1/t*log(1+t))^(1/(t))
        if(x[j]/(i*cohortsize)<target){
          a[j] <- 1
        }else{
          a[j] <- -1
        }
      }else{
        a[j] <- 0
      }
    }

    if (any(a==1)){
      lower[i] <- max(which(a==1))-1
    } else {
      lower[i] <- NA
    }
    if (any(a==-1)){
      upper[i] <- min(which(a==-1))-1
    } else {
      upper[i] <- NA
    }

    ex <- 0
    for(j in 1:length(x)){
      if(y[[j]]<cutoff_e){
        if(x[j]/(i*cohortsize)<target){
          ex[j] <- 1 # exclude for being subtherapeutic
        }else{
          ex[j] <- -1 # exclude for being too toxic
        }
      }else{
        ex[j] <- 0
      }
    }
    if (any(ex==1)){
      lower.ex[i] <- max(which(ex==1))-1
    } else {
      lower.ex[i] <- NA
    }
    if (any(ex==-1)){
      upper.ex[i] <- min(which(ex==-1))-1
    } else {
      upper.ex[i] <- NA
    }
  }
  b <- rbind(lower,upper,lower.ex,upper.ex)
  return(b)
}
