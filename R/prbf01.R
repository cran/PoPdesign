prbf01 <- function(n,y,target){
  p <- exp(lbeta(y+2,n-y+1)-lbeta(1+y,n-y+1))
  dbinom(y,n,prob=target)*exp(1)/dbinom(y,n,prob=p)
}
