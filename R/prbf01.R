prbf01 <- function(n,y,target){
  p <- beta(y+2,n-y+1)/beta(1+y,n-y+1)
  dbinom(y,n,prob=target)*exp(1)/dbinom(y,n,prob=p)
}
