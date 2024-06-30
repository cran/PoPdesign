#'
#' Maximum tolerated dose (MTD) selection for single-agent trials
#'
#' Select the maximum tolerated dose (MTD) when the single-agent trial is completed
#'
#' @usage select.mtd.pop(target, n.pts, n.tox)
#'
#' @param target the target DLT rate
#' @param n.pts a vector containing the number of patients treated at each dose level
#' @param n.tox a vector containing the number of patients who experienced dose-limiting
#'              toxicity at each dose level
#'
#'
#' @return  \code{select.mtd.pop()} returns (1) selected MTD (\code{$MTD}),
#' (2) isotonic estimate of the DLT probablity at each dose and associated
#'
#' @references Brunk, H., Barlow, R. E., Bartholomew, D. J. & Bremner, J. M (1972, ISBN-13: 978-0471049708).
#'
#' @examples
#'
#' ### select the MTD for PoP trial
#' n <- c(4, 4, 16, 8, 0)
#' y <- c(0, 0, 5, 5, 0)
#' selmtd <- select.mtd.pop(target=0.3,n.pts=n, n.tox=y)
#' summary(selmtd)
#' plot(selmtd)
#'
#' @export


select.mtd.pop <- function(target, n.pts, n.tox){

  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1)
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol)))
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }

  l <- which(n.pts>0)
  p <- (n.tox[l]+0.05)/(n.pts[l]+0.1)
  p.var = (n.tox[l] + 0.05) * (n.pts[l] - n.tox[l] + 0.05)/((n.pts[l] +
                                                               0.1)^2 * (n.pts[l] + 0.1 + 1))
  p.iso <- pava(p, wt = 1/p.var)
  p.iso = p.iso + (1:length(p.iso)) * 1e-10

  ## eliminate dose based on posterior probability
  K = length(n.pts)
  elimi = rep(0, K)
  for (i in 1:K) {
    if (1 - pbeta(target, n.tox[i] + 1, n.pts[i] - n.tox[i] + 1) > 0.95) {
      elimi[i:K] = 1
      break
    }
  }
  m <- which(elimi!=1)

  l <- l[m]
  p.iso <- p.iso[m]
  if(length(l)==0) {return(99)}
  # l[sort(abs(p.iso - target), index.return = T)$ix[1]]

  d <- abs(p.iso-target)
  mtd <- l[max(which(d==min(d)))]
  p_est <- p.iso

  out <- list(target = target,
              MTD = mtd,
              p_est = p_est)
  class(out)<-"pop"
  return(out)
}
