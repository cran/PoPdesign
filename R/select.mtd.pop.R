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

  fit.isoreg <- function(iso, x0)
  {
    if(length(x0)==1){
      return(iso$yf)
    }
    o = iso$o
    if (is.null(o))
      o = 1:length(x0)
    x = unique(iso$x[o])
    y = iso$yf
    ind = cut(x0, breaks = x, labels = FALSE, include.lowest = TRUE)
    min.x <- min(x)
    max.x <- max(x)
    adjusted.knots <- iso$iKnots[c(which(iso$yf[iso$iKnots] > 0))]
    fits = sapply(seq(along = x0), function(i) {
      j = ind[i]

      # Find the upper and lower parts of the step
      upper.step.n <- min(which(adjusted.knots > j))
      upper.step <- adjusted.knots[upper.step.n]
      lower.step <- ifelse(upper.step.n==1, 1, adjusted.knots[upper.step.n -1] )

      # Perform a liner interpolation between the start and end of the step
      denom <- x[upper.step] - x[lower.step]
      denom <- ifelse(denom == 0, 1, denom)
      val <- y[lower.step] + (y[upper.step] - y[lower.step]) * (x0[i] - x[lower.step]) / (denom)
    })
    fits
  }

  l <- which(n.pts>0)
  p <- n.tox[l]/n.pts[l]
  if(sum(p)==0){
    return(max(l))
  }
  iso.model <- isoreg(p)
  p.iso <- fit.isoreg(iso.model,1:length(l))
  d <- abs(p.iso-target)
  mtd <- l[max(which(d==min(d)))]
  p_est <- p.iso
  out <- list(target = target,
              MTD = mtd,
              p_est = p_est)
  class(out)<-"pop"
  return(out)
}
