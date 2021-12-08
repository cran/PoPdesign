#'
#' Generate descriptive summary for objects returned by other functions in PoPdesign
#'
#' Generate descriptive summary for objects returned by other functions.
#'
#' @param object the object returned by other functions.
#' @param ... ignored arguments
#'
#' @import knitr
#'
#' @return \code{summary()} prints the objects returned by other functions.
#'
#' @examples
#' ## summarize the results returned by get.boundary.pop()
#' bound <- get.boundary.pop(n.cohort = 10, cohortsize = 3, target=0.3,
#'                           cutoff=exp(1), K=3,cutoff_e=exp(-1))
#' summary(bound)
#'
#' ## summarize the results returned by get.oc.pop()
#' oc <- get.oc.pop(target=0.3,n.cohort=10,cohortsize=3,titration=TRUE,
#'                  cutoff=TRUE,cutoff_e=exp(-1),skeleton=c(0.3,0.4,0.5,0.6),n.trial=1000,
#'                  risk.cutoff=0.8,earlyterm=TRUE,start=1)
#' summary(oc)
#'
#' ### summarize the results returned by select.mtd.pop()
#' n <- c(3, 3, 15, 9, 0)
#' y <- c(0, 0, 4, 4, 0)
#' selmtd <- select.mtd.pop(target=0.3,n.pts=n, n.tox=y)
#' summary(selmtd)
#'
#' @export

summary.pop<-  function (object, ...) {

  ## get.boundary -----
  if (!is.null(object$out.boundary)) {
    cat("The decision boundaries are\n")
    print(object$out.boundary)
    cat("\n")
    cat("A more completed version of the decision boundaries is given by\n")
    print(object$out.full.boundary)
  }

  ## get.oc -----
  if (!is.null(object$sel.pct)) {
    cat("selection percentage at each dose level (%):\n")
    cat(formatC(object$sel.pct*100, digits = 1, format = "f"),
        sep = "  ", "\n")
    cat("average number of patients treated at each dose level:\n")
    cat(formatC(object$num.p, digits = 1, format = "f"),
        sep = "  ", "\n")
    cat("average number of toxicity observed at each dose level: ")
    cat(formatC(object$num.tox, digits = 1, format = "f"),
        sep = "  ", "\n")
    cat("average number of toxicities:", formatC(sum(object$num.tox),
                                                 digits = 1, format = "f"), "\n")
    cat("average number of patients:", formatC(sum(object$num.p),
                                               digits = 1, format = "f"), "\n")
    cat("percentage of early stopping due to toxicity:",
        formatC(object$early*100, digits = 1, format = "f"),
        "% \n")
    cat("risk of underdosing (>80% of patients treated below the MTD):",
        formatC(object$risk.under*100, digits = 1, format = "f"),
        "% \n")
    cat("risk of overdosing (>80% of patients treated above the MTD):",
        formatC(object$risk.over*100, digits = 1, format = "f"),
        "% \n")
  }


  ## select.mtd -----
  if (!is.null(object$MTD)) {
    tbl <- data.frame(dose = 1:length(object$p_est),
                      dlt = round(object$p_est,2))
    colnames(tbl) <- c("Dose level", "Posterior DLT")
    cat("The MTD is dose level ", object$MTD, "\n\n")
    # print(tbl)
    # kable(tbl, "simple")
    cat("Dose    Posterior DLT             \n",
        sep = "")
    cat("Level     Estimate \n", sep = "")
    for (i in 1:length(object$p_est)) {
      cat(" ", i, "        ", as.character(round(object$p_est[i],2)),
          "\n")
    }

  }



}





