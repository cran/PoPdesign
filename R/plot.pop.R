#'
#' Plot the flowchart and simulation results for PoP designs
#'
#' Plot the objects returned by other functions, including (1) flowchart of PoP design;
#' (2) operating characteristics of the design, including selection percentage and the
#' number of patients treated at each dose;
#' (3) the estimate of toxicity probability for each dose and corresponding 95\% credible interval
#'
#' @param x the object returned by other functions
#' @param ... ignored arguments
#'
#'
#' @return \code{plot()} returns a figure or a series of figures depending on the object entered
#'
#'
#' @importFrom graphics abline barplot axis
#' @importFrom magick image_read
#' @export

plot.pop<- function (x, ...){
  ## get.boundary -----
  if (!is.null(x$out.boundary)) {
    link = system.file("Flowchart", "PoP_flowchart.png", package = "PoPdesign")
    flowchart <- image_read(link)
    print(flowchart)

    tbl <- x$out.boundary
    print(tbl)
    }

  ## get.oc.pdf -----
  if (!is.null(x$sel.pct)) {
    dt.plot <- data.frame(dose = 1:length(x$num.p),
                          pct = x$sel.pct[1:length(x$num.p)])
    p <- barplot(dt.plot$pct,names.arg=as.character(1:length(x$num.p)),
                 xlab="Dose level", ylab="Selection percentage (%)")
    p
  }


  ## select.mtd.pdf -----
  if (!is.null(x$MTD)) {
    plot(x$p_est, ylab="DLT rate", xlab="Dose level",pch=20,xaxt = "n",
         col=ifelse(x$p_est==x$p_est[x$MTD], "red", "black"))
    axis(1, at=1:length(x$p_est))
    abline(h=x$target, col="Red", lty=3)
  }
}


