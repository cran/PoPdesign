## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# install.packages("PoPdesign")
library(PoPdesign)

## -----------------------------------------------------------------------------
bd <-  get.boundary.pop(n.cohort = 10, cohortsize = 3, target=0.3, 
                        cutoff=exp(1), K=3,cutoff_e=exp(-1))
summary(bd)

## -----------------------------------------------------------------------------
plot(bd)

## ----echo=FALSE, out.width="400px"--------------------------------------------
link = system.file("Flowchart", "PoP_flowchart.png", package = "PoPdesign")
knitr::include_graphics(link)


## -----------------------------------------------------------------------------
oc <- get.oc.pop(target=0.3,n.cohort=10,cohortsize=3,titration=TRUE,
                 cutoff=TRUE,cutoff_e=exp(-1),
                 skeleton=c(0.3,0.4,0.5,0.6),n.trial=1000,
                     risk.cutoff=0.8,earlyterm=TRUE,start=1)

summary(oc) # summarize design operating characteristics
plot(oc)

## -----------------------------------------------------------------------------
n <- c(4, 4, 16, 8, 0) 
y <- c(0, 0, 5, 5, 0)
selmtd <- select.mtd.pop(target = 0.2, n.pts=n, n.tox=y)
summary(selmtd)
plot(selmtd) 


