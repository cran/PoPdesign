% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.oc.pop.R
\name{get.oc.pop}
\alias{get.oc.pop}
\title{Operating characteristics for single-agent trials}
\usage{
get.oc.pop(target,n,cohortsize,titration,skeleton,n.trial,cutoff,cutoff_e,
                     risk.cutoff,earlyterm,start,seed)
}
\arguments{
\item{target}{the target DLT rate}

\item{n}{total sample size}

\item{cohortsize}{the cohort size}

\item{titration}{default is TRUE. Set \code{titration=TRUE} to perform dose
escalation with cohort size = 1 to accelerate dose escalation
at the beginning of the trial.}

\item{skeleton}{a vector containing the true toxicity probabilities of the
investigational dose levels.}

\item{n.trial}{the total number of trials to be simulated}

\item{cutoff}{the cutoff for the predictive Bayes Factor (PrBF). Users can specify either a value or a function
for cutoff. If PrBF < cutoff, we assign the next cohort of patients to an adjacent dose based on observed DLT.
Otherwise, the evidence is in favor of \eqn{H_{0j}} and we need to retain the current dose.}

\item{cutoff_e}{the cutoff for the dose exclusion rule. If \eqn{PrBF_{0,1}<E(n_j)}, the evidence is in favor of \eqn{H_{1j}}. If \eqn{\hat{\pi}_j < \phi},
the current dose is deemed as subtherapeutic and we exclude the current dose and lower doses; If \eqn{\hat{\pi}_j > \phi}, the current dose
is overly toxic and we exclude the current dose and higher doses.}

\item{risk.cutoff}{the cutoff to eliminate an over/under toxic dose.
We recommend the default value of (\code{risk.cutoff=0.8}) for general use.}

\item{earlyterm}{the early termination parameter.}

\item{start}{specify the starting dose level. Default value is 1.}

\item{seed}{the seed for random number generation. Default is 123.}
}
\value{
\code{get.oc.pop()} returns the operating characteristics of the PoP design as a list,
       including:

       (1) selection percentage at each dose level (\code{$sel.pct}),

       (2) the number of patients treated at each dose level (\code{$num.p}),

       (3) the number of toxicities observed at each dose level (\code{$num.tox}),

       (4) the average number of toxicities,

       (5) the average number of patients,

       (6) the percentage of early stopping without selecting the MTD (\code{$early}),

       (7) risk of underdosing 80\% or more of patients (\code{$risk.under}),

       (8) risk of overdosing 80\% or more of patients (\code{$risk.over})
}
\description{
Generate the operating characteristics of the PoP design by simulating trials.
}
\details{
TBD
}
\examples{

## get the operating characteristics for single-agent trials
oc <- get.oc.pop(target=0.3,n=15,cohortsize=3,titration=TRUE,
                 cutoff=2.5,cutoff_e=5/24,
                 skeleton=c(0.3,0.4,0.5,0.6),n.trial=1000,
                     risk.cutoff=0.8,earlyterm=TRUE,start=1, seed=123)

summary(oc) # summarize design operating characteristics
plot(oc)

}
\references{
Brunk, H., Barlow, R. E., Bartholomew, D. J. & Bremner, J. M (1972, ISBN-13: 978-0471049708).
}
