#'
#' Generate the dose escalation and de-escalation boundaries for single-agent trials.
#'
#' Use this function to generate the dose escalation and deescalation boundaries for single-agent trials.
#'
#' @usage get.boundary.pop(target, n, cohortsize, cutoff, K, cutoff_e)
#'
#' @param target the target DLT rate
#' @param n total sample size
#' @param cohortsize the cohort size
#' @param cutoff the cutoff for the predictive Bayes Factor (PrBF). Users can specify either a value or a function
#'               for cutoff. If PrBF < cutoff, we assign the next cohort of patients to an adjacent dose based on observed DLT.
#'               Otherwise, the evidence is in favor of \eqn{H_{0j}} and we need to retain the current dose.
#' @param K number of dose levels. It is required when argument cutoff is a function that requires K.
#' @param cutoff_e the cutoff for the dose exclusion rule. If \eqn{PrBF_{0,1}<E(n_j)}, the evidence is in favor of \eqn{H_{1j}}. If \eqn{\hat{\pi}_j < \phi},
#'                 the current dose is deemed as subtherapeutic and we exclude the current dose and lower doses; If \eqn{\hat{\pi}_j > \phi}, the current dose
#'                 is overly toxic and we exclude the current dose and higher doses.
#'
#' @details We assume that there are \eqn{J} pre-specified dose levels of the drug of interest. Let \eqn{d_1,d_2,\ldots,d_J} denote
#'          these dose levels. The dose-limiting toxicity (DLT) is assessed as a binary outcome, experiencing
#'          toxicity or not. The true dose toxicity is monotonically increasing as the dose level increases.
#'          Let \eqn{\phi} be the target toxicity rate and \eqn{\pi_j} be the true dose-toxicity of dose level \eqn{d_j}, for \eqn{j=1,2,\ldots,J}.
#'
#'          We formulate our hypothesis as:
#'                                 \deqn{H_{0j}: \pi_j=\phi}
#'                                 \deqn{H_{1j}: \pi_j\ne\phi}
#'
#'          \eqn{H_{0j}} indicates that \eqn{d_j} is the desired MTD so that we should stay; \eqn{H_{1j}} reflects the current dose is
#'          either below or above the MTD so that we should transit to a lower or upper dose level. Whether
#'          escalate or de-escalate the dose is straightforward: if the observed toxicity rate is above the target
#'          toxicity rate \eqn{\phi}, we de-escalate the dose; if the observed toxicity rate is below \eqn{\phi}, we escalate the dose.
#'
#'          With the hypothesis, the predictive Bayes factor comparing \eqn{H_{0j}} and \eqn{H_{1j}} is given by
#'                                \deqn{PrBF_{0,1}=\frac{\phi^{y_i}(1-\phi)^{n_j-y_j}B(y_j+1,n_j-y_j+1)^{n_j}exp(1)}{B(y_j+2,n_j-y_j+1)^{y_j}B(y_j+1,n_j-y_j+2)^{n_j-y_j}}}
#'          where \eqn{x_j} is the toxicity response of the ith subject among \eqn{n_j} subjects that received dose \eqn{d_j}, for \eqn{j=1,2,\ldots,J}.
#'          \eqn{y_j} denotes the sum of toxicity response. We assume that
#'                                \deqn{y_j \sim Bin(n_j,\pi_j)}
#'
#'          According to the calibration of the PrBF, a decision rule based on \eqn{PrBF_{0,1}} is:
#'          1. If \eqn{PrBF_{0,1}>C(n_j)}, the evidence is in favor of \eqn{H_{0j}} and we need to retain the current dose;
#'          2. Otherwise, we assign the next cohort of patients to an adjacent dose according to the observed
#'          DLT \eqn{\hat{\pi}_j = y_j/n_j}, such as:
#'
#'          (a) If \eqn{\hat{\pi}_j < \phi}, we escalate the dose;
#'
#'          (b) If \eqn{\hat{\pi}_j > \phi}, we de-escalate the dose.
#'
#'          For patient safety and trial efficiency, the PoP design employs a dose exclusion rule. On
#'          the one hand, if the PrBF based on the observed DLT indicates a dose is above the MTD with a
#'          certain evidence, we exclude the current dose and doses above to avoid treating patients at an overly
#'          toxic dose; on the other hand, if the PrBF implies that a dose is substantially below the MTD, we
#'          eliminate the current dose and doses below to prevent wasting patients at a subtherapeutic dose.
#'          Such a dose exclusion rule is as follow:
#'
#'          If \eqn{PrBF_{0,1}<E(n_j)}, the evidence is in favor of \eqn{H_{1j}} and:
#'
#'          1. If \eqn{\hat{\pi}_j < \phi}, the current dose is deemed as subtherapeutic and we exclude the current dose and lower doses;
#'
#'          2. If \eqn{\hat{\pi}_j > \phi}, the current dose is overly toxic and we exclude the current dose and higher doses.
#'
#'          Once all the doses are eliminated from further investigation, the trial is terminated early.
#'          The selection of the cut-off value for the dose exclusion is critical for the performance of the PoP
#'          design, because it ensure the safety of the patients and efficiency of the design by influencing
#'          the early termination rule. The exclusion boundaries in the table above were determined using
#'          \eqn{E(n_j)=exp(-1)}.
#'
#'
#' @return  \code{get.boundary.pop()} returns a list object, including the corresponding decision tables
#'          \code{$out.boundary} and \code{$out.full.boundary}.
#'
#' @examples
#'
#' ## get the dose escalation and deescalation boundaries for PoP design with
#' ## the target DLT rate of 0.3, maximum sample size of 30, and cohort size of 3
#' bound <- get.boundary.pop(target=0.5, n = 15, cohortsize = 3,
#'                           cutoff=2.5,K=4,cutoff_e=5/24)
#' summary(bound) # get the descriptive summary of the boundary
#' plot(bound)    # plot the flowchart of the design along with decision boundaries
#'
#' @import stats
#' @export


get.boundary.pop = function(target,n,cohortsize,cutoff=2.5,K=4,cutoff_e=5/24){

  if (target < 0.05) {
    stop("the target is too low! ")
  }
  if (target > 0.6) {
    stop("the target is too high!")
  }

  n.cohort <- n/cohortsize
  if (n.cohort%%1){
    stop("The total sample size is not a multiple of cohort size!")
  }
  out.boundary <- rbind(cohortsize * (1:n.cohort),
                        bound(target=target,n.cohort=n.cohort,cohortsize = cohortsize,cutoff=cutoff,K=K,cutoff_e=cutoff_e))
  out.full.boundary <- rbind(1:n,
                             bound(target=target,n.cohort=n,cohortsize = 1,cutoff=cutoff,K=K,cutoff_e=cutoff_e))

  rownames(out.boundary) <- c("Number of patients treated",
                              "Escalation if # of DLT (U1) <=","Deescalation if # of DLT (U2) >=",
                              "Subtherapeutic exclusion if # of DLT (V1) <=",
                              "Overly toxic exclusion if # of DLT (V2) >=")

  rownames(out.full.boundary) <- c("Number of patients treated",
                                   "Escalation if # of DLT (U1) <=","Deescalation if # of DLT (U2) >=",
                                   "Subtherapeutic exclusion if # of DLT (V1) <=",
                                   "Overly toxic exclusion if # of DLT (V2) >=")

  out = list(out.boundary = out.boundary,
             out.full.boundary = out.full.boundary)
  class(out)<-"pop"
  return(out)
}


