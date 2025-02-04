% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/std_diff.R
\name{std_diff_plot}
\alias{std_diff_plot}
\title{Compare weighted and unweighted (naive analysis) standardized difference in plot}
\usage{
std_diff_plot(
  diff.list,
  legend.pos = "right",
  prefix.title = "In strata:",
  xlim.low = 0,
  xlim.high = 1
)
}
\arguments{
\item{diff.list}{data list returned by function \link{std_diff}.}

\item{legend.pos}{legend position: "left", "top", "right", "bottom".}

\item{prefix.title}{prefix for title}

\item{xlim.low}{(\code{numeric}) lower bound of xlim}

\item{xlim.high}{(\code{numeric}) upper bound of xlim}
}
\description{
Compare weighted and unweighted (naive analysis) standardized difference in plot
}
\examples{
\dontrun{
library(dplyr)
clinical_1 <- clinical \%>\% mutate( 
  indicator = case_when(
    STRATUM == "strata_1" ~ 0, 
    STRATUM == "strata_2" ~ 1,
    is.na(STRATUM) & ARM == "experimental" ~ 1,
    TRUE ~ -1 
  ),
  ARM  = factor(ARM, levels = c("control","experimental")),
  BNLR = case_when(
    is.na(BNLR) ~ median(BNLR, na.rm = TRUE),
    TRUE ~ BNLR
  )
)
ipw_res1 <- ipw_strata(
  data.in = clinical_1, formula = indicator ~ BECOG + SEX + BNLR,
  indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", trt = "ARM",
  class.of.int = list("strata_1" = 1, "strata_2" = 0)
 )
ipw_diff <- std_diff(
  data.in = ipw_res1$data, vars = c("BECOG", "SEX", "BNLR"),
  indicator.var = "indicator", trt = "ARM",
  class.of.int = list("strata_1" = 1, "strata_2" = 0),
  usubjid.var = "SUBJID"
)
std_diff_plot(ipw_diff)
}
}
