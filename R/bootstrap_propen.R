#' Calculate bootstrap CI for treatment effect estimate
#' @inheritParams ipw_strata
#' @inheritParams ps_match_strata
#' @param seed seed
#' @param n.boot number of bootstraps to run; note only runs without warning/error msg will be considered
#' when calculating bootstrap CI; so it is possible more that n.boot runs are performed (but capped by max.num.run)
#' @param max.num.run max number of bootstraps to run (include both valid and not valid runs)
#' @param estimate.res result object from ipw_strata() or ps_match_strata()
#' @param method "ipw" or "ps". If "ipw", ipw_strata() will be used. If "ps", ps_match_strata() will be used.
#' @param hr.ratio.ref no longer to be used, please use pairs instead
#' @param ref.denom no longer to be used, please use pairs instead
#' @param pairs pairs of interest when calculating ratio of HR (delta of delta for OR).
#' this should be a matrix whose rows are names of strata, 1st column indicates the stratum to be used as numerator (HR or ORR diff);
#' 2nd column indicates denominator.
#' If pairs is NULL, ratio of HR (difference of OR difference) will not be calculated.
#' @param wild.boot whether wild bootstrap should be used. If so, weights will be generated using rexp(1)
#' @param non.converge.check whether to output number of time each level of each categorical variable for each stratum specified in indicator having N<non.converge.check.thresh when non-convergence occurs
#' @param non.converge.check.thresh see above
#' @note only estimates from runs without error or warning will be considered when calculating bootstrap CI.
#' If none of the bootstrap runs is error/warning free, CI of est.ci.mat will be NA
#' 
#' @importFrom stats rexp qnorm var median
#' @importFrom dplyr `%>%`
#' @export
#' @examples 
#' library(dplyr)
#' clinical_1 <- clinical %>% mutate( 
#'   indicator = case_when(
#'     STRATUM == "strata_1" ~ 0, 
#'     STRATUM == "strata_2" ~ 1,
#'     is.na(STRATUM) & ARM == "experimental" ~ 1,
#'     TRUE ~ -1 
#'   ),
#'   ARM  = factor(ARM, levels = c("control","experimental")),
#'   BNLR = case_when(
#'     is.na(BNLR) ~ median(BNLR, na.rm = TRUE),
#'     TRUE ~ BNLR
#'   )
#' )
#' # ipw: default model
#' ipw_res <- ipw_strata(
#'   data.in = clinical_1, formula = indicator ~ BECOG + SEX + BNLR,
#'   indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", trt = "ARM",
#'   class.of.int = list("strata_1" = 1, "strata_2" = 0)
#'  )
#' boot_ipw <- bootstrap_propen(
#'   data.in = clinical_1, formula = indicator ~ BECOG + SEX + BNLR,
#'   indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", trt = "ARM",
#'   class.of.int = list("strata_1" = 1, "strata_2" = 0),
#'   estimate.res = ipw_res, method = "ipw", n.boot = 5
#' )
#' boot_ipw$est.ci.mat
#' boot_ipw$boot.out.est 
#' # ps: DWC model
#' clinical_2 <- clinical %>% mutate( 
#'   indicator = case_when(
#'     STRATUM == "strata_1" ~ 0, 
#'     STRATUM == "strata_2" ~ 1,
#'     is.na(STRATUM) & ARM == "experimental" ~ 2,
#'     TRUE ~ -1
#'   ),
#'   ARM  = factor(ARM, levels = c("control","experimental")),
#'   BNLR = case_when(
#'     is.na(BNLR) ~ median(BNLR, na.rm = TRUE),
#'     TRUE ~ BNLR
#'   )
#' )
#' ps_res <- ps_match_strata(
#'   data.in = clinical_2, formula = indicator ~ BECOG + SEX + BNLR, model = "dwc",
#'   indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", trt = "ARM",
#'   class.of.int = list("strata_1" = 0, "strata_2" = 1, "missing" = 2)
#'  ) 
#' boot_ps <- bootstrap_propen(
#'   data.in = clinical_2, formula = indicator ~ BECOG + SEX + BNLR, model = "dwc",
#'   indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", trt = "ARM",
#'   class.of.int = list("strata_1" = 0, "strata_2" = 1, "missing" = 2),
#'   estimate.res = ps_res, method = "ps", n.boot = 5
#' )
#' boot_ps$est.ci.mat
#' boot_ps$boot.out.est 
#' 
#' @return return a `list` containing the following components:
#' \itemize{
#' \item boot.out.est a `matrix` with rows as estimates such as Coefficient and Variance in strata 
#'  and columns as summary statistics such as Mean and Median of the estimates.
#' \item est.ci.mat a `matrix` with rows as strata and columns as Estimate and Bootstrap CIs.
#' \item eff.comp.ci.mat a `matrix` with rows as strata comparisons and columns as Estimate and Bootstrp CIs.
#' \item conv.est a `logical` vector to indicate whether model in each bootstrap converges.
#' \item error.est `numeric` to indicate the total number of models in bootstrap which gives errors.
#' \item boot.results a `matrix` with rows as each bootstrap and columns as model results such as Coefficient in strata.
#' \item glm.warn.est a `logical` vector to indicate whether glm model gives warning in each bootstrap.
#' \item num.valid.boots `numeric` to indicate the total number of valid bootstraps.
#' \item num.total.boots `numeric` for the total number of bootstrap runs.
#' \item warning.msg a `list` to capture warning messages from models.
#' \item error.msg a `list` to capture error messages from models.
#' \item non.converge.dat a `matrix` 
#'  with rows as each level of each categorical variable for each stratum specified in 
#'  indicator having N less than `non.converge.check.thresh` and columns as treatment groups
#' }
#' 


bootstrap_propen <- function(data.in, indicator.var = "indicator",
                             formula, indicator.next = NULL,
                             seed = 100, class.of.int, estimate.res,
                             n.boot = 1000, method = "ipw", wild.boot = FALSE,
                             tte = "AVAL", event = "event", trt = "trt",
                             response = NULL, caliper = NULL, pairs = NULL, hr.ratio.ref = NULL, ref.denom = TRUE,
                             model = "plain", max.num.run = 5000, non.converge.check = FALSE, multinom.maxit = 100,
                             non.converge.check.thresh = 1) {
  set.seed(seed)

  # bootstrap
  conv.est <- stat.est <- warn.est <- NULL
  error.est <- 0
  ci.term <- "Coefficient"
  if (!is.null(response)) ci.term <- "Diff"
  warning.msg <- estimate.res$warning.msg
  for (j in warning.msg) j <- NULL
  error.msg <- NULL

  non.converge.dat <- NULL

  boot.names <- paste(rep(names(class.of.int), each = ncol(estimate.res$stat)),
    rep(colnames(estimate.res$stat), length(class.of.int)),
    sep = ";"
  )
  i <- i.run <- 0

  if (non.converge.check) {
    formula.vars <- (formula %>% as.character())[3] %>%
      strsplit(split = " \\+ ") %>%
      unlist()
    tmp <- SummaryVars(data.in, var = formula.vars, trt = indicator.var)
    non.converge.dat <- matrix(0,
      nrow = nrow(tmp), ncol = ncol(tmp),
      dimnames = list(row.names(tmp), colnames(tmp))
    )
    non.converge.dat <- non.converge.dat[rownames(non.converge.dat) != "Total (non-NA)", ]
    rownames(non.converge.dat)[rownames(non.converge.dat) == "NA's"] <- "Non-converge"
  }

  while (i < n.boot & i.run < max.num.run) {
    i <- i + 1
    i.run <- i.run + 1
    if (!wild.boot) {
      ns <- tapply(1:nrow(data.in), data.in[, indicator.var], c)
      id <- sapply(ns, function(j) sample(j, length(j), replace = TRUE))
      data.boot <- data.in[unlist(id), ]
      weights <- NULL
    } else {
      data.boot <- data.in
      weights <- rexp(nrow(data.in), rate = 1)
    }
    if (method == "ipw") {
      est.boot <- try(ipw_strata(
        data.in = data.boot, formula = formula, trt = trt,
        class.of.int = class.of.int, tte = tte, event = event,
        response = response, model = model,
        indicator.var = indicator.var,
        indicator.next = indicator.next,
        weights = weights, multinom.maxit = multinom.maxit
      ), silent = TRUE)
    }
    if (method == "ps") {
      est.boot <- try(ps_match_strata(
        data.in = data.boot, formula = formula, trt = trt,
        class.of.int = class.of.int, tte = tte, event = event, response = response,
        model = model, indicator.var = indicator.var, caliper = caliper,
        weights = weights, multinom.maxit = multinom.maxit
      ), silent = TRUE)
    }
    if (class(est.boot) == "try-error" || any(is.na(est.boot$stat[, 3]))) {
      error.est <- error.est + 1
      error.msg <- c(error.msg, est.boot[[1]])
      i <- i - 1
    } else {
      conv.est <- c(conv.est, est.boot$converged)
      if (!est.boot$converged) i <- i - 1
      if (!est.boot$converged & non.converge.check) {
        # browser()
        tmp0 <- SummaryVars(est.boot$data, var = formula.vars, trt = indicator.var)
        tmp0 <- tmp0[rownames(tmp0) != "Total (non-NA)", ]
        rownames(tmp0)[rownames(tmp0) == "NA's"] <- "Non-converge"
        tmp0[rownames(tmp0) == "Non-converge"] <- 0

        tmp <- sapply(tmp0, function(x) (unlist(strsplit(x, split = " "))[1] %>% as.numeric()) < non.converge.check.thresh) %>% matrix(nrow = nrow(tmp0), byrow = FALSE)

        non.converge.dat <- non.converge.dat + tmp
      }

      if (!is.null(unlist(est.boot$warning.msg))) {
        if (est.boot$converged) i <- i - 1 # if already reset counter by convergence criteria, dont reset count again
        for (j in 1:length(warning.msg)) if (!is.null(est.boot$warning.msg[[j]])) warning.msg[[j]] <- c(warning.msg[[j]], est.boot$warning.msg[[j]])
        warn.est <- c(warn.est, est.boot$any_warning_glm)
      } else {
        stat.est <- rbind(stat.est, as.vector(t(est.boot$stat)))
      }
    }
  }
  if (is.null(stat.est)) stat.est <- matrix(NA, ncol = length(boot.names), nrow = 1)
  colnames(stat.est) <- boot.names
  num.valid <- i
  if (i.run == max.num.run) {
    warning(paste(
      "Among", max.num.run, "bootstraps, only",
      i, "are valid"
    ))
  }
  if (is.null(hr.ratio.ref)) hr.ratio.ref <- names(class.of.int) [1] # if RHR reference stratum is not defined, use the 1st level
  if (is.null(pairs)) { # make sure updated codes with param "pairs" can wont break old wrappers
    hr.ratio.nonref <- names(class.of.int) # setdiff(names(class.of.int), hr.ratio.ref)
    hr.ratio.nonref <- matrix(hr.ratio.nonref, ncol = 1)
    pairs <- cbind(hr.ratio.nonref, hr.ratio.ref)
    if (!ref.denom) pairs <- cbind(pairs[, 2], pairs[, 1]) # if switch numerator and denominator
  }


  # Ratio of HR or diff of ORR diff
  if (!is.null(pairs)) { # AN
    # stopifnot(hr.ratio.ref%in%names(class.of.int))
    stopifnot(all(pairs %in% names(class.of.int))) # AN
    eff.name <- ifelse(is.null(response), "HR", "Diff") # AN
    # if(!ref.denom) pairs <- cbind(pairs[,2],pairs[,1]) # AN
    eff.terms <- boot.names[which(sapply(boot.names, function(k) strsplit(k, split = ";")[[1]][2]) == eff.name)]
    eff.ref <- paste0(pairs[, 2], ";", eff.name) # paste0(hr.ratio.ref,";HR")
    eff.noref <- paste0(pairs[, 1], ";", eff.name) # hr.noref <- setdiff(hr.terms, hr.ref)
    # eff.ref   <- pairs[,2]
    # eff.noref <- pairs[,1]
    # bootstrap
    # if(ref.denom)hr.ratio.calc <- sapply(hr.noref, function(k)stat.est[,k]/stat.est[,hr.ref], simplify = FALSE)
    # if(!ref.denom)hr.ratio.calc <- sapply(hr.noref, function(k)stat.est[,hr.ref]/stat.est[,k], simplify = FALSE)
    if (eff.name == "HR") eff.comp.calc <- sapply(1:nrow(pairs), function(k) stat.est[, eff.noref[k]] / stat.est[, eff.ref[k]], simplify = FALSE)
    if (eff.name == "Diff") eff.comp.calc <- sapply(1:nrow(pairs), function(k) stat.est[, eff.noref[k]] - stat.est[, eff.ref[k]], simplify = FALSE)
    eff.comp.mat <- do.call(cbind, eff.comp.calc)
    # if(ref.denom)colnames(hr.ratio.mat) <-  paste(sapply(hr.noref, function(k)strsplit(k, split=";")[[1]][1]), "over", hr.ratio.ref)
    # if(!ref.denom)colnames(hr.ratio.mat) <- paste( hr.ratio.ref, "over", sapply(hr.noref, function(k)strsplit(k, split=";")[[1]][1]))
    # colnames(eff.comp.mat) <- apply(pairs, 1, function(x) paste0(x[1], ';', eff.name, ' over ', x[2],';', eff.name))
    colnames(eff.comp.mat) <- apply(pairs, 1, function(x) paste0(x[1], " over ", x[2]))
    # NL: need to keep "over" here since Sautoir wrapper looks for it
    stat.est <- cbind(stat.est, eff.comp.mat)

    # point estimate, calculate every one
    est.eff <- estimate.res$stat[, eff.name]
    # point estimate
    # if(ref.denom)hr.ratio.point <- sapply(names(est.hr), function(k)est.hr[k]/est.hr[hr.ratio.ref])
    # if(!ref.denom)hr.ratio.point <- sapply(names(est.hr), function(k)est.hr[hr.ratio.ref]/est.hr[k])
    if (eff.name == "HR") {
      eff.comp.point <- sapply(1:nrow(pairs), function(k) est.eff[pairs[k, 1]] / est.eff[pairs[k, 2]])
    } # AN
    if (eff.name == "Diff") {
      eff.comp.point <- sapply(1:nrow(pairs), function(k) est.eff[pairs[k, 1]] - est.eff[pairs[k, 2]])
    } # AN
    names(eff.comp.point) <- colnames(eff.comp.mat)


    boot.out.est <- t(apply(stat.est, 2, function(i) {
      c(
        mean(i), var(i), median(i), range(i),
        var(suppressWarnings(log(i)))
      )
    })) # expect to see warnings since some are always 0 by nature
    colnames(boot.out.est) <- c("Mean", "Var", "Median", "Min", "Max", "VarLog")


    # if(model=="wri"){
    #  for(cc in 1:length(class.of.int)) class.of.int[[cc]] <- setdiff(class.of.int[[cc]],2)
    #  check.class <- sapply(class.of.int, length)
    #  class.of.int <- class.of.int[check.class>0] # cannot produce any class where level 2 is included
    # }
    est.ci.mat <- NULL
    for (i in names(class.of.int)) {
      tmp <- c(
        estimate.res$stat[i, ci.term],
        estimate.res$stat[i, ci.term] + c(-1, 1) * qnorm(.975) * sqrt(boot.out.est[paste0(i, ";", ci.term), "Var"])
      )
      est.ci.mat <- rbind(est.ci.mat, tmp)
    }
    # browser()
    rownames(est.ci.mat) <- names(class.of.int)
    colnames(est.ci.mat) <- c("Estimate", "Bootstrap CI low", "Bootstrap CI high")

    # if(!is.null(hr.ratio.ref)){
    #   est.ci.mat <- cbind(est.ci.mat, hr.ratio.point)
    #   if(ref.denom) colnames(est.ci.mat)[4] <- paste0("HR ratio: group over", hr.ratio.ref) else
    #     colnames(est.ci.mat)[4] <- paste0("HR ratio: ", hr.ratio.ref, "over group")
    # }

    if (!is.null(pairs)) { # AN
      eff.comp.ci.mat <- eff.comp.point
      eff.comp.ci.mat <- cbind(eff.comp.ci.mat, t(sapply(names(eff.comp.point), function(x) {
        if (eff.name == "Diff") tmp.re <- eff.comp.point[x] + c(-1, 1) * qnorm(.975) * sqrt(boot.out.est[x, "Var"])
        if (eff.name == "HR") tmp.re <- exp(log(eff.comp.point[x]) + c(-1, 1) * qnorm(.975) * sqrt(boot.out.est[x, "VarLog"]))
        tmp.re
      })))
    }
    eff.comp.ci.mat <- matrix(eff.comp.ci.mat, ncol = 3)
    rownames(eff.comp.ci.mat) <- names(eff.comp.point)
    colnames(eff.comp.ci.mat) <- c("Estimate", "Bootstrap CI low", "Bootstrap CI high")

    if (!is.null(hr.ratio.ref)) { # for backward compatibility
      est.ci.mat <- cbind(est.ci.mat, eff.comp.point)
      if (ref.denom) {
        colnames(est.ci.mat)[4] <- paste0("HR ratio: strata over ", hr.ratio.ref)
      } else {
        colnames(est.ci.mat)[4] <- paste0("HR ratio: ", hr.ratio.ref, "over strata")
      }
    }
  }

  out <- list(
    boot.out.est = boot.out.est, est.ci.mat = est.ci.mat, eff.comp.ci.mat = eff.comp.ci.mat, # AN
    conv.est = conv.est, error.est = error.est, boot.results = stat.est,
    glm.warn.est = warn.est, num.valid.boots = num.valid, num.total.boots = i.run,
    warning.msg = warning.msg, error.msg = error.msg, non.converge.dat = non.converge.dat
  )
}
