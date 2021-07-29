#' Propensity Score Matching of strata (two or more classes, survival or binary endpoint)
#'
#' This function perfroms propensity score matching of two or more strata. \cr
#' It could be used when arm1 has 2 or more strata, while strata information is unknown in arm2. \cr
#' The function will fit a logistic regression (when 2 classes) or multinomial logistic
#' regression (when > 2 classes) based on strata labels in arm1 (model: label~features), then
#' predict strata labels in both arm1 and arm2 based on the fitted model.
#' The predicted probability of being stratum X will be used for propensity score matching.
#' The matching results will then be used to estimate treatment difference of two arms (Hazard ratio
#' for survival endpoint; response rate difference for binary endpoint). When ties are allowed,
#' weights from the ties will be used to calculate the HR or response rate difference.
#' 
#' @md
#' @param ties (`logical`) `TRUE` allows for ties. \cr
#' \describe{
#'   \item{ties}{is `TRUE`, all samples in the tie will be included. When calculating summary statistics, samples in
#' ties will be assigned a smaller weight. for example, if two samples ties, these two samples
#' will both be included in the summary statistics calculation with weight 0.5.}
#'  \item{ties}{is `FALSE`, one random sample will be draw from the tied samples when calculating summary
#' statistics. In this case, it is recommended to run `ps_match_strata` multiple times with
#' different seeds and take the average or median summary statistics from multiple runs.
#' Note when ties is `FALSE`, codes were tested less thoroughly and extra caution may be needed.}
#' }
#' @inheritParams ipw_strata
#' @inheritParams Matching::Match
#' @note Different from the original version, iter is no longer a parameter
#' if tie = `FALSE` is specified, user need to run for loops of sapply outside of this function
#' to get results from multiple seeds.
#' Three elements in the output list - the data element is a data frame that contains input data and
#' estimated probabilities. The stat element contains estimated treatment difference between 2 arms, in each of the strata of interest.
#' The converge element indicates whether the model converged (taking from $converged from glm and $convergency from multinom)
#' if return.data is `FALSE`, data won't be returned.
#' model = "wri" is not supported in ps_match_strata
#' @return return a `list` containing the following components:
#' \itemize{
#' \item stat a `matrix` with rows as strata and columns as Estimate and CIs.
#' \item converged `logical` to indicate whether model converges.
#' \item any_warning_glm `logical` to indicate whether there's warning from glm model.
#' \item warning.msg a `list` to capture any warning message from the modeling process.
#' \item models a `list` to capture the glm model results.
#' \item roc.list a `list` to capture information about Area under the curve from glm model.
#' \item data a `data.frame` which is the original input data plus predicted probabilities. 
#' }
#' @examples  
#' library(dplyr)
#'  # example 1:  Impute NA as one stratum in experimental arm; default model 
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
#' ps_res1 <- ps_match_strata(
#'   data.in = clinical_1, formula = indicator ~ BECOG + SEX + BNLR,
#'   indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", trt = "ARM",
#'   class.of.int = list("strata_1" = 1, "strata_2" = 0)
#'  )
#'  ## Weighted HRs
#'  ps_res1$stat
#'  
#'  # example 2: "doubly weighted control" model
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
#' 
#' ps_res2 <- ps_match_strata(
#'   data.in = clinical_2, formula = indicator ~ BECOG + SEX + BNLR, model = "dwc",
#'   indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", trt = "ARM",
#'   class.of.int = list("strata_1" = 0, "strata_2" = 1, "missing" = 2)
#'  ) 
#'  ps_res2$stat
#'  ps_res2$converged 
#'   
#' 
#' @import survival 
#' @importFrom dplyr mutate filter case_when bind_rows `%>%`
#' @importFrom  stats glm predict pnorm as.formula coef binomial quasibinomial
#' @importFrom rlang `:=` quo_name `!!`
#' @importFrom plyr `.`
#' @importFrom Matching Match
#' @importFrom pROC roc
#' @importFrom nnet multinom 
#' 
#' @export
#' 
ps_match_strata <- function(data.in, formula, indicator.var = "indicator", ties = TRUE,
                              class.of.int = NULL, tte = "AVAL", event = "event", trt = "trt",
                              response = NULL, caliper = NULL, model = "plain", weights = NULL, multinom.maxit = 100,
                              return.data = TRUE) {
  bi <- !is.null(response)

  if (!is.null(weights)) {
    stopifnot(length(weights) == nrow(data.in))
    data.in$weights <- weights
  } else {
    data.in$weights <- 1
  }

  stopifnot(class(data.in[, trt]) == "factor")
  stopifnot(nlevels(data.in[, trt]) == 2)

  data.withlabel <- data.in[which(data.in[, indicator.var] != -1), ]
  data.nolabel <- data.in[which(data.in[, indicator.var] == -1), ]
  data.pred <- data.in
  data.pred$labelnum <- ifelse(data.pred[, indicator.var] == -1, 0, 1)
  # Match() needs 0 and 1; 1: reference; 0: matching pool

  n.class <- length(unique(data.withlabel[, indicator.var, drop = TRUE]))
  stopifnot(n.class >= 2)
  stopifnot(max(data.withlabel[, indicator.var]) == n.class - 1) # need to be consequent 0,1,2...
  stopifnot(min(data.withlabel[, indicator.var]) == 0) # need to start from 0
  if (bi) stopifnot(all(data.in[, response] %in% c(0, 1)))
  if (!is.null(class.of.int)) stopifnot(class(class.of.int) == "list")
  if (is.null(class.of.int)) class.of.int <- as.list(sort(unique(data.withlabel[, indicator.var])))
  if (is.null(names(class.of.int))) names(class.of.int) <- sapply(class.of.int, function(i) paste0(i, collapse = ",")) # those will be row names to show in output table

  fitted.model <- roc.list <- vector("list", 0)
  # family <- ifelse(is.null(weights), binomial, quasibinomial)
  family <- if (is.null(weights))  binomial() else quasibinomial()


  if (model == "dwc") stopifnot(length(unique(data.in[, indicator.var, drop = TRUE])) != 3)

  if (n.class == 2) {
    # fit.glm <- glm(formula,family=binomial,data=data.withlabel) # within trt arm fitting
    # if two classes, logistic regression
    fit.glm <- tryCatch(glm(formula, family = family, data = data.withlabel, weights = weights),
      warning = function(w) list(glm(formula, family = family, data = data.withlabel, weights = weights), w)
    )
    if ("glm" %in% class(fit.glm)) {
      glm_warn <- FALSE
      warning.msg <- list(fit.glm = NULL)
    } else {
      glm_warn <- TRUE
      warning.msg <- list(fit.glm = fit.glm[[2]]$message)
      fit.glm <- fit.glm[[1]]
    }
    conv <- ifelse(fit.glm$converged == TRUE, TRUE, FALSE)

    data.pred$pred1 <- predict(fit.glm, newdata = data.pred, type = "response")
    data.pred$pred0 <- 1 - data.pred$pred1
    fitted.model$fit.glm <- fit.glm
    roc.list$fit.glm <- roc(data.withlabel[, indicator.var] ~ predict(fit.glm))
  }

  if (n.class > 2) {
    if (model == "dwc") {
      glm_warn <- FALSE
      warning.msg <- list(fit.layer1 = NULL, fit.layer2 = NULL)

      # first fit level1 vs others
      data.withlabel.layer1 <- data.withlabel %>%
        mutate(!!quo_name(indicator.var) := ifelse(.[, indicator.var] == 2, 1, 0))

      fit.layer1 <- tryCatch(fit.layer1 <- glm(formula, family = family, data = data.withlabel.layer1, weights = weights),
        warning = function(w) list(fit.layer1 <- glm(formula, family = family, data = data.withlabel.layer1, weights = weights), w)
      )
      if (!"glm" %in% class(fit.layer1)) {
        glm_warn <- TRUE
        warning.msg$fit.layer1 <- fit.layer1[[2]]$message
        fit.layer1 <- fit.layer1[[1]]
      }
      # subset to non level1
      data.withlabel.layer2 <- data.withlabel %>% dplyr::filter(.[, indicator.var] != 2)

      fit.layer2 <- tryCatch(fit.layer2 <- glm(formula, family = family, data = data.withlabel.layer2, weights = weights),
        warning = function(w) list(fit.layer2 <- glm(formula, family = family, data = data.withlabel.layer2, weights = weights), w)
      )
      if (!"glm" %in% class(fit.layer2)) {
        glm_warn <- TRUE
        warning.msg$fit.layer2 <- fit.layer2[[2]]$message
        fit.layer2 <- fit.layer2[[1]]
      }


      data.pred$pred.layer1 <- predict(fit.layer1, newdata = data.pred, type = "response")
      data.pred$pred.layer2 <- predict(fit.layer2, newdata = data.pred, type = "response")
      data.pred <- data.pred %>%
        mutate(
          pred0 = (1 - .data$pred.layer1) * (1 - .data$pred.layer2),
          pred1 = (1 - .data$pred.layer1) * .data$pred.layer2,
          pred2 = .data$pred.layer1
        ) %>%
        dplyr::select(-.data$pred.layer1, -.data$pred.layer2)
      conv <- ifelse(fit.layer1$converged == TRUE && fit.layer2$converged == TRUE, TRUE, FALSE)
      fitted.model$fit.layer1 <- fit.layer1
      fitted.model$fit.layer2 <- fit.layer2
      roc.list$fit.layer1 <- roc(data.withlabel.layer1[, indicator.var] ~ predict(fit.layer1))
      roc.list$fit.layer2 <- roc(data.withlabel.layer2[, indicator.var] ~ predict(fit.layer2))
    }

    if (model == "plain") {
      # fit.glm <- multinom(formula, data=data.withlabel, trace=FALSE) # within trt arm fitting;
      # if more than two classes, multinomial regression
      fit.glm <- tryCatch(fit.glm <- multinom(formula, data = data.withlabel, trace = FALSE, weights = weights, maxit = multinom.maxit),
        warning = function(w) list(fit.glm <- multinom(formula, data = data.withlabel, trace = FALSE, weights = weights, maxit = multinom.maxit), w)
      )
      if ("multinom" %in% class(fit.glm)) {
        glm_warn <- FALSE
        warning.msg <- list(fit.glm = NULL)
      } else {
        glm_warn <- TRUE
        warning.msg <- list(fit.glm = fit.glm[[2]]$message)
        fit.glm <- fit.glm[[1]]
      }
      conv <- ifelse(fit.glm$convergence == 0, TRUE, FALSE)

      multinom.summary <- vector("list", 0)
      fit.glm <- do.call("multinom", list(formula = formula, data = data.withlabel, trace = FALSE, weights = data.withlabel$weights, maxit = multinom.maxit))
      # for some reason need do.call to make summary() works
      multinom.summary$coef <- stats::coef(fit.glm)
      multinom.summary$standard.error <- summary(fit.glm)$standard.errors
      multinom.summary$pval <- (1 - pnorm(abs(multinom.summary$coef / multinom.summary$standard.error), 0, 1)) * 2

      pred.mat <- predict(fit.glm, newdata = data.pred, "probs")
      for (i in colnames(pred.mat)) {
        data.pred[, paste0("pred", i)] <- pred.mat[, i]
      }
      fitted.model$fit.multinom <- fit.glm
    }
  }



  stat.mat <- matrix(NA, nrow = length(class.of.int), ncol = ifelse(bi, 3, 5)) # if binary endpoint, only one summary statistics will be provided
  if (!bi) colnames(stat.mat) <- c("Coefficient", "Variance", "HR", "CI.Lower", "CI.Upper")
  if (bi) colnames(stat.mat) <- c("Diff", "ctrlORR", "trtORR")
  rownames(stat.mat) <- names(class.of.int)
  data.out <- NULL
  warning.msg$cox <- NULL
  for (i in 1:length(class.of.int)) {
    data.tmp <- data.pred
    class.nums <- class.of.int[[i]]
    class.name.ofint <- ifelse(length(class.nums) == 1, paste0("pred", class.nums), paste0("pred", paste0(class.nums, collapse = "or")))
    if (length(class.nums) == 1) data.tmp$pred <- data.tmp[, class.name.ofint]
    if (length(class.nums) > 1) {
      data.tmp$pred <- rowSums(data.tmp[, paste0("pred", class.nums)])
      # if it is a combined class, e.g. prob of being either class2 or 3
      # is of interest, calculate the prob(being either class) by taking
      # prob(class2)+prob(class3)
      data.pred[[class.name.ofint]] <- rowSums(data.pred[, paste0("pred", class.nums)])
    }
    if (!is.null(weights)) {
      stopifnot(length(weights) == nrow(data.tmp))
      data.tmp$weights <- weights
    }
    data.use <- data.tmp[which(data.tmp[, indicator.var] %in% c(-1, class.of.int[[i]])), ] # take all
    # the no label pts, and class i pts in trt arm
    y.var <- ifelse(bi, response, tte)
    if (ties) {
      rr1 <- Match(Y = data.use[, y.var], Tr = data.use[, "labelnum"], X = data.use$pred, M = 1, ties = TRUE, caliper = caliper)
      # M=number of matches to be found,
      # Match known using pts with unknown label (regardless which arm they are in)
      mid1 <- c(rr1$index.treated, rr1$index.control)
      m1 <- data.use[mid1, ] # Enhance data with additionala ttributes
      wt1 <- rr1$weights # Vector of weights. There is one weight for each matched-pair in the matched dataset.
      m1[, paste0(class.name.ofint, "_glm")] <- m1[, "pred"] # save the predicted probability from pred_glm
      m1[, class.name.ofint] <- rep(wt1, 2) # Dubplicate to assign weight for treated and control
      if (!is.null(weights)) {
        idx.trt <- rr1$index.treated
        names(idx.trt) <- idx.trt
        trt.ori.wts <- data.use[rr1$index.treated, "weights"]
        names(trt.ori.wts) <- idx.trt
        trt.wts <- wt1 * trt.ori.wts # reweight Match weights by input weights
        ctrl.wts.in <- data.use[rr1$index.control, "weights"] # ctrl samples input weight
        ctrl.wts.sum <- tapply(ctrl.wts.in, rr1$index.treated, sum) # denominator to use, reweight multipl control samples mateched to the same trt
        ctrl.wts.sum.each <- ctrl.wts.sum[names(idx.trt)]
        ctrl.wts <- ctrl.wts.in / ctrl.wts.sum.each * trt.ori.wts # input weight / sum input weight for same matched trt * match weight
        m1[, class.name.ofint] <- c(trt.wts, ctrl.wts)
      }
      if (!bi) {
        cox.fit <- tryCatch(cox.fit <- coxph(as.formula(paste("Surv(", tte, ",", event, ")~", trt)),
          weights = m1[, class.name.ofint], ties = "breslow", data = m1, robust = TRUE
        ),
        warning = function(w) {
          list(cox.fit <- coxph(as.formula(paste("Surv(", tte, ",", event, ")~", trt)),
            weights = m1[, class.name.ofint], ties = "breslow", data = m1, robust = TRUE
          ),
          warning = w
          )
        }
        )
        if ("warning" %in% names(cox.fit)) {
          warning.msg$cox <- c(warning.msg$cox, cox.fit[[2]]$message)
          cox.fit <- cox.fit[[1]]
        } else {
          warning.msg$cox <- warning.msg$cox
        }
      }
      if (bi) {
        m1$wresp <- m1[, class.name.ofint] * m1[, response]
        d1 <- m1[which(m1[, indicator.var] == -1), ]
        d2 <- m1[which(m1[, indicator.var] != -1), ]
        stat.mat[i, ] <- c(sum(d2$wresp) / sum(d2[, class.name.ofint]) - sum(d1$wresp) / sum(d1[, class.name.ofint]), sum(d1$wresp) / sum(d1[, class.name.ofint]), sum(d2$wresp) / sum(d2[, class.name.ofint]))
      }
      data.out <- bind_rows(data.out, m1[, c(colnames(data.in), class.name.ofint, paste0(class.name.ofint, "_glm"))]) # only keep matched samples in the output data
    }
    if (!ties) {
      rr1 <- Match(Y = data.use[, y.var], Tr = data.use[, "labelnum"], X = data.use$pred, M = 1, ties = FALSE, caliper = caliper)
      mid1 <- c(rr1$index.treated, rr1$index.control)
      # data.use$weights <- 1
      m1 <- data.use[mid1, ] # Enhance data with additionala ttributes

      if (!bi) {
        cox.fit <- tryCatch(cox.fit <- coxph(as.formula(paste("Surv(", tte, ",", event, ")~", trt)),
          weights = m1[, "weights"], ties = "breslow", data = m1
        ),
        warning = function(w) {
          list(cox.fit <- coxph(as.formula(paste("Surv(", tte, ",", event, ")~", trt)),
            weights = m1[, "weights"], ties = "breslow", data = m1
          ),
          warning = w
          )
        }
        )
        if ("warning" %in% names(cox.fit)) {
          warning.msg$cox <- c(warning.msg$cox, cox.fit[[2]]$message)
          cox.fit <- cox.fit[[1]]
        } else {
          warning.msg$cox <- warning.msg$cox
        }
      }
      if (bi) {
        m1$wresp <- m1[, response]
        d1 <- data.use[which(data.use[, trt] == levels(data.use[, trt])[1]), ]
        d2 <- data.use[which(data.use[, trt] != levels(data.use[, trt])[1]), ]
        stat.mat[i, ] <- c(mean(d2$wresp) - mean(d1$wresp), mean(d1$wresp), mean(d2$wresp))
      }
    }
    if (!bi) {
      cox.sum <- summary(cox.fit)
      stat.mat[i, ] <- c(as.numeric(cox.sum$coefficients[1]), as.numeric(cox.sum$coefficients[3])^2, as.numeric(cox.sum$conf.int[c(1, 3, 4)]))
    }
  }
  output <- list("stat" = stat.mat, "converged" = conv, "any_warning_glm" = glm_warn, warning.msg = warning.msg, models = fitted.model, roc.list = roc.list)
  if (return.data) output <- c(output, list("data" = data.out))
  if (exists("multinom.summary")) output <- c(output, list("multinom.summary" = multinom.summary))
  return(output)
}
