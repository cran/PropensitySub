#' Inverse Probability weighting of strata (two or more strata, survival or binary endpoint)
#' 
#' This function performs inverse probability weighting of two or more strata. \cr 
#' It could be used when arm1 has 2 or more strata, while stratum information is unknown in arm2. \cr
#' The function will fit a logistic regression (when 2 classes) or multinomial logistic
#' regression (when > 2 classes) based on stratum labels in arm1 (model: label ~ features), then
#' predict stratum labels for pts in arm2 based on the fitted model (as well as pts in arm1 who have missing labels,
#' if there is any).
#' The predicted probability of being stratum X will be used as weights when estimating treatment difference of two arms
#' (Hazard ratio for survival endpoint; response rate difference for binary endpoint)
#'
#' @md
#' @param data.in (`data.frame`) input data
#' @param formula (`formula`) to input to the logistic or multinomial logistic model (in the form of strata~features)
#' @param indicator.var (`string`) column name of the strata indicator variable which must be numeric.
#' Assume arm1 has strata labeling and arm2 does not have strata labeling.
#' pts without strata labeling should be indicated as -1 (e.g. pts in the arm1, or pts in arm2 but with missing label).
#' within arm1 (the arm with strata labeling), subclasss should be indicated as 0,1,2...
#' @param tte (`string`) column name of the time to event variable
#' @param event (`string`) column name of the event variable (1: event, 0: censor)
#' @param trt (`string`) column name of the treatment variable. The variable is expected to be in `factor` format and the first level
#' will be considered as reference (control) arm when calculating summary statistics.
#' @param response (`string`) column name of the response variable. 1 indicates responder and 0 indicates non responder.
#' if response is not NULL, tte and event will be ignored and the function will assume binary outcome.
#' @param class.of.int (`list`) classes (stratum) of interest. Request to be in list format.
#' It could be subset of classes in arm1; it could also define combined classes.
#' For example: class.of.int = list("class1"=0, "class2"=1, "class3"=2, "class2or3"="c(1,2)").
#' for "class2or3", Prob(class 2 or 3) will be calculated as Prob(class2) + Prob(class3)
#' @param return.data (`logical`) whether to return data with estimated probabilities.
#' @param model (`string`) one of (`plain`, `dwc`, `wri`).
#' \describe{
#'   \item{"plain"}{when 2 levels are specified in indicator variable, a binomial glm will be fitted;
#' when more than 2 levels are specified, a multinomial glm will be fitted;}
#'   \item{"dwc"}{Doubly weighted control: Two separated models will be fitted: one is binomial glm of 2 vs. (1, 0), the
#' other one is bionomial glm of 1 vs 0. The probability of being each class is then calculated by aggregating
#' these two models. Note this is similar to the plain method but with different (less rigid) covariance assumptions.}
#'   \item{"wri"}{Weight regression imputation: the current status is going to be learned from the next status.
#' Indicator of the next status should be specified using indicator.next.
#' Currently "wri" only support the case where there are only two non-missing strata.
#' In indicator variable, the two nonmissing strata should be coded as 0 and 1, the missing group should
#' be coded as 2.}
#' }
#' @param indicator.next (`string`) column name of the column which indicates status at a different measurement.
#' It should be coded in the same way as in indicator.var (e.g. -1, 0, 1). Patients who have both
#' missing current status and missing next status should be excluded in the modeling.
#' @param weights (`numeric`) weights of each subject. If not `NULL`, the estimated probabilities will be reweightsed to
#' ensure sum(probability) of a subject = the subject's weights. If weights is not `NULL`, quasibinomial model will be used.
#' @param multinom.maxit see parameter `maxit` in [nnet::multinom], default is 100
#' 
#' 
#' @import survival 
#' @importFrom dplyr mutate filter case_when bind_rows `%>%`
#' @importFrom  stats glm predict pnorm as.formula coef binomial quasibinomial
#' @importFrom rlang `:=` quo_name `!!`
#' @importFrom plyr `.`
#' @importFrom pROC roc
#' @importFrom nnet multinom 
#' 
#' @export
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
#'  # example 1: Impute NA as one stratum in experimental arm; default model 
#'  library(dplyr)
#'  clinical_1 <- clinical %>% mutate( 
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
#' ipw_res1 <- ipw_strata(
#'   data.in = clinical_1, formula = indicator ~ BECOG + SEX + BNLR,
#'   indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", trt = "ARM",
#'   class.of.int = list("strata_1" = 1, "strata_2" = 0)
#'  )
#'  ## Weighted HRs
#'  ipw_res1$stat
#'  
#'  # example 2: "Weight regression imputation" model
#' clinical_2 <- clinical %>% mutate( 
#'   indicator = case_when(
#'     STRATUM == "strata_1" ~ 0, 
#'     STRATUM == "strata_2" ~ 1, 
#'     is.na(STRATUM) & ARM == "experimental" ~ 2,
#'     TRUE ~ -1
#'   ),
#'   indicator_next = case_when(
#'     STRATUM_NEXT == "strata_1" ~ 0, 
#'     STRATUM_NEXT == "strata_2" ~ 1, 
#'     is.na(STRATUM_NEXT) & ARM == "experimental" ~ 2,
#'     TRUE ~ -1
#'   ),
#'   ARM  = factor(ARM, levels = c("control","experimental")),
#'   BNLR = case_when(
#'     is.na(BNLR) ~ median(BNLR, na.rm = TRUE),
#'     TRUE ~ BNLR
#'   )
#' )
#' 
#' ipw_res2 <- ipw_strata(
#'   data.in = clinical_2, formula = indicator ~ BECOG + SEX + BNLR, model = "wri",
#'   indicator.var = "indicator", indicator.next = "indicator_next",
#'    tte = "OS_MONTH", event = "OS_EVENT", trt = "ARM",
#'   class.of.int = list("strata_1" = 1, "strata_2" = 0)
#'  )
#'  ## Weighted HRs
#'  ipw_res2$stat 
#'  
#' @note Three elements in the output list - the data element is a data frame that contains input data and
#' estimated probabilities. The stat element contains estimated treatment difference between 2 arms, \cr 
#' in each of the strata of interest.
#' The converge element indicates whether the model converged (taking from $converged from [stats::glm] \cr 
#' and $convergency from [nnet::multinom]).
#' if return.data is `FALSE`, data won't be returned.
#' 

ipw_strata <- function(data.in, formula, indicator.var = "indicator",
                         class.of.int = NULL, tte = "AVAL", event = "event", trt = "trt",
                         response = NULL, model = "plain", indicator.next = NULL,
                         weights = NULL, multinom.maxit = 100,
                         return.data = TRUE) {
  bi <- !is.null(response)
  if (!is.null(weights)) {
    stopifnot(length(weights) == nrow(data.in))
    data.in$weights <- weights
  } else {
    data.in$weights <- 1
  }

  data.withlabel <- data.in[which(data.in[, indicator.var] != -1), ]
  data.nolabel <- data.in[which(data.in[, indicator.var] == -1), ]

  n.class <- length(unique(data.withlabel[, indicator.var, drop = TRUE]))
  stopifnot(n.class >= 2)
  stopifnot(class(data.in[, trt, drop = TRUE]) == "factor")
  stopifnot(nlevels(data.in[, trt, drop = TRUE]) == 2)
  stopifnot(max(data.withlabel[, indicator.var]) == n.class - 1) # need to be consequent 0,1,2...
  stopifnot(min(data.withlabel[, indicator.var]) == 0) # need to start from 0
  if (bi) stopifnot(all(data.in[, response] %in% c(0, 1)))
  if (!is.null(class.of.int)) stopifnot(class(class.of.int) == "list")
  if (is.null(class.of.int)) class.of.int <- as.list(sort(unique(data.withlabel[, indicator.var])))
  if (is.null(names(class.of.int))) names(class.of.int) <- sapply(class.of.int, function(i) paste0(i, collapse = ",")) # those will be row names to show in output table
  fitted.model <- roc.list <- vector("list", 0)
  # family <- ifelse(is.null(weights), binomial, quasibinomial)
  family <- if (is.null(weights))  binomial() else quasibinomial()
  
  f1 <- as.character(formula)
  vars.in.formula <- strsplit(f1[[3]], split = " \\+ ")[[1]]
  stopifnot(all(vars.in.formula %in% colnames(data.in)))
  if (any(is.na(data.in[, vars.in.formula]))) stop("Matrix of dependent variables should not contain NAs")

  if (model == "dwc") stopifnot(length(unique(data.in[, indicator.var])) != 3)

  if (n.class == 2) {
    # fit.glm <- glm(formula,family=binomial,data=data.withlabel) # within trt arm fitting
    fit.glm <- tryCatch(fit.glm <- glm(formula, family = family, data = data.withlabel, weights = weights),
      warning = function(w) list(fit.glm <- glm(formula, family = family, data = data.withlabel, weights = weights), w)
    )
    if ("glm" %in% class(fit.glm)) {
      glm_warn <- FALSE
      warning.msg <- list(fit.glm = NULL)
    } else {
      glm_warn <- TRUE
      warning.msg <- list(fit.glm = fit.glm[[2]]$message)
      fit.glm <- fit.glm[[1]]
    }
    # if two classes, logistic regression
    conv <- ifelse(fit.glm$converged == TRUE, TRUE, FALSE)

    data.nolabel$pred1 <- predict(fit.glm, newdata = data.nolabel, type = "response") # Predict subclasss
    # for control patients (prob to be class nonzero);
    # not binary since we want the weights

    data.nolabel$pred0 <- 1 - data.nolabel$pred1
    data.withlabel$pred1 <- data.withlabel[, indicator.var] # binary in trt arm
    data.withlabel$pred0 <- 1 - data.withlabel[, indicator.var]
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
      data.withlabel.layer2 <- data.withlabel %>% filter(.[, indicator.var] != 2)
      # fit strata within level1
      fit.layer2 <- tryCatch(fit.layer2 <- glm(formula, family = family, data = data.withlabel.layer2, weights = weights),
        warning = function(w) list(fit.layer2 <- glm(formula, family = family, data = data.withlabel.layer2, weights = weights), w)
      )
      if (!"glm" %in% class(fit.layer2)) {
        glm_warn <- TRUE
        warning.msg$fit.layer2 <- fit.layer2[[2]]$message
        fit.layer2 <- fit.layer2[[1]]
      }

      
      data.nolabel$pred.layer1 <- predict(fit.layer1, newdata = data.nolabel, type = "response")
      data.nolabel$pred.layer2 <- predict(fit.layer2, newdata = data.nolabel, type = "response")
      data.nolabel <- data.nolabel %>%
        mutate(
          pred0 = (1 - .data$pred.layer1) * (1 - .data$pred.layer2),
          pred1 = (1 - .data$pred.layer1) * .data$pred.layer2,
          pred2 = .data$pred.layer1
        ) %>%
        dplyr::select(-.data$pred.layer1, -.data$pred.layer2)
      data.withlabel <- data.withlabel %>% mutate(
        pred0 = case_when(
          .[, indicator.var] == 0 ~ 1,
          .[, indicator.var] == 1 ~ 0,
          TRUE ~ 0
        ),
        pred1 = case_when(
          .[, indicator.var] == 0 ~ 0,
          .[, indicator.var] == 1 ~ 1,
          TRUE ~ 0
        ),
        pred2 = ifelse(.[, indicator.var] == 2, 1, 0)
      )

      conv <- ifelse(fit.layer1$converged == TRUE && fit.layer2$converged == TRUE, TRUE, FALSE)
      fitted.model$fit.layer1 <- fit.layer1
      fitted.model$fit.layer2 <- fit.layer2
      roc.list$fit.layer1 <- roc(data.withlabel.layer1[, indicator.var] ~ predict(fit.layer1))
      roc.list$fit.layer2 <- roc(data.withlabel.layer2[, indicator.var] ~ predict(fit.layer2))
    }
    if (model == "wri") {
      glm_warn <- FALSE
      warning.msg <- list(fit.current.by.baseNnext = NULL, fit.next = NULL)
      f1 <- as.character(formula)
      formula.current.by.baseNnext <- as.formula(paste(f1[2], "~", f1[3], "+", indicator.next))

      formula.next <- as.formula(paste(indicator.next, "~", f1[3]))

      # to build model of current status ~ baseline + next status,
      # using data that are non missing in both current and next
      data.current.by.baseNnext.valid <- data.withlabel[which(data.withlabel[, indicator.var] != 2 & data.withlabel[, indicator.next] != 2), ]
      # How does FPSTATnum coded??? 0 is neg and 1 is pos???
      fit.current.by.baseNnext <- tryCatch(fit.current.by.baseNnext <- glm(formula.current.by.baseNnext, family = family, data = data.current.by.baseNnext.valid, weights = weights),
        warning = function(w) {
          list(fit.current.by.baseNnext <- glm(formula.current.by.baseNnext,
            family = family,
            data = data.current.by.baseNnext.valid, weights = weights
          ), w)
        }
      )
      if (!"glm" %in% class(fit.current.by.baseNnext)) {
        glm_warn <- TRUE
        warning.msg$fit.current.by.baseNnext <- fit.current.by.baseNnext[[2]]$message
        fit.current.by.baseNnext <- fit.current.by.baseNnext[[1]]
      }
      data.current.missing <- data.withlabel[which(data.withlabel[, indicator.var] == 2), ]
      data.current.missing$pred1 <- predict(fit.current.by.baseNnext, newdata = data.current.missing, type = "response")
      data.current.missing$pred0 <- 1 - data.current.missing$pred1
      # for current missing population, pred(current | next) will be used

      # model of next status ~ baseline
      # for control population, pred(current| next = +) pred(next=+) + pred(current|next=-) Pred(next=-) will be used
      data.next.evaluable <- data.withlabel[which(data.withlabel[, indicator.next] != 2), ]
      fit.next <- tryCatch(fit.next <- glm(formula.next, family = family, data = data.next.evaluable, weights = weights),
        warning = function(w) list(fit.next <- glm(formula.next, family = family, data = data.next.evaluable, weights = weights), w)
      )
      if (!"glm" %in% class(fit.next)) {
        glm_warn <- TRUE
        warning.msg$fit.next <- fit.next[[2]]$message
        fit.next <- fit.next[[1]]
      }
      data.nolabel$pred.current.by.baseNnext.0 <- predict(fit.current.by.baseNnext, newdata = data.nolabel %>% mutate(indicator.next = 0), type = "response")
      data.nolabel$pred.current.by.baseNnext.1 <- predict(fit.current.by.baseNnext, newdata = data.nolabel %>% mutate(indicator.next = 1), type = "response")
      data.nolabel$pred.current <- predict(fit.next, newdata = data.nolabel, type = "response")
      data.nolabel <- data.nolabel %>%
        mutate(pred1 = .data$pred.current.by.baseNnext.0 * (1 - .data$pred.current) + 
                 .data$pred.current.by.baseNnext.1 * .data$pred.current) %>%
        mutate(pred0 = 1 - .data$pred1) %>%
        mutate(pred2 = 0) # need pred2 to make sure downstream functions won break

      data.current.nonmissing <- data.withlabel %>% filter(.[, indicator.var] %in% c(0, 1))
      data.current.nonmissing <- data.current.nonmissing %>% mutate(
        pred0 = case_when(
          .[, indicator.var] == 0 ~ 1,
          .[, indicator.var] == 1 ~ 0
        ),
        pred1 = case_when(
          .[, indicator.var] == 0 ~ 0,
          .[, indicator.var] == 1 ~ 1
        )
      )

      data.withlabel <- bind_rows(data.current.missing, data.current.nonmissing)
      data.withlabel$pred2 <- 0

      conv <- ifelse(fit.current.by.baseNnext$converged == TRUE && fit.next$converged == TRUE, TRUE, FALSE)


      fitted.model$fit.current.by.baseNnext <- fit.current.by.baseNnext
      fitted.model$fit.next <- fit.next
      fitted.model$fit.current.by.nextind <- glm(as.formula(paste(indicator.var, "~", indicator.next)), family = family, data = data.current.by.baseNnext.valid, weights = weights) # for checking only

      roc.list$fit.current.by.baseNnext <- roc(data.current.by.baseNnext.valid[, indicator.var] ~ predict(fit.current.by.baseNnext))
      roc.list$fit.next <- roc(data.next.evaluable[, indicator.next] ~ predict(fit.next))
      roc.list$fit.current.by.nextind <- roc(data.current.by.baseNnext.valid[, indicator.var] ~ predict(fitted.model$fit.current.by.nextind))

      # for(cc in 1:length(class.of.int)) class.of.int[[cc]] <- setdiff(class.of.int[[cc]],2)
      # check.class <- sapply(class.of.int, length)
      # class.of.int <- class.of.int[check.class>0] # cannot produce any class where level 2 is included
    }

    if (model == "plain") {
      # fit.glm <- multinom(formula, data=data.withlabel, trace=FALSE) # within trt arm fitting;
      # if more than two classes, multinomial regression
      fit.glm <- tryCatch(
        multinom(formula, data = data.withlabel, trace = FALSE, weights = weights, maxit = multinom.maxit),
        warning = function(w) {
          list(multinom(
            formula,
            data = data.withlabel, trace = FALSE, weights = weights, maxit = multinom.maxit
          ), w)
        }
      )
      if ("multinom" %in% class(fit.glm)) {
        glm_warn <- FALSE
        warning.msg <- list(fit.glm = NULL)
      } else {
        glm_warn <- TRUE
        warning.msg$fit.glm <- fit.glm[[2]]$message
        fit.glm <- fit.glm[[1]]
      }
      conv <- ifelse(fit.glm$convergence == 0, TRUE, FALSE)
      multinom.summary <- vector("list", 0)
      fit.glm <- do.call("multinom", list(formula = formula, data = data.withlabel, trace = FALSE, weights = data.withlabel$weights, maxit = multinom.maxit))
      # for some reason need do.call to make summary() works
      multinom.summary$coef <- stats::coef(fit.glm)
      multinom.summary$standard.error <- summary(fit.glm)$standard.errors
      multinom.summary$pval <- (1 - pnorm(abs(multinom.summary$coef / multinom.summary$standard.error), 0, 1)) * 2
      pred.mat <- predict(fit.glm, newdata = data.nolabel, "probs")
      for (i in colnames(pred.mat)) {
        data.nolabel[, paste0("pred", i)] <- pred.mat[, i]
        data.withlabel[, paste0("pred", i)] <-
          ifelse(data.withlabel[, indicator.var] == as.numeric(i), 1, 0)
      }
      fitted.model$fit.multinom <- fit.glm
    }
  }

  data.pred <- bind_rows(data.nolabel, data.withlabel) # full population; with predicted value
  if (!is.null(weights)) {
    new.column.names <- setdiff(colnames(data.pred), colnames(data.in))
    pred.column.names <- grep("pred", new.column.names, value = TRUE)
    data.pred[, pred.column.names] <- data.pred[, pred.column.names] * data.pred$weights
  }


  stat.mat <- matrix(NA, nrow = length(class.of.int), ncol = ifelse(bi, 3, 5)) # if binary endpoint, only one summary statistics will be provided
  if (!bi) colnames(stat.mat) <- c("Coefficient", "Variance", "HR", "CI.Lower", "CI.Upper")
  if (bi) colnames(stat.mat) <- c("Diff", "ctrlORR", "trtORR")
  rownames(stat.mat) <- names(class.of.int)
  warning.msg$cox <- NULL
  for (i in 1:length(class.of.int)) {
    data.tmp <- data.pred
    class.nums <- class.of.int[[i]]
    if (length(class.nums) == 1) data.tmp$pred <- data.tmp[, paste0("pred", class.nums)]
    if (length(class.nums) > 1) {
      data.tmp$pred <- rowSums(data.tmp[, paste0("pred", class.nums)])
      # if it is a combined class, e.g. prob of being either class2 or 3
      # is of interest, calculate the prob(being either class) by taking
      # prob(class2)+prob(class3)
      data.pred[[paste0("pred", paste0(class.nums, collapse = "or"))]] <- rowSums(data.pred[, paste0("pred", class.nums)])
    }
    data.use <- data.tmp[which(data.tmp[, indicator.var] %in% c(-1, class.of.int[[i]])), ] # take all
    if (model == "wri") data.use <- data.tmp[which(data.tmp[, indicator.var] %in% c(-1, 2, class.of.int[[i]])), ] # if wri, missing in trt are also predicted
    # the nolabel pts, and class i pts
    tmp <- .Machine
    data.use[which(data.use$pred == 0), "pred"] <- tmp$double.xmin # impute numerical underflow to smallest number

    if (!bi) {
      cox.fit <- tryCatch(cox.fit <- coxph(as.formula(paste("Surv(", tte, ",", event, ")~", trt)),
        weights = data.use$pred, ties = "breslow", data = data.use
      ),
      warning = function(w) {
        list(cox.fit <- coxph(as.formula(paste("Surv(", tte, ",", event, ")~", trt)),
          weights = data.use$pred, ties = "breslow", data = data.use
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
      cox.sum <- summary(cox.fit)
      stat.mat[i, ] <- c(as.numeric(cox.sum$coefficients[1]), as.numeric(cox.sum$coefficients[3])^2, as.numeric(cox.sum$conf.int[c(1, 3, 4)]))
    }

    if (bi) {
      data.use$wresp <- data.use$pred * data.use[, response]
      d1 <- data.use[which(data.use[, trt] == levels(data.use[, trt])[1]), ]
      d2 <- data.use[which(data.use[, trt] != levels(data.use[, trt])[1]), ]
      stat.mat[i, ] <- c(sum(d2$wresp) / sum(d2$pred) - sum(d1$wresp) / sum(d1$pred), sum(d1$wresp) / sum(d1$pred), sum(d2$wresp) / sum(d2$pred))
    }
  }
  output <- list("stat" = stat.mat, "converged" = conv, "any_warning_glm" = glm_warn, warning.msg = warning.msg, models = fitted.model, roc.list = roc.list)
  if (return.data) output <- c(output, list("data" = data.pred))
  if (exists("multinom.summary")) output <- c(output, list("multinom.summary" = multinom.summary))

  return(output)
}
