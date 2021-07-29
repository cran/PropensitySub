#' Compare weighted and unweighted (naive analysis) standardized difference
#' 
#' @md
#' @inheritParams calc_std_diff
#' @inheritParams km_plot_weight 
#' @param data.in.unadj data set to use for the unadjusted analysis. For example, if PSM is used, the
#' adjusted analysis should be done on the matched population but the unadjusted analysis should be done on
#' the original population
#' @param return.levels whether to return levels of each factor within each class.
#' @param subj.aggr whether aggregate multiple entries from the same patients to one record
#' @param usubjid.var column name indiacts subjuect id
#' @return return a `list`, each `list` element is a `data.frame` 
#'  containing absolute standardized difference for each variable.
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
#' ipw_res1 <- ipw_strata(
#'   data.in = clinical_1, formula = indicator ~ BECOG + SEX + BNLR,
#'   indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", trt = "ARM",
#'   class.of.int = list("strata_1" = 1, "strata_2" = 0)
#'  )
#'  std_diff(
#'   data.in = ipw_res1$data, vars = c("BECOG", "SEX", "BNLR"),
#'   indicator.var = "indicator", trt = "ARM",
#'   class.of.int = list("strata_1" = 1, "strata_2" = 0),
#'   usubjid.var = "SUBJID"
#' )
#' @importFrom rlang `:=` quo_name `!!` sym
#' @importFrom dplyr group_by_at summarise distinct  
#' @export
#' @note Calculation from Austin and Stuart (2015)
std_diff <- function(data.in, data.in.unadj = NULL, trt, vars, indicator.var = "indicator",
                     class.of.int = NULL, prob.names = NULL, return.levels = FALSE, subj.aggr = TRUE,
                     usubjid.var = "USUBJID") {
  data.trt.name <- unique(data.in[which(data.in[, indicator.var] != -1), trt])
  stopifnot(length(data.trt.name) == 1)
  data.trt <- data.in[which(data.in[, trt] == data.trt.name), ] # take trt only
  data.ctrl <- data.in[which(data.in[, trt] != data.trt.name), ] # take control only
  data.trt.unadj <- data.in.unadj[which(data.in.unadj[, trt] == data.trt.name), ] # take trt only
  data.ctrl.unadj <- data.in.unadj[which(data.in.unadj[, trt] != data.trt.name), ] # take control only

  all.std.diff <- calc_std_diff(
    vars = vars, data0 = data.ctrl, weight0 = rep(1, nrow(data.ctrl)),
    data1 = data.trt, weight1 = rep(1, nrow(data.trt))
  ) # this is used to get all variable levels
  # in strata plots we will plot out all levels appeared in full population
  if (!is.null(data.in.unadj)) {
    all.std.diff <- calc_std_diff(
      vars = vars, data0 = data.ctrl.unadj,
      weight0 = rep(1, nrow(data.ctrl.unadj)), data1 = data.trt.unadj,
      weight1 = rep(1, nrow(data.trt.unadj))
    )
  }

  n.class <- length(unique(data.trt[, indicator.var]))
  stopifnot(n.class >= 2)
  stopifnot(max(data.trt[, indicator.var]) == n.class - 1) # need to be consequent 0,1,2...
  stopifnot(min(data.trt[, indicator.var]) == 0) # need to start from 0
  if (!is.null(class.of.int)) stopifnot(class(class.of.int) == "list")
  if (is.null(class.of.int)) class.of.int <- as.list(sort(unique(data.trt[, indicator.var])))
  if (is.null(names(class.of.int))) names(class.of.int) <- sapply(class.of.int, function(i) paste0(i, collapse = ",")) # those will be row names to show in output table
  if (is.null(prob.names)) prob.names <- sapply(class.of.int, function(i) paste0("pred", paste0(i, collapse = "or")))
  stopifnot(all(vars %in% colnames(data.in)))

  diff.list <- vector("list", length(class.of.int))
  factor.levels.list <- vector("list", length(class.of.int))
  names(diff.list) <- names(factor.levels.list) <- names(class.of.int)
  for (i in 1:length(class.of.int)) {
    grp_label <- names(class.of.int)[i]
    prob.var <- prob.names[i]
    data.trt.use <- data.trt[which(data.trt[, indicator.var] %in% c(class.of.int[[i]], -1) & !is.na(data.trt[, prob.var])), ] # selected class and unk
    data.ctrl.use <- data.ctrl[which(!is.na(data.ctrl[, prob.var])), ] # in PSM, only pts matched to this class have valid prob
    factor.levels.list[[i]]$trt <- sapply(data.trt.use[, vars], function(ii) unique(ii))
    factor.levels.list[[i]]$ctrl <- sapply(data.ctrl.use[, vars], function(ii) unique(ii))
    if (subj.aggr) {
      data.trt.aggr <- data.trt.use %>%
        group_by_at(usubjid.var) %>%
        summarise(!!prob.var := sum(!!sym(prob.var)))
      data.ctrl.aggr <- data.ctrl.use %>%
        group_by_at(usubjid.var) %>%
        summarise(!!prob.var := sum(!!sym(prob.var)))
      data.trt.use <- merge(data.trt.aggr, data.trt.use[, setdiff(colnames(data.trt.use), prob.var)], by = usubjid.var, all.x = TRUE, all.y = FALSE) %>% distinct()
      data.ctrl.use <- merge(data.ctrl.aggr, data.ctrl.use[, setdiff(colnames(data.ctrl.use), prob.var)], by = usubjid.var, all.x = TRUE, all.y = FALSE) %>% distinct()
    }

    if (is.null(data.in.unadj)) {
      res.adj.raw <- calc_std_diff(
        vars = vars, data0 = data.ctrl.use,
        weight0 = data.ctrl.use[, prob.var],
        data1 = data.trt.use,
        weight1 = data.trt.use[, prob.var]
      )
    }
    # If PSM, sometimes factor levels in matched data set is fewer than factor levels in original data set
    # still want to show all levels in original data set
    # so input original data set with weight 0 (just to get the level showing)
    if (!is.null(data.in.unadj)) {
      res.adj.raw <- calc_std_diff(
        vars = vars, data0 = bind_rows(data.ctrl.use, data.ctrl.unadj),
        weight0 = c(data.ctrl.use[, prob.var], rep(0, nrow(data.ctrl.unadj))),
        data1 = bind_rows(data.trt.use, data.trt.unadj),
        weight1 = c(data.trt.use[, prob.var], rep(0, nrow(data.trt.unadj)))
      )
    }
    res.adj <- rep(0, length(all.std.diff))
    names(res.adj) <- names(all.std.diff)
    res.adj[names(res.adj.raw)] <- res.adj.raw
    if (!is.null(data.in.unadj)) {
      data.trt.use <- data.trt.unadj[which(data.trt.unadj[, indicator.var] %in% c(class.of.int[[i]], -1) & !is.na(data.trt.unadj[, prob.var])), ] # selected class and unk
      data.ctrl.use <- data.ctrl.unadj[which(!is.na(data.ctrl.unadj[, prob.var])), ] # in PSM, only pts matched to this class have valid prob
    }
    res.unadj.raw <- calc_std_diff(
      vars = vars, data0 = data.ctrl.use, weight0 = rep(1, nrow(data.ctrl.use)),
      data1 = data.trt.use, weight1 = rep(1, nrow(data.trt.use))
    )
    res.unadj <- rep(0, length(all.std.diff))
    names(res.unadj) <- names(all.std.diff)
    res.unadj[names(res.unadj.raw)] <- res.unadj.raw
    wmd <- NULL
    wmd$var <- rep(names(res.adj), 2)
    wmd$absolute_std_diff <- abs(c(res.adj, res.unadj))
    wmd$type <- c(
      rep("adjusted difference", length(res.adj)),
      rep("unadjusted difference", length(res.unadj))
    )
    wmd <- as.data.frame(wmd)
    diff.list[[i]] <- wmd
  }
  if (return.levels) {
    output <- list(diff.list = diff.list, factor.levels.list = factor.levels.list)
  } else {
    output <- diff.list
  }
  return(output)
}


#' Compare weighted and unweighted (naive analysis) standardized difference in plot
#' 
#' @md
#' @param diff.list data list returned by function [std_diff].
#' @param legend.pos legend position: "left", "top", "right", "bottom".
#' @inheritParams km_plot_weight
#' @param xlim.low (`numeric`) lower bound of xlim
#' @param xlim.high (`numeric`) upper bound of xlim 
#' 
#' @import ggplot2 
#' @importFrom  grDevices dev.off png
#' @export
#' @examples 
#' \dontrun{
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
#' ipw_res1 <- ipw_strata(
#'   data.in = clinical_1, formula = indicator ~ BECOG + SEX + BNLR,
#'   indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", trt = "ARM",
#'   class.of.int = list("strata_1" = 1, "strata_2" = 0)
#'  )
#' ipw_diff <- std_diff(
#'   data.in = ipw_res1$data, vars = c("BECOG", "SEX", "BNLR"),
#'   indicator.var = "indicator", trt = "ARM",
#'   class.of.int = list("strata_1" = 1, "strata_2" = 0),
#'   usubjid.var = "SUBJID"
#' )
#' std_diff_plot(ipw_diff)
#' }

std_diff_plot <- function(diff.list, legend.pos = "right",
                          prefix.title = "In strata:", xlim.low = 0, xlim.high = 1) {
  for (i in 1:length(diff.list)) diff.list[[i]]$absolute_std_diff[which(diff.list[[i]]$absolute_std_diff == Inf)] <- NA
  Map(function(x) {
    ggplot(diff.list[[x]], aes(.data$absolute_std_diff, .data$var)) +
      geom_point(aes(color = .data$type)) +
      scale_color_manual(values = c("blue", "red")) +
      ggtitle(paste(prefix.title, x)) +
      theme(
        legend.position = legend.pos,
        text = element_text(size = 16),
      ) +
      scale_x_continuous(
        breaks = round(seq(xlim.low, xlim.high, by = 0.1), 1),
        limits = c(xlim.low, xlim.high)
      )
  }, names(diff.list))
}
