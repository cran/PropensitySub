#' Weighted KM plot

#' @param data.in input data, patients in rows and variables in columns.
#' This could be an output from ipw_strata() or ps_match_strata(). 
#' @param prob.names column names for the probability scores to be used as weights.
#' The order of probnames should match the order of class.of.int.
#' if probnames is NULL, the function will assume that the probnames are
#' pred0, pred1, prod2, prod1or2 in the example in class.of.int.
#' @param filename if it is not NULL, a png file will be generated
#' @param time.unit time unit to be marked in x axis
#' @param prefix.title prefix for title
#' @inheritParams ipw_strata
#' @return return a list of plots with class of `ggsurvplot`
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
#'  km_plot_weight(ipw_res1$data,   
#'    indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", 
#'    trt = "ARM", class.of.int = list("strata_2" = 0))
#' ps_res1 <- ps_match_strata(
#'   data.in = clinical_1, formula = indicator ~ BECOG + SEX + BNLR,
#'   indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", trt = "ARM",
#'   class.of.int = list("strata_1" = 1, "strata_2" = 0)
#'  )
#'  km_plot_weight(ps_res1$data,   
#'    indicator.var = "indicator", tte = "OS_MONTH", event = "OS_EVENT", 
#'    trt = "ARM", class.of.int = list("strata_1" = 1, "strata_2" = 0)) 
#'    
#' @import survminer
#' 
#' @export
#' 

km_plot_weight <- function(data.in, indicator.var = "indicator", class.of.int = NULL,
                           prob.names = NULL, filename = NULL, tte = "AVAL", event = "event",
                           trt = "trt", time.unit = "month", prefix.title = "In strata:") {
  data.trt.name <- unique(data.in[which(data.in[, indicator.var] != -1), trt])
  stopifnot(length(data.trt.name) == 1)
  data.trt <- data.in[which(data.in[, trt] == data.trt.name), ] # take trt only
  data.ctrl <- data.in[which(data.in[, trt] != data.trt.name), ] # take control only

  n.class <- length(unique(data.trt[, indicator.var]))
  stopifnot(n.class >= 2)
  stopifnot(max(data.trt[, indicator.var]) == n.class - 1) # need to be consequent 0,1,2...
  stopifnot(min(data.trt[, indicator.var]) == 0) # need to start from 0
  if (!is.null(class.of.int)) stopifnot(class(class.of.int) == "list")
  if (is.null(class.of.int)) class.of.int <- as.list(sort(unique(data.trt[, indicator.var])))
  if (is.null(names(class.of.int))) names(class.of.int) <- sapply(class.of.int, function(i) paste0(i, collapse = ",")) # those will be row names to show in output table
  if (is.null(prob.names)) prob.names <- sapply(class.of.int, function(i) paste0("pred", paste0(i, collapse = "or")))

  g <- vector("list", length(class.of.int))
  for (i in 1:length(class.of.int)) {
    grp_label <- names(class.of.int)[i]
    prob.var <- prob.names[i]
    data.use <- data.in
    # data.use <- data.in[which(data.in[,indicator.var] %in% c(-1, class.of.int[[i]])),] # take all
    # the ctrl pts, and class i pts in trt arm
    data.use[, c("tte", "event", "trt")] <- data.use[, c(tte, event, trt)]
    data.use <- data.use[which(!is.na(data.use[, prob.var])), ] # take out pts with missing weights
    data.use <- data.use[which(data.use[, prob.var] != 0), ] # take out pts with weight = 0
    # (https://github.com/kassambara/survminer/issues/382)
    data.use$trt <- paste(grp_label, data.use$trt)
    km.fit <- survfit(Surv(tte, event) ~ trt,
      weights = data.use[, prob.var], data = data.use
    )
    g[[i]] <- ggsurvplot(km.fit,
      data = data.use, risk.table = FALSE,
      pval = FALSE, conf.int = FALSE, axes.offset = FALSE,
      xlim = c(0, max(data.use[, tte])), break.time.by = 3,
      xlab = paste0("Time (", time.unit, ")"),
      palette = c("blue", "red"),
      legend.title = "",
      ggtheme = theme_survminer(
        base_size = 16,
        font.x = c(16, "plain", "black"),
        font.y = c(16, "plain", "black"),
        font.caption = c(18, "plain", "black"),
        font.tickslab = c(16, "plain", "black"),
        font.legend = c(14, "plain", "black")
      ),
      title = paste0(prefix.title, grp_label)
    )


    if (!is.null(filename)) png(paste(filename, "_", grp_label, ".png", sep = ""), width = 720, height = 480)
    # print(g[[i]])
    if (!is.null(filename)) dev.off()
  }

  g
}
