## ---- message=FALSE, warning=FALSE, comment=""--------------------------------
library(PropensitySub) 
library(dplyr)
head(biomarker, 3)
#  Control arm subjects have no STRATUM measurments
#  Experimental arm subjects partially miss STRATUM measurement.
table(biomarker$Arm, biomarker$STRATUM, useNA = "ifany")

## ---- message=FALSE-----------------------------------------------------------
biomarker <- biomarker %>%
  mutate(
    indicator_1 = case_when(
      STRATUM == "Positive" ~ 1,
      STRATUM == "Negative" ~ 0,
      # impute missing as "Negative" in Experimental Arm
      Arm == "Experimental" & is.na(STRATUM) ~ 0,
      is.na(STRATUM) & Arm == "Control" ~ -1
    ),
    indicator_2 = case_when(
      STRATUM == "Positive" ~ 1,
      STRATUM == "Negative" ~ 0,
      # keep missing as a seperate stratum in Experimental Arm
      Arm == "Experimental" & is.na(STRATUM) ~ 2,
      is.na(STRATUM) & Arm == "Control" ~ -1
    ),
    # treatment group needs to be factor
    Arm = factor(Arm, levels = c("Control", "Experimental"))
  
  )

## ---- message=FALSE, comment=""-----------------------------------------------
# `plain` model with IPW method and aggressive imputation where missing is imputated as negative
ipw_plain_str2 <- ipw_strata(
  data.in = biomarker, formula = indicator_1 ~ ECOG + Sex + Age, model = "plain",
  indicator.var = "indicator_1", tte = "OS", event = "OS.event", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0)
)
# get weighted HRs
ipw_plain_str2$stat
# check model converge
ipw_plain_str2$converged
# get weights
ipw_plain_str2$data %>% 
  dplyr::select(Patient.ID, Arm, STRATUM, ECOG, Sex, Age, indicator_1, pred1, pred0) %>%
  head()


# `plain` model with IPW method and missing is kept as another stratum
ipw_plain_str3 <- ipw_strata(
  data.in = biomarker, formula = indicator_2 ~ ECOG + Sex + Age, model = "plain",
  indicator.var = "indicator_2", tte = "OS", event = "OS.event", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0, "Unknown" = 2)
)
# get weighted HRs
ipw_plain_str3$stat

# `plain` model with PSM method and aggressive imputation
ps_plain_str2 <- ps_match_strata(
  data.in = biomarker, formula = indicator_1 ~ ECOG + Sex + Age, model = "plain",
  indicator.var = "indicator_1", tte = "OS", event = "OS.event", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0)
)
# get weighted HRs
ps_plain_str2$stat 

## ---- message=FALSE, comment=""-----------------------------------------------
# dwc model with IPW method
ipw_dwc <- ipw_strata(
  data.in = biomarker, formula = indicator_2 ~ ECOG + Sex + Age, model = "dwc",
  indicator.var = "indicator_2", tte = "OS", event = "OS.event", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0)
)
# get weighted HRs
ipw_dwc$stat 

# dwc model with PSM method
ps_dwc <- ps_match_strata(
  data.in = biomarker, formula = indicator_2 ~ ECOG + Sex + Age, model = "dwc",
  indicator.var = "indicator_2", tte = "OS", event = "OS.event", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0, "Missing" = 2)
)
# get weighted HRs
ps_dwc$stat 

## ---- message=FALSE, comment=""-----------------------------------------------
# data process: create a numeric version of next STRATUM to learn from  
biomarker <- biomarker %>%
  mutate(
    indicator_next_2 = case_when(
      STRATUM.next == "Positive" ~ 1,
      STRATUM.next == "Negative" ~ 0, 
      Arm == "Experimental" & is.na(STRATUM.next) ~ 2,
      is.na(STRATUM.next) & Arm == "Control" ~ -1
    )
  )
ipw_wri <- ipw_strata(
  data.in = biomarker, formula = indicator_2 ~ ECOG + Sex + Age, model = "wri",
  indicator.var = "indicator_2", indicator.next = "indicator_next_2",
  tte = "OS", event = "OS.event", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0)
)
# get weighted HRs
ipw_wri$stat 

## ---- message=FALSE, fig.width=6, fig.height=5, comment=""--------------------
# for ipw_plain model results
km_plot_weight(
  data.in = ipw_plain_str2$data,
  indicator.var = "indicator_1", tte = "OS", event = "OS.event", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0)
)

# to get weights from model result for further usage
ipw_plain_str2$data %>% 
  dplyr::select(Patient.ID, Arm, STRATUM, ECOG, Sex, Age, indicator_1, pred1, pred0) %>%
  head()
# for ipw_wri model results
km_plot_weight(
  data.in = ipw_wri$data,
  indicator.var = "indicator_2", tte = "OS", event = "OS.event", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0)
)
# for ps_dwc model results
km_plot_weight(
  data.in = ps_dwc$data,
  indicator.var = "indicator_2", tte = "OS", event = "OS.event", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0, "Missing" = 2)
)


## ---- message=FALSE, fig.width=6, fig.height=5, comment=""--------------------
ipw_plain_diff <- std_diff(
  data.in = ipw_plain_str2$data, vars = c("ECOG", "Sex", "Age"),
  indicator.var = "indicator_1", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0),
  usubjid.var = "Patient.ID"
)
ipw_plain_diff$Positive
ipw_plain_diff$Negative
# Visualize differences 
std_diff_plot(ipw_plain_diff, legend.pos = "bottom")


## ---- message=FALSE, fig.width=6, fig.height=5, comment=""--------------------

ps_dwc_diff <- std_diff(
  data.in = ps_dwc$data, vars = c("ECOG", "Sex", "Age"),
  indicator.var = "indicator_2", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0, "Missing" = 2),
  usubjid.var = "Patient.ID"
)
ps_dwc_diff$Missing
 
# Visualize differences 
std_diff_plot(ps_dwc_diff)

## ---- message=FALSE, comment=""-----------------------------------------------
vars <- c("ECOG", "Sex", "Age")
thresholds <- c(0.15, 0.2)
class_int_list <- list("Positive" = 1, "Negative" = 0)
rand_ratio <- 2 # randomization ratio
n_arms <- lapply(class_int_list, function(x) nrow(biomarker %>% filter(indicator_1 %in% x)))
exp_diff <- sapply(n_arms, function(x) {
  expected_feature_diff(
    n.feature = length(vars),
    n.arm1 = x,
    n.arm2 = x / rand_ratio,
    threshold = thresholds
  )
}) %>% t()

# Expected imbalanced features
exp_diff

# Calculate the observed imbalanced features in model ipw_plain_diff
obs_diff_cnt <- sapply(thresholds, function(th){
  sapply(ipw_plain_diff, function(gp){
    ft <- as.character(subset(gp, type=="adjusted difference" & absolute_std_diff>=th)$var)
    length(unique(sapply(ft, function(ff)strsplit(ff, split="\\.")[[1]][1])))
  })
})
colnames(obs_diff_cnt) <- paste0("Observed # features > ", thresholds)
rownames(obs_diff_cnt) <- names(ipw_plain_diff)
# the number of observed imbalanced features for each threshold
# Compare expected to observed # of imbalanced features to check model fit
obs_diff_cnt

obs_diff_fac <- sapply(thresholds, function(th) {
  sapply(ipw_plain_diff, function(gp) {
    ft <- as.character(subset(gp, type=="adjusted difference" & absolute_std_diff>=th)$var)
    paste(ft, collapse=",")
  })
})
colnames(obs_diff_fac) <- paste0("features > ", thresholds)
rownames(obs_diff_fac) <- names(ipw_plain_diff)
# the observed individual features that are imbalanced for each threshold
obs_diff_fac



## ---- message=FALSE, comment=""-----------------------------------------------
boot_ipw_plain <- bootstrap_propen(
  data.in = biomarker, formula = indicator_1 ~ ECOG + Sex + Age,
  indicator.var = "indicator_1", tte = "OS", event = "OS.event", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0),
  estimate.res = ipw_plain_str2, method = "ipw", n.boot = 100
)
# get bootstrap CI
boot_ipw_plain$est.ci.mat
# summary statistics from bootstraps
boot_ipw_plain$boot.out.est
# error status and convergence status
boot_ipw_plain$error.est
boot_ipw_plain$conv.est

## ---- message=FALSE, comment=""-----------------------------------------------
boot_ipw_wri <- bootstrap_propen(
  data.in = biomarker, formula = indicator_2 ~ ECOG + Sex + Age,
  indicator.var = "indicator_2", indicator.next = "indicator_next_2",
  tte = "OS", event = "OS.event", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0),
  estimate.res = ipw_wri, method = "ipw", n.boot = 100
)
# get bootstrap CI
boot_ipw_wri$est.ci.mat  

## ---- message=FALSE, comment=""-----------------------------------------------
boot_ps_dwc <- bootstrap_propen(
  data.in = biomarker, formula = indicator_2 ~ ECOG + Sex + Age,
  indicator.var = "indicator_2",  
  tte = "OS", event = "OS.event", trt = "Arm",
  class.of.int = list("Positive" = 1, "Negative" = 0),
  estimate.res = ps_dwc, method = "ps", n.boot = 100
)
# get bootstrap CI
boot_ps_dwc$est.ci.mat  

## ---- message=FALSE, fig.width=7, fig.height=6, comment=""--------------------
boot_models <- list(
  ipw_plain = boot_ipw_plain, 
  ipw_wri = boot_ipw_wri, 
  ps_dwc = boot_ps_dwc
)
cols <- c("Estimate", "Bootstrap CI low", "Bootstrap CI high")
# get HRs and bootstrap CIs
boots_dp <- lapply(boot_models, function(x){
  cis <- x$est.ci.mat[ , cols, drop = FALSE] %>%
    exp() %>% round(2) %>% as.data.frame() 
  colnames(cis) <- c("HR", "LOWER", "UPPER")
  cis %>% mutate(
    Group = rownames(cis)
  )
})
boots_dp <- do.call(`rbind`, boots_dp) %>%
  mutate(
    Methods = rep(c("IPW plain", "IPW wri", "PS dwc"), 2),
    Methods_Group = paste(Methods, Group), 
    Group = factor(Group, levels = c("Positive", "Negative")),
    Methods = factor(Methods, levels = c("IPW plain", "IPW wri", "PS dwc"))
  ) %>%
  arrange(Methods, Group)
forest_bygroup(
  data = boots_dp, summarystat = "HR", upperci = "UPPER", lowerci = "LOWER",
  population.name = "Methods_Group", group.name = "Methods",
  color.group.name = "Group", 
  stat.label = "Hazard Ratio", 
  stat.type.hr = TRUE, log.scale = FALSE,  
  endpoint.name = "OS", study.name = "Example Study", draw = TRUE
)

