calculate_pval_randomizations <- function(obs, ran) {
  # Calculating p.value as it would be calculated using permTest from regioneR (version 1.6.2)
  orig.ev <- sum(obs)
  rand.ev <- ran
  num.nas <- length(which(is.na(rand.ev)))
  ntimes  <- length(rand.ev)
  alternative <-  "auto"
  num.valid.values <- ntimes - num.nas
  if (num.valid.values < ntimes) {
    if (num.valid.values > 0) {
      warning(paste0(num.nas, " iterations returned NA or NaN. Only ", 
                     , " iterations have been used to compute the p-value."))
    }
    else {
      warning(paste0("All ", num.nas, " iterations returned NA or NaN. No valid values returned. It is not possible to compute the p-value nor z-score."))
    }
  }
  if (num.valid.values > 0) {
    if (alternative == "auto") {
      alt <- ifelse(orig.ev < mean(rand.ev, na.rm = TRUE), 
                    "less", "greater")
    }
    else {
      alt <- alternative
    }
    if (alt == "less") {
      pval <- (sum(orig.ev >= rand.ev, na.rm = TRUE) + 
                 1)/(num.valid.values + 1)
    }
    else {
      pval <- (sum(orig.ev <= rand.ev, na.rm = TRUE) + 
                 1)/(num.valid.values + 1)
    }
    if (alternative == "greater" & orig.ev < mean(rand.ev, 
                                                  na.rm = TRUE)) 
      message("Alternative is greater and the observed statistic is less than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
    if (alternative == "less" & orig.ev > mean(rand.ev, 
                                               na.rm = TRUE)) 
      message("Alternative is less and the observed statistic is greater than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
    if (orig.ev == 0 & all(rand.ev == 0)) {
      warning(paste0("All permuted values and the original evaluation value are equal to 0. Z-score cannot be computed."))
      pval <- 1
      zscore <- NA
    }
    else {
      zscore <- round((orig.ev - mean(rand.ev, na.rm = TRUE))/sd(rand.ev, 
                                                                 na.rm = TRUE), 4)
    }
    
  } else {
    pval <- NA
    zscore <- NA
    alt <- alternative
  }
  if (!is.na(pval)) {
    if (!is.finite(zscore)) {
      warning(paste0("All permuted values are equal to ", 
                     rand.ev[1], ". Z-score is infinite."))
    }
  }
  res <- list(pval = pval, 
              ntimes = ntimes, 
              alternative = alt, 
              observed = orig.ev, 
              permuted = rand.ev, 
              zscore = zscore) 
  res
}