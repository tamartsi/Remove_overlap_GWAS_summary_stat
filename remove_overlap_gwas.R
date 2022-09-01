## Goal of function: A function that takes GWAS summary statistics from a meta-analysis (beta_all, se_all), 
## and GWAS summary stat from a study that contributed to this meta-analysis (beta_rm, se_rm), and computes
## what would be the meta-analysis results without the participation of the study. 
##
## Note vectorized form: The function uses vectorized data to perform computation on many 
## variants at the same time. 
##
## QC: currently no variant ID is given, it is assumed that the data are cleaned prior to using this function
## therefore no checks for matches/correct order between the vectors corresponding to tested variants, 
## and not checks for missing values. 

remove_overlap_gwas <- function(beta_all, se_all, beta_rm, se_rm){
  
  # we need variances rather than SEs	
  v_all <- se_all^2
  v_rm <- se_rm^2
  
  
  # sometimes it is convenient to use the inverse-variance weights 
  # rather than the variance: compute weights of study to remove
  w_rm <- 1/v_rm
  
  # compute variances of beta estimates of study "to keep" after removing 
  # "study to remove"
  v_keep <- v_all/(1-v_all/v_rm)
  se_keep <- sqrt(v_keep)
  
  # the inverse-variance weights of study "to keep":
  w_keep <- 1/v_keep
  
  # compute betas after removing overlap:
  beta_keep <- (beta_all*(w_rm + w_keep) - w_rm*beta_rm)/w_keep
  
  # return a data.frame
  res_keep <- data.frame(beta_keep, se_keep)
  res_keep
  
}
