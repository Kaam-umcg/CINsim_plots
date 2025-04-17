#' Function for calculating the Copy number Frequency Score (CnFS)
#' Takes a singular karyosim object as input and inputs the CnFS of that
#' Karyosim object
get_CnFS <- function(sim_results) {
    if (!class(sim_results) == "karyoSim"){
    message("Incorrect type for given object, needs to be karyoSim")
    return(FALSE)
    }
  # inits the inverse score to 0, smaller reverse score = better
  inverse_CnFS <- 0
  
  # gets the observed (target) karyotype frequencies
  karyo_obs <- sim_results$selection_metric

  # and the copy number frequencies of the final generation
  final_g <- max(sim_results$pop_measures$g)
  karyo_sim_final <- sim_results$pop_measures[sim_results$pop_measures$g == final_g, ]
  
  # calculates the CnFS (Copy number Frequency Score) for the final generation
  for (chrom in colnames(karyo_obs)){
    # early exit for X, since we skip that chrom for the CnFS
    if (chrom == "X"){
      next
    }

    # cast to int as we need it for comparisons
    chrom <- as.integer(chrom)
    chrom_df <- karyo_sim_final[karyo_sim_final$chromosome == chrom, ]
    
    for (cn in rownames(karyo_obs)){
      freq_obs <- karyo_obs[as.character(cn), as.character(chrom)]
      
      # test whether there is a value for cn in chrom_df,
      # if not we assign it as 0.
      if (!(as.integer(cn) %in% chrom_df$copy)){
        freq_sim <- 0
      }else{
        freq_sim <- chrom_df[chrom_df$copy == as.integer(cn), ]$cn_freq[[1]][["freq"]]
      }
      
      # takes the square of difference between observed and real frequency
      chrom_cn_CnFS <- (freq_sim - freq_obs)^2
      # and adds that to the inverse CnFS term
      inverse_CnFS <- inverse_CnFS + chrom_cn_CnFS
    }
  }
  # calculates the final CnFS from the inverse
  CnFS <- 1 / inverse_CnFS
  return(CnFS)
}