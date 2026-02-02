library(CINsim)
library(tidyverse)
library(ggplot2)

# sets wd to location of the script - this is dirty, and I do not 
# recommend it for anyone else following in my footsteps.
setwd(Sys.getenv("TMPDIR"))

ITERATIONS <- 100

if (!dir.exists("plots/figure_1")){
  dir.create("plots/figure_1", recursive = TRUE)
}

# plotting theme
cinsim_theme <- function() {
  theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          aspect.ratio = 1)
}

# set the seed for the file, as we do simulations.
# this makes future reruns reproducible (I think)!
set.seed(42)

# we need to load in a few sample simulation results and plot their
# CnFS score over generations. Since these are old sims, we will also need to
# recalculate their CnFS to include the bounding. 

# recalc the CnFS to the new bounded limits
recalc_CnFS <- function(CnFS_score){
  old_denominator <- 1 / CnFS_score
  new_CnFS <- 1 / (old_denominator + 1)
  return(new_CnFS)
}

# we take the data from plot 2B, as the statement of no evolution is made in that context
BASE_DATA_PATH = Sys.getenv("DATA_SOURCE")

# define the data we want to use for the plot, 2 p_missegs and 2 surv_FCs
p_missegs <- list("1e-06", "0.0001743")
surv_FCs <- list("2.14", "7.5")
to_load <- expand.grid(misseg = p_missegs, surv_FC = surv_FCs)

# composes the loading path so that we can loop over this structure in a sec
to_load$loading_path <- file.path(BASE_DATA_PATH,
                                  paste0("misseg_", to_load$misseg,
                                         "_surv_FC_", to_load$surv_FC, ".Rds"))

for (i in 1:nrow(to_load)){
  row <- to_load[i,]
  simulation_obj <- readRDS(row$loading_path)  
  
  CnFS_of_sims <- sapply(1:ITERATIONS, function(i) {
    simulation_obj[[i]]$gen_measures$CnFS # since the old sims didn't track the CnFS I'm rerunning them
    # we can run this script afterwards to make the supplemental figure that shows the lack of evolution
  })

  names(CnFS_of_sims) <- 1:ITERATIONS

  # now we've got the sims we need to extract the CnFS for each
  df <- purrr::map_df(names(CnFS_of_sims), ~{
    data.frame(
      Generations = seq_along(CnFS_of_sims[[.x]]),
      old_CnFS = CnFS_of_sims[[.x]],
      rep = .x
    )
  
  df$new_CnFS <- recalc_CnFS(df$old_CnFS)

  #TODO maybe make a barplot showing the average number of generations per combination of parameters,
  # then show the plot of CnFS progression over simulations?

  #TODO
  # now want to plot that CnFS progression of the sims by plotting the generations against the CnFS
  p1 <- ggplot(df, aes(x = Generations, y = new_CnFS, colour = rep)) +
      geom_line(size = 1.2) +
      geom_point() +
      theme_minimal() #TODO add a title + ggtitle()
  })
}
