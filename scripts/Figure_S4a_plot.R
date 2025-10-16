# standard boilerplate code for HPC
library(CINsim)
library(tidyverse)
library(ggplot2)
setwd(Sys.getenv("TMPDIR"))

if (!dir.exists("plots/figure_S4")){
  dir.create("plots/figure_S4", recursive = TRUE)
}
if (!dir.exists("results")){
  dir.create("results", recursive = TRUE)
}

# sets seed for reproducibility 
set.seed(42)

# Gets paths from associated .sh scripts
BASE_DATA_PATH <- Sys.getenv("BASE_PATH")
DIPLOID_PATH <- Sys.getenv("DIPLOID_PATH")
WGD_PATH <- Sys.getenv("WGD_PATH")
OUTPUT_PATH <- Sys.getenv("OUTPUT_PATH")
NUM_CELLS <- Sys.getenv("START_CELLS")

# these should be consistent though!
ORGANOID_NAMES <- c("14T", "16T", "9T", "24TB")

# we defines some min/max ranges for the CnFS colour scale to allow for comparisons 
# between the heatmaps of the same organoids, but keep a unique scale for each organoid.
# we do this because the inter organoid comparisons are much more interesting than between for 
# this context!
scale_lims <- list(
        "9T" = c(0, 2.5), 
        "14T" = c(0, 240), 
        "16T" = c(0, 0.5), 
        "24TB" = c(0, 0.3))

# loops through all our relevant variables
for(load_path in c(DIPLOID_PATH, WGD_PATH)){
        if (grepl("WGD", load_path)){
                sim_type <- "WGD"
        }else{
                sim_type <- "diploid"
        }
#        if(sim_type == "diploid"){next} #for skipping diploid plotting

        for (organoid in ORGANOID_NAMES){
                scale_lim <- scale_lims[[organoid]]
                # compiles and loads the data
                sim_data <- readRDS(file.path(BASE_DATA_PATH, load_path, paste0(organoid, ".Rds")))
                sim_data$division_FC <- factor(sim_data$division_FC)

                # we need to make plots per survival FC
                for (surv_FC in unique(sim_data$survival_FC)){
                        plot_data <- sim_data[sim_data$survival_FC == surv_FC,]

                        p1 <- ggplot(plot_data, aes(x = division_FC, y = p_misseg)) +
                                geom_tile(aes(fill = CnFS), color = "white", lwd = 0.2, linetype = 1) +
                                scale_fill_gradient(low = "blue", high = "yellow", name = "CnFS", limits = scale_lim) +
                                scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                                        breaks = 10^(-6:-1)) +
                                scale_x_discrete(name = "Division FC", breaks = c(10, 5, 3.33, 2.5, 2, 1.67, 1.43, 1.25, 1.11),
                                        limits = unique(plot_data$division_FCs)) +
                                labs(x = "Division FC", y = "p misseg") +
                                theme(aspect.ratio = 1, panel.grid = element_blank()) +
                                coord_fixed()
                        ggsave(filename = file.path("plots/figure_S4", paste0(organoid, "_", sim_type, "_", surv_FC,".pdf")))
                        saveRDS(plot_data, file = file.path("results", paste0(organoid, "_", sim_type, "_", surv_FC,".Rds")))
                        # now that we have the plot, we want to store it before composing the entire patchwork
        }
    }
    # creates the patchwork plot for either WGD or diploid
}