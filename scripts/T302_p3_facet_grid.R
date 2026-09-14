library(ggplot2)
library(dplyr)

setwd(Sys.getenv("TMPDIR"))

data_locations <- list.dirs("/scratch/p319788/CINsim/single_sweep_T302_p3", recursive = FALSE)
plots_dir_path <- paste0("plots")

if (!dir.exists(plots_dir_path)){
  dir.create(plots_dir_path, recursive = TRUE)
}

compiled_data <- NULL

for (loc in data_locations){
    read_path <- file.path(loc, 
                            "tmp", 
                            "results", 
                            "single_sweep_T302_p3_",
                            "full_sweep_T302_p3_results.Rds")
    if (is.null(compiled_data)){
        compiled_data <- readRDS(read_path)
    }else{
        new_data <- readRDS(read_path)
        compiled_data <- rbind(compiled_data, new_data)
    }
}

# convert some vars to factor for the heatmap
compiled_data$division_FCs <- round(compiled_data$division_FCs, 2)
compiled_data$division_FCs <- as.factor(compiled_data$division_FCs)

saveRDS(compiled_data, file.path("/scratch/p319788/CINsim/single_sweep_T302_p3", "combined_results.Rds"))

# gets the location of highest CnFS
red_star <- compiled_data %>%
  filter(new_CnFS == max(new_CnFS)) %>%
  select(division_FCs, pMisseg, surv_FC)

# plotting the facet_grid heatmap of p_misseg, surv_FC, div_FC, viability and CnFS
p1 <- ggplot(compiled_data, aes(x = division_FCs, y = pMisseg)) +
    geom_tile(aes(fill = new_CnFS, alpha = viability), color = "white", lwd = 0.2, linetype = 1) +
    scale_fill_gradient(low = "blue", high = "yellow", name = "CnFS") +
    scale_alpha_continuous(range = c(0.1, 1), name = "Viability") +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                    breaks = 10^(-6:-1)) +
    scale_x_discrete(name = "Division FC", breaks = c(10, 5, 3.33, 2.5, 2, 1.67, 1.43, 1.25, 1.11),
                        limits = unique(compiled_data$division_FCs)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
            axis.text.y = element_text(size = 15),
            axis.title = element_text(size = 15),
            aspect.ratio = 1, panel.grid = element_blank()) +
    facet_grid(vars(surv_FC)) + 
    geom_point(data = red_star, aes(x = division_FCs, y = pMisseg), color = "red",
             size = 3, shape = 4)

ggsave(filename = file.path(plots_dir_path, "facet_grid_T302_p3.pdf"),
        plot = p1,
        width = 7,
        height = 6 * length(data_locations),
        unit = "in",
        limitsize = FALSE) # since we're saving a very long plot
