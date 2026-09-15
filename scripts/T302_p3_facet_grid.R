library(ggplot2)
library(dplyr)

TMP_DIR <- Sys.getenv("TMPDIR")
setwd(TMP_DIR)

data_dirs <- list.dirs(
    file.path(
        Sys.getenv("SCRATCH"),
        "CINsim",
        "T302_p3_single_sweep"),
    recursive = FALSE,
    full.names = TRUE
)

plot_dir <- file.path(
    TMP_DIR,
    "plots"
)

if (!dir.exists(plot_dir)){
  dir.create(plot_dir, recursive = TRUE)
}

results_dir <- file.path(
    TMP_DIR,
    "results"
)

if (!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE)
}

# inits the var to save compiled data in
compiled_data <- NULL

pattern <- "^T302_p3_single_sweep_surv_FC_[0-9]+(?:\\.[0-9]+)?_sim_df\\.Rds$"

for (loc in data_dirs){
    # gets the full path of the (should be singular) file in loc
    read_path <- list.files(
        path = file.path(
            loc,
            "tmp",
            "results"),
        pattern = pattern,
        recursive = TRUE,
        full.names = TRUE
    )

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

saveRDS(
    file.path(
        results_dir,
        "merged_data.Rds"
    )
)

# gets the location of highest CnFS
red_star <- compiled_data %>%
  filter(CnFS == max(CnFS)) %>%
  select(division_FCs, pMisseg, surv_FC)

# plotting the facet_grid heatmap of p_misseg, surv_FC, div_FC, viability and CnFS
p1 <- ggplot(compiled_data, aes(x = division_FCs, y = pMisseg)) +
    geom_tile(aes(fill = CnFS, alpha = viability), color = "white", lwd = 0.2, linetype = 1) +
    scale_fill_gradient(low = "blue", high = "yellow", name = "CnFS") +
    scale_alpha_continuous(range = c(0.1, 1), name = "Viability") +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                    breaks = 10^(-4:-1)) +
    scale_x_discrete(name = "Division FC", breaks = c(10, 5, 3.33, 2.5, 2, 1.67, 1.43, 1.25, 1.11),
                        limits = unique(compiled_data$division_FCs)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
            axis.text.y = element_text(size = 15),
            axis.title = element_text(size = 15),
            aspect.ratio = 1, panel.grid = element_blank()) +
    facet_grid(vars(surv_FC)) + 
    geom_point(data = red_star, aes(x = division_FCs, y = pMisseg), color = "red",
             size = 3, shape = 4)

ggsave(filename = file.path(plot_dir, "T302_p3_facet_grid.pdf"),
        plot = p1,
        width = 7,
        height = 6 * length(data_dirs),
        unit = "in",
        limitsize = FALSE) # since we're saving a very long plot
