library(CINsim)
library(ggplot2)
library(dplyr)
library(ggbreak)

print(getwd())
# constants for collating the data together
ORGANOID_NAMES <- c("9T", "14T", "16T", "24TB")
SIM_TYPES <- c("diploid", "WGD")
ITERATIONS <- 6

# prepares the separate cols for some more readability
df_organoid_col <- rep(ORGANOID_NAMES, each = ITERATIONS * length(SIM_TYPES))
df_sim_col <- rep(rep(SIM_TYPES, each = ITERATIONS), times = length(ORGANOID_NAMES))
df_repeat_col <- rep(seq(1, ITERATIONS), times = length(ORGANOID_NAMES) * length(SIM_TYPES))
df_CnFS_col <- rep(0, times = ITERATIONS * length(SIM_TYPES) * length(ORGANOID_NAMES))

# creates the dataframe we'll use for plotting
best_CnFS <- data.frame(
  organoid = df_organoid_col,
  sim_type = df_sim_col,
  rep = df_repeat_col,
  CnFS = df_CnFS_col
)

# loads the required data for filling the dataframe with CnFS values
for (organoid in ORGANOID_NAMES){
    for (sim_type in SIM_TYPES){
        # this path management is going to be a bit creative again
        sim_summary_path <- file.path(
            paste0("/scratch/p319788/CINsim/Figure_S4a_", sim_type),
            "tmp",
            "results",
            paste0(organoid, ".Rds"))
        
        sim_summary <- readRDS(sim_summary_path)

        best_row <- sim_summary[which.max(sim_summary$CnFS), ]
        print(organoid)
        print(best_row)

        # we've now got the row with the best params, we can then load the relevant sims
        sim_best_params_fn <- paste0(sim_type, "_", organoid, "_misseg_", round(best_row$p_misseg, 7),
                                    "_div_FC_", best_row$division_FC,
                                    "_surv_FC_", best_row$survival_FC, ".Rds")

        sim_best_params_path <- file.path(
            paste0("/scratch/p319788/CINsim/Figure_S4a_", sim_type),
            "tmp",
            "results",
            sim_best_params_fn
        )
        best_sims <- readRDS(sim_best_params_path)

        # then we can finally start extracting the CnFS from the simulations 
        # and filling it in the df for the barchart in a second
        for (rep in seq(1, ITERATIONS)){
            CnFS_list <- best_sims[[rep]]$gen_measures$CnFS
            CnFS_final <- CnFS_list[[length(CnFS_list)]]

            # only changes the value in the matching row
            best_CnFS[
                best_CnFS$organoid == organoid &
                best_CnFS$sim_type == sim_type &
                best_CnFS$rep == rep, 
                "CnFS"] <- CnFS_final
        }
    }
}

# recalc the CnFS to the new bounded limits
recalc_CnFS <- function(CnFS_score){
  old_denominator <- 1 / CnFS_score
  new_CnFS <- 1 / (old_denominator + 1)
  return(new_CnFS)
}

# recalculates the CnFS to be bounded between 0-1
best_CnFS$CnFS <- recalc_CnFS(best_CnFS$CnFS)

# summarize your data by organoid and sim_type
df_summary <- best_CnFS %>%
  group_by(organoid, sim_type) %>%
  summarise(
    mean_CnFS = mean(CnFS, na.rm = TRUE),
    sd_CnFS = sd(CnFS, na.rm = TRUE),
    n = n(),
    se_CnFS = sd_CnFS / sqrt(n)
  )

# plot low limit
p1 <- ggplot(df_summary, aes(x = organoid, y = mean_CnFS, fill = sim_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean_CnFS - se_CnFS, ymax = mean_CnFS + se_CnFS),
                position = position_dodge(width = 0.9),
                width = 0.2) +
  labs(
    title = "CnFS by organoid and simulation type",
    y = "Mean CnFS (± SE)",
    x = "Organoid"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("WGD" = "#f8766d", "diploid" = "#00bfc4"))

plot_fn <- "organoid_CnFS_comp.pdf"
ggsave(file.path("plots/figure_5", plot_fn), plot = p1, useDingbats = FALSE, create.dir = TRUE)
