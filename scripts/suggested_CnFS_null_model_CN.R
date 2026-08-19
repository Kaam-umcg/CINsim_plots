# testing a possible null-model for the CnFS metric.
# main concern currently is that the CnFS metric represents an abstract
# value without an expectated value, making it non-interpretable against
# background. To address this issue, we need to make a H0-model, or an 
# expected value distribution for the CnFS to take that we can compare
# the acquired value against.
library(CINsim)
`%>%` <- magrittr::`%>%`
scratch_dir <- Sys.getenv("SCRATCH")
if (scratch_dir == ""){
    print("Did not extract scratch path - exiting")
    q("n")
}

save_dir <- file.path(scratch_dir, "CINsim", "null_model_test_modal_CN_permutation")
plot_dir <- file.path(save_dir, "plots")
output_dir <- file.path(save_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# this one doesnt need recur since the parent was just made
dir.create(plot_dir, showWarnings = FALSE, recursive = FALSE)

sims_to_shuffle <- list.files(
    path = "/scratch/p319788/CINsim/best_sims",
    full.names = TRUE)

names(sims_to_shuffle) <- tools::file_path_sans_ext(basename(sims_to_shuffle))

create_modal_karyotype <- function(selection_metric, n = 1000, p_perm = 0.0001){
    # here we define the modal karyotype of the selection metric
    modal_karyotype <- apply(
        X = selection_metric,
        MARGIN = 2, # apply over columns
        FUN = function(col){
            idx <- which.max(col)
            rn <- rownames(selection_metric)[idx]
            return(as.integer(rn))
        }
    )

    modal_cells <- matrix(
        data = modal_karyotype,
        nrow = n,
        ncol = ncol(selection_metric),
        byrow = TRUE)
    colnames(modal_cells) <- colnames(selection_metric)
    rownames(modal_cells) <- paste0("cell_", 1:n)

    # now we want to permutate CNs a little bit
    permutated_karyos <- apply(
        X = modal_cells,
        MARGIN = 1, 
        FUN = swap_two_in_row
    )
    # transpose since apply flips the matrix
    permutated_karyos <- t(permutated_karyos)

    return(permutated_karyos)
}

swap_two_in_row <- function(row) {
  # row: numeric (or character) vector representing one row of a matrix
  
  n_cols <- base::length(row)
  
  # If fewer than 2 columns, nothing to swap
  if (n_cols < 2) {
    return(row)
  }
  
  # Randomly choose two distinct positions
  idx <- base::sample.int(n_cols, size = 2, replace = FALSE)
  
  # Swap values at these positions
  tmp <- row[idx[1]]
  row[idx[1]] <- row[idx[2]]
  row[idx[2]] <- tmp
  
  return(row)
}

modal_shuffle_CnFS <- function(sim){
    selection_metric <- sim$selection_metric
    shuffled_modal_karyos <- create_modal_karyotype(
        selection_metric = selection_metric,
        n = 1000,
        p_perm = 0.0001
    )
    CnFS <- calc_CnFS(
        karyotypes = shuffled_modal_karyos,
        selection_metric = selection_metric
    )
    return(CnFS)
}

for (sim_to_shuffle_name in names(sims_to_shuffle)){
    sim_to_shuffle <- readRDS(sims_to_shuffle[sim_to_shuffle_name])
    # recalculates the final CnFS of all simulations, since not all used the new CnFS mapping
    CnFS_of_sims <- unlist(lapply(sim_to_shuffle, FUN = function(sim){
        CnFS <- calc_CnFS(
            karyotypes = sim$karyotypes,
            selection_metric = sim$selection_metric)
        return(CnFS)
    }))
    mean_sim_CnFS <- mean(CnFS_of_sims)
    sim_CnFS <- data.frame(CnFS_of_sims)
    colnames(sim_CnFS) <- "value"

    shuffled_modal_CnFS_values <- lapply(
        X = sim_to_shuffle, 
        FUN = modal_shuffle_CnFS
    )
    # makes the CnFS values into a format that is easier for plotting
    CnFS_combined <- purrr::map_df(
        shuffled_modal_CnFS_values,
        .f = function(x) {
            df <- data.frame(value = x)
            return(df)},
        .id = "source")

    p1 <- ggplot2::ggplot(CnFS_combined, ggplot2::aes(x = value)) +
    ggplot2::geom_density(
        mapping = ggplot2::aes(y = ggplot2::after_stat(density / max(density)), colour = "Permutated"),
        alpha = 0.5,
        fill = "lightblue",
        data = CnFS_combined) +
    ggplot2::geom_density(
        mapping = ggplot2::aes(y = ggplot2::after_stat(density / max(density)), colour = "Normal"),
        alpha = 0.2,
        fill = "red",
        data = sim_CnFS) +
    ggplot2::geom_vline(
        xintercept = mean_sim_CnFS,
        linetype = "dashed",
        linewidth = 0.5,
        colour = "red") +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::scale_colour_manual(
        name = "Data type",
        values = c(
            "Permutated" = "lightblue",
            "Normal" = "red"
        )
    ) +
    ggplot2::labs(
        x = "CnFS value of simulation",
        y = "Scaled density",
        title = glue::glue("Shuffled modal karyotype CnFS against real CnFS for {sim_to_shuffle_name}")) +
    CINsim:::cinsim_theme()

    plot_name <- file.path(plot_dir, glue::glue("permuted_modal_CnFS_density_", sim_to_shuffle_name,".png"))

    ggplot2::ggsave(filename = plot_name,
                    plot = p1)

    file.copy(from = plot_name,
                to = "~",
                overwrite = TRUE)

}
