# testing a possible null-model for the CnFS metric.
# main concern currently is that the CnFS metric represents an abstract
# value without an expectated value, making it non-interpretable against
# background. To address this issue, we need to make a H0-model, or an 
# expected value distribution for the CnFS to take that we can compare
# the acquired value against.
rm(list = ls())

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

sample_cn_for_column <- function(frac_vec, cn_states, n_cells) {
  # frac_vec: numeric vector of fractions for one chromosome (length = number of CN states)
  # cn_states: vector with the CN states corresponding to rows (same length as frac_vec)
  # n_cells: number of cells to sample
  
  # Normalise fractions to sum to 1 (safeguard)
  prob_vec <- frac_vec / base::sum(frac_vec)
  
  sampled_cn <- base::sample(
    x = cn_states,
    size = n_cells,
    replace = TRUE,
    prob = prob_vec
  )
  
  return(sampled_cn)
}

generate_cell_population <- function(frac_mat, n_cells, cn_states = NULL, seed = NULL) {
  # frac_mat: matrix, rows = CN states, columns = chromosomes, values = fractions
  # n_cells: number of cells (rows) to generate
  # cn_states: optional vector with CN states; by default rownames(frac_mat)
  # seed: optional RNG seed for reproducibility
  
  if (!base::is.null(seed)) {
    base::set.seed(seed)
  }
  
  n_rows <- base::nrow(frac_mat)
  n_cols <- base::ncol(frac_mat)
  
  # Determine CN states (assumed to match row order)
  if (base::is.null(cn_states)) {
    cn_states <- base::rownames(frac_mat)
    # If rownames are NULL, fall back to 1:n_rows
    if (base::is.null(cn_states)) {
      cn_states <- base::seq_len(n_rows)
    }
  }
  
  # Preallocate output matrix: cells × chromosomes
  out_mat <- base::matrix(
    data = NA_integer_,
    nrow = n_cells,
    ncol = n_cols
  )
  
  # Preserve column names
  base::colnames(out_mat) <- base::colnames(frac_mat)
  
  # Sample per column
  for (j in base::seq_len(n_cols)) {
    frac_vec <- frac_mat[, j]
    
    out_mat[, j] <- sample_cn_for_column(
      frac_vec = frac_vec,
      cn_states = cn_states,
      n_cells = n_cells
    )
  }
    # Ensure entire matrix is integer type
  out_mat <- base::matrix(
    data = base::as.integer(out_mat),
    nrow = n_cells,
    ncol = n_cols,
    dimnames = base::dimnames(out_mat)
  )

  return(out_mat)
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

get_CnFS_shuffled_perfect_pop <- function(sim, n_cells = 1000, n_perms = NULL){
  selection_metric <- sim$selection_metric
  permutated_karyos <- generate_cell_population(
    frac_mat = selection_metric,
    n_cells = n_cells,
    cn_states = NULL,
    seed = NULL)
  n_perms <- stats::rnbinom(
    n = 1, 
    mu = 2,
    size = 5)

  # now we want to permutate CNs a little bit
  #TODO - we might want some extra shuffling here?
  # could do multiple calls,
  while (n_perms > 0){
    permutated_karyos <- apply(
        X = permutated_karyos,
        MARGIN = 1, 
        FUN = swap_two_in_row)

    # apply flips dimensions so we transpose them back  
    permutated_karyos <- t(permutated_karyos)
    n_perms <- n_perms - 1
  }

  shuffled_CnFS <- calc_CnFS(
    karyotypes = permutated_karyos,
    selection_metric = selection_metric)

  return(shuffled_CnFS)
}

# we now have functions that can create a perfectly matching cell population, 
# but we need to do some permutation to spread the distribution.
# Doing that in a loop for every set of simulations
for (sim_to_shuffle_name in names(sims_to_shuffle)){
  sim_to_shuffle <- readRDS(sims_to_shuffle[sim_to_shuffle_name])

  CnFS_of_sims <- unlist(lapply(sim_to_shuffle, FUN = function(sim){
      CnFS <- calc_CnFS(
          karyotypes = sim$karyotypes,
          selection_metric = sim$selection_metric)
      return(CnFS)
  }))

  mean_sim_CnFS <- mean(CnFS_of_sims)
  sim_CnFS <- data.frame(CnFS_of_sims)
  colnames(sim_CnFS) <- "value"

  permutated_CnFS <- lapply(
    X = sim_to_shuffle,
    FUN = get_CnFS_shuffled_perfect_pop)

  # now we want to make the distribution plot for the real CnFS and the permutated ones
  CnFS_combined <- purrr::map_df(
    permutated_CnFS,
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
    ggplot2::xlim(c(0, 1)) +
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
        title = glue::glue("Shuffled perfect karyotype CnFS against real CnFS for {sim_to_shuffle_name}")) +
    CINsim:::cinsim_theme()

  plot_name <- file.path(plot_dir, glue::glue("permuted_perfect_CnFS_density_", sim_to_shuffle_name,".png"))

  ggplot2::ggsave(filename = plot_name,
                  plot = p1)
  file.copy(from = plot_name,
            to = "~",
            overwrite = TRUE)
}