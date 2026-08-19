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

# we now have functions that can create a perfectly matching cell population, but we need to do
# some permutation to spread the distribution.