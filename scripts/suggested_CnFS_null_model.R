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

save_dir <- file.path(scratch_dir, "CINsim", "null_model_test")
plot_dir <- file.path(save_dir, "plots")
output_dir <- file.path(save_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# this one doesnt need recur since the parent was just made
dir.create(plot_dir, showWarnings = FALSE, recursive = FALSE)

# loads in a best-in-set simulation to shuffle CNs for
sim_to_shuffle <- readRDS("/scratch/p319788/CINsim/best_sims/full_T_ALL.Rds")

permute_row <- function(row){
    # shuffles the values of the row and returns the result
    permuted_row <- sample(row, size = length(row), replace = FALSE)
    return(permuted_row)
}

shuffle_cn_assignment <- function(cn_matrix, seed = FALSE){
    # function that takes in a copy number matrix as created from CINsim and
    # shuffles the copy number over all the chromosomes. This effectively
    # keeps the aneuploidy the same, 

    # sets seed for reproducibility, if wanted
    if (is.integer(seed)){
        set.seed(seed)
    }
    shuffled_cn_matrix <- apply(X = cn_matrix, 
                                MARGIN = 1, 
                                FUN = permute_row)
    # the matrix comes out transposed, so we transpose it back to cells (rows) x chromosomes (cols)
    shuffled_cn_matrix <- base::t(shuffled_cn_matrix)
    colnames(shuffled_cn_matrix) <- colnames(cn_matrix)
    return(shuffled_cn_matrix)
}

permute_simulation <- function(sim, permutations = 50){
    # takes in a CINsim object and does permutations on the cn_matrix, calculating a CnFS for all permutations
    # and returning those values
    # seeding is currently not taken into account - since it's a permutation test it shouldn't matter
    # too much, but we'll see after some testing
    permuted_CnFS <- vector("numeric", permutations)
    cn_matrix <- sim$karyotypes
    selection_metric <- sim$selection_metric

    # this is an extremely Pythonic way of coding this, possibly better (faster) way to do this
    for (i in seq_len(permutations)){
        permuted_cn_matrix <- shuffle_cn_assignment(cn_matrix = cn_matrix, seed = FALSE)
        permuted_CnFS[i] <- calc_CnFS(karyotypes = permuted_cn_matrix,
                                        selection_metric = selection_metric)
    }

    return(permuted_CnFS)
}

CnFS_permutated <- lapply(sim_to_shuffle, 
                            FUN = permute_simulation)

# makes the CnFS values into a format that is easier for plotting
CnFS_combined <- purrr::map_df(
  CnFS_permutated,
  .f = function(x) {
    df <- data.frame(value = x)
    return(df)
  },
  .id = "source"
)

# gets a p-value by checking how many simulation have a higher CnFS than the mean from the simulation


# now we should be able to make a density plot, highlighting where our actual CnFS maps
p1 <- ggplot2::ggplot(CnFS_combined, ggplot2::aes(x = value)) +
        ggplot2::geom_density() +
        # adds in the actual CnFS for this simulation         
        ggplot2::labs() +
        # triple dot since its a hidden variable :^)
        CINsim:::cinsim_theme()

ggplot2::ggsave(filename = file.path(plot_dir, "permuted_CnFS_density.png"),
                plot = p1)