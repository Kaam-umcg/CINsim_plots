library(ggplot2)
library(patchwork)
library(CINsim)
library(purrr)
library(tidyr)
library(dplyr)
library(pheatmap)

# plotting theme and scale colours for copy numbers
# copy number colors
copy_num_cols <- c("gray90", "darkorchid3", "springgreen2", "red3", "gold2", "navy",
                   "lemonchiffon", "dodgerblue", "chartreuse4")
names(copy_num_cols) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")

setwd(Sys.getenv("TMPDIR"))

if (!dir.exists("plots/figure_S2")){
  dir.create("plots/figure_S2", recursive = TRUE)
}

# some vars used in functions later
set.seed(42)
chrom_names <- c(1:19, "X")

# and for simulations
ITERATIONS <- 6
GENERATIONS <- 250
MAX_NUM_CELLS <- Inf # want the full # of generations possible, so cells are Inf
CORES_AVAIL <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
# and plotting - these should be derived in the script but updating the
# formulas for calculating these atm, and want to keep it consistent
# if you read this as a reviewer - good find! I'll fix it if you name it 
# as an issue.
BASELINE_KMS <- 0.08493877  # from diploid 
BASELINE_CNFS <- 0.07890889 # from diploid

make_CN_heatmap <- function(karyo, chrom_names = c(1:19, "X")){
  karyo <- data.frame(karyo)
  colnames(karyo) <- chrom_names
  karyo$Cell <- paste0("Cell", c(1:length(rownames(karyo))))
  
  karyo <- pivot_longer(karyo, cols = all_of(chrom_names),
                               names_to = "Chromosome", values_to = "CN")
  karyo$Chromosome <- factor(karyo$Chromosome, levels = chrom_names)
  
  p <- ggplot(karyo, aes(x = Chromosome, y = Cell, fill = factor(CN))) +
    geom_tile() +
    scale_fill_manual(values = copy_num_cols) +
    #  scale_fill_manual(values = copy_num_cols) + 
    labs(x = "Chromosome", y = "Cell", fill = "Copy Number") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) 
  return(p)
}

make_pref_gains <- function(n_cells = 50, cn_probs = CINsim::Mps1_X){
  # inits the object where we store all the CNs of the cells
  karyo_mat <- matrix(data = NA, 
                      ncol = length(colnames(cn_probs)), nrow = n_cells, 
                      dimnames = list(c(paste0("cell_", 1:n_cells)), 
                                      c(1:(length(colnames(cn_probs)) - 1), "X")))
  
  #TODO redo this for loop into something vectorised
  # I'm bad at R, I write it like Python :(
  for (i in 1:n_cells){
    for (chrom in colnames(cn_probs)){
      # woe be upon whoever named the cumulative sum function
      cum_probs <- cumsum(cn_probs[, chrom])
      
      # takes a random value and gets the associated CN value from the
      # probability matrix, then assigns it to our new karyo_mat
      random_value <- runif(1)
      selected_CN <- as.integer(names(which(cum_probs >= random_value)[1]))
      karyo_mat[i, chrom] <- selected_CN
      }
  }
  return(karyo_mat)
}

calc_KMS_between_karyos <- function(karyo_1, karyo_2){
  # we use triple colon indexing for the aneuploidy & heterogeneity functions
  # since those are included in CINsim, based on the AneuFinder formulas
  aneu_1 <- CINsim:::qAneuploidy(karyo_1)
  aneu_2 <- CINsim:::qAneuploidy(karyo_2)
  het_1 <- CINsim:::qHeterogeneity(karyo_1)
  het_2 <- CINsim:::qHeterogeneity(karyo_2)
  
  # inits the inverse KMS, since order matters for inverting
  inverse_KMS <- 0
  
  #TODO this can easily be vectorised
  for (chrom in names(aneu_1)){
    # skipping X for KMS calc, as we do in the paper
    if (chrom == "X"){next} 
    # inverse sum of squares of difference for aneuploidy & heterogeneity
    inverse_KMS <- inverse_KMS + (aneu_1[chrom] - aneu_2[chrom])^2 + 
                                  (het_1[chrom] - het_2[chrom])^2
  }
  KMS <- 1 / inverse_KMS
  return(KMS)
}

calc_CnFS_between_karyos <- function(karyo_1, karyo_2, chrom_names = c(1:19, "X")){
  # we make the assumption that the karyos have the same amount of chromosomes
  # ensures colnames are consistent
  karyo_1 <- data.frame(karyo_1)
  colnames(karyo_1) <- chrom_names
  karyo_2 <- data.frame(karyo_2)
  colnames(karyo_2) <- chrom_names
  
  # inits the ivnerse CnFS for calcing the actual CnFS later (order matters)
  inverse_CnFS <- 0
  
  # goes over every chromosome  
  for (chrom in chrom_names){
    # skipping X for CnFS calc, as we do in the paper
    if (chrom == "X"){next}
    
    # simplifies our karyos to a single chrom and gets their frequencies
    karyo_1_prop <- proportions(table(karyo_1[, chrom]))
    karyo_2_prop <- proportions(table(karyo_2[, chrom]))
    
    # then calculates the inverse CnFS for every unique CN
    for (cn in union(unique(names(karyo_1_prop)), unique(names(karyo_2_prop)))){
      # gets the designated ratio or 0 otherwise
      freq_1 <- pluck(karyo_1_prop, cn, .default = 0)
      freq_2 <- pluck(karyo_2_prop, cn, .default = 0)
      
      # applies the inverse CnFS formula, invese sum of squares
      inverse_CnFS <- inverse_CnFS + (freq_1 - freq_2)^2
    }
  }
  CnFS <- 1 / inverse_CnFS
  return(CnFS)
}
  
# creates 7 random karyotypes
random_karyos <- rep(50, 7) %>%
  map(makeKaryotypes, species = "mouse", copies = 2, random = TRUE)
# and their plots
random_karyos_plots <- random_karyos %>%
  map(make_CN_heatmap)

# makes the 7 preferential gain karyotypes, which we sample from 
# the Mps1_X object embedded into the CINsim package
preferential_karyos <- rep(50, 7) %>%
  map(make_pref_gains, cn_probs = CINsim::Mps1_X)
preferential_karyos_plots <- preferential_karyos %>%
  map(make_CN_heatmap)

############
# TODO
# amount of cells doesn't line up for: T158, T257 missing 84 cells total
# Not yet entirely resolved, might want to remake the scWGS karyos.
# Know exactly where to find all the relevant information now
############

# reads Mps1 karyos from saved object
Mps1_karyos_per_cell <- readRDS("data/chrom_counts_per_cell.Rds")

# identifies the unique samples from the loaded object
Mps1_names <- strsplit(as.character(Mps1_karyos_per_cell$sample), split = "_")
Mps1_names <- sapply(Mps1_names, function(x) x[[1]])
Mps1_karyos_per_cell$sample <- Mps1_names
Mps1_karyos <- vector("list", length = length(unique(Mps1_karyos_per_cell$sample)))
Mps1_tumor_count <- 1

# splits the unique tumours into their individual items in a list
for (tumor_name in unique(Mps1_karyos_per_cell$sample)){
  Mps1_karyos[[Mps1_tumor_count]] <- Mps1_karyos_per_cell[Mps1_karyos_per_cell$sample == tumor_name, ]
  Mps1_karyos[[Mps1_tumor_count]]$sample <- NULL
  colnames(Mps1_karyos[[Mps1_tumor_count]]) <- chrom_names
  rownames(Mps1_karyos[[Mps1_tumor_count]]) <- paste0("cell_", 1:length(rownames(Mps1_karyos[[Mps1_tumor_count]])))
  Mps1_tumor_count <- Mps1_tumor_count + 1
}
# makes the CN plots - might look a little funky as they're not 
# 50 cells on the dot.
Mps1_karyos_plots <- Mps1_karyos %>%
  map(make_CN_heatmap)

names(Mps1_karyos) <- unique(Mps1_karyos_per_cell$sample)

# puts all the karyos and plots into a single variable
# thinking of, would be better to do the mapping of karyos into plots
# here in a single line. Oops!
entries <- paste0(rep(c("random", "preferential", "Mps1"), each = 7), "_", rep(1:7, 3))
all_karyos <- append(random_karyos, append(preferential_karyos, Mps1_karyos))
rm(random_karyos, preferential_karyos, Mps1_karyos)
all_plots <- append(random_karyos_plots, append(preferential_karyos_plots,
                                                Mps1_karyos_plots))
rm(random_karyos_plots, preferential_karyos_plots, Mps1_karyos_plots)
names(all_karyos) <- entries
names(all_plots) <- entries

# creates the first plot  
row_title <- vector("list", length = 3)

row_title[[1]] <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Random", 
           size = 6, angle = 0) + theme_void()
row_title[[2]] <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Preferential\ngains", 
           size = 6, angle = 0) + theme_void()
row_title[[3]] <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Mps1\ntumours", 
           size = 6, angle = 0) + theme_void()

# continues with making the patchwork for the first plot
patch_p <- row_title[[1]] +
  all_plots[[1]] + all_plots[[2]] + all_plots[[3]] + all_plots[[4]] + 
  all_plots[[5]] + all_plots[[6]] + all_plots[[7]] +
  row_title[[2]] + 
  all_plots[[8]] + all_plots[[9]] + all_plots[[10]] + all_plots[[11]] + 
  all_plots[[12]] + all_plots[[13]] + all_plots[[14]] +
  row_title[[3]] + 
  all_plots[[15]] + all_plots[[16]] + all_plots[[17]] + all_plots[[18]] +
  all_plots[[19]] + all_plots[[20]] + all_plots[[21]] +
  plot_layout(ncol = 8, nrow = 3, width = c(2, rep(1, 7)),
                             guides = "collect", axes = "collect")
ggsave("plots/figure_S2/fig_s2a.pdf", plot = patch_p)

# Continues with plot S2b, which requires we do pairwise comparisons between
# all the karyotypes
# inits the frame for saving all the CnFS comparisons - we also do this with
# the KMS
matrix_CnFS <- data.frame(matrix(NA, nrow = 21, ncol = 21))
matrix_KMS <- data.frame(matrix(NA, nrow = 21, ncol = 21))

# can 100% do this more efficiently
colnames(matrix_CnFS) <- entries
rownames(matrix_CnFS) <- entries
colnames(matrix_KMS) <- entries
rownames(matrix_KMS) <- entries

# fills the distance matrix with the CnFS distances
pairing <- combn(entries, 2)

# I use the `a` variable to mute the output here
a <- apply(pairing, 2, function(idx) {
  i <- idx[1]
  j <- idx[2]
  # first we do the CnFS
  score <- calc_CnFS_between_karyos(all_karyos[i], all_karyos[j])
  matrix_CnFS[i, j] <<- score  # double arrow due to function scope from apply
  matrix_CnFS[j, i] <<- score # the score is mirrored, so we assign both

  # also sets the diagonal to 0, this does double work but makes it easier later
  matrix_CnFS[i, i] <<- 0
  matrix_CnFS[j, j] <<- 0
  
  # and then the KMS
  # these have the double indexing required `[[`, oopsie
  score <- calc_KMS_between_karyos(all_karyos[[i]], all_karyos[[j]])
  matrix_KMS[i, j] <<- score
  matrix_KMS[j, i] <<- score
  # and diagonal to 0
  matrix_KMS[i, i] <<- 0
  matrix_KMS[j, j] <<- 0
})
rm(a) # deletes the temp var a

# log2 for prettier plotting
dist_CnFS <- log2(matrix_CnFS)
dist_KMS <- log2(matrix_KMS)

# assigns the group factor for plot dressing
group_row <- factor(rep(c("random", "preferential", "Mps1"), each = 7))
row_annotation <- data.frame(Group = group_row)
rownames(row_annotation) <- rownames(dist_CnFS)
col_annotation <- row_annotation

# resets the diagonals to 0 as the log2 makes them -Inf
for (i in 1:nrow(dist_CnFS)){
  dist_CnFS[i, i] <- 0  
  dist_KMS[i, i] <- 0
}
# CnFS plot
p1 <- pheatmap(dist_CnFS,
         color = colorRampPalette(c("blue", "yellow", "red"))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_row = row_annotation,
         annotation_col = col_annotation,
         annotation_colors = list(Group = c("random" = "blue", "preferential" = "green", "Mps1" = "red")),
         scale = "none",
         treeheight_row = 50,
         treeheight_col = 50,
         legend = TRUE)

# KMS plot
p2 <- pheatmap(dist_KMS,
               color = colorRampPalette(c("blue", "yellow", "red"))(100),
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               show_rownames = FALSE,
               show_colnames = FALSE,
               annotation_row = row_annotation,
               annotation_col = col_annotation,
               annotation_colors = list(Group = c("random" = "blue", "preferential" = "green", "Mps1" = "red")),
               scale = "none",
               treeheight_row = 50,
               treeheight_col = 50,
               legend = TRUE)

# patchwork doesn't work for this as these aren't ggplots, so we save to disk
ggsave(filename = "plots/figure_S2/fig_s2b1.pdf", plot = p1)
ggsave(filename = "plots/figure_S2/fig_s2b2.pdf", plot = p2)

# figure S2c
# sets up the dataframe we'll use for plotting the figure later
plotting_frame <- data.frame(
  sim_id = NA,
  selection = rep(c("no", "yes"), each = GENERATIONS * ITERATIONS),
  iteration = rep(rep(1:ITERATIONS, each = GENERATIONS), times = 2),
  CnFS = NA,
  KMS = NA,
  g = rep(1:GENERATIONS, times = ITERATIONS * 2))

# first we do the simulation without selection
no_selection_sims <- parallelCinsim(iterations = ITERATIONS,
                                    g = GENERATIONS,
                                    max_num_cells = MAX_NUM_CELLS,
                                    cores = CORES_AVAIL,
                                    selection_metric = Mps1_X,
                                    selection_mode = NULL,
                                    karyotypes = makeKaryotypes(numCell = 1),
                                    euploid_ref = 2,
                                    pMisseg = 0.0025,
                                    CnFS = TRUE,
                                    KMS = TRUE,
                                    final_aneu_het_scores = karyotype_measures)

# then with selection
selection_sims <- parallelCinsim(iterations = ITERATIONS,
                                 g = GENERATIONS,
                                 max_num_cells = MAX_NUM_CELLS,
                                 cores = CORES_AVAIL,
                                 selection_metric = Mps1_X,
                                 selection_mode = "cn_based",
                                 karyotypes = makeKaryotypes(numCell = 1),
                                 euploid_ref = 2,
                                 pMisseg = 0.0025,
                                 CnFS = TRUE,
                                 KMS = TRUE,
                                 final_aneu_het_scores = karyotype_measures)

# add all the CnFS and KMS values to the plotting frame
# we skip the first val and add 1 to skip the initialisation value (euploid baseline)
CnFS_no_select <- map(no_selection_sims, ~ .x$gen_measures$CnFS[2:(GENERATIONS + 1)])
CnFS_select <- map(selection_sims, ~ .x$gen_measures$CnFS[2:(GENERATIONS + 1)])
plotting_frame$CnFS <- append(unlist(CnFS_no_select), unlist(CnFS_select))

KMS_no_select <- map(no_selection_sims, ~ .x$gen_measures$KMS[2:(GENERATIONS + 1)])
KMS_select <- map(selection_sims, ~ .x$gen_measures$KMS[2:(GENERATIONS + 1)])
plotting_frame$KMS <- append(unlist(KMS_no_select), unlist(KMS_select))

# determines the highest CnFS & KMS from the selection repeats and marks them
final_CnFS <- unlist(map(selection_sims, ~ .x$gen_measures$CnFS[GENERATIONS + 1]))
highest_CnFS_sim <- names(sort(final_CnFS, decreasing = TRUE)[1])

final_KMS <- unlist(map(selection_sims, ~ .x$gen_measures$KMS[GENERATIONS + 1]))
highest_KMS_sim <- names(sort(final_KMS, decreasing = TRUE)[1])

# takes the second highest sim as high_KMS if it's the same as high_CnFS
if (highest_KMS_sim == highest_CnFS_sim){
  highest_KMS_sim <- names(sort(final_KMS, decreasing = TRUE)[2])
}

# assigns the sim_ids for plotting and converts to factor
high_CnFS_number <- strsplit(highest_CnFS_sim, "_")[[1]][2]
high_KMS_number <- strsplit(highest_KMS_sim, "_")[[1]][2]

plotting_frame$sim_id <- "Other"
plotting_frame[plotting_frame$selection == "yes" &
                 plotting_frame$iteration == high_CnFS_number, "sim_id"] <- "High CnFS"

plotting_frame[plotting_frame$selection == "yes" &
                 plotting_frame$iteration == high_KMS_number, "sim_id"] <- "High KMS"

# factor for plotting
plotting_frame$sim_id <- factor(plotting_frame$sim_id)

# because of a bug in CINsim the first KMS value can be Inf in the case of a 
# perfect diploid cell, causing to divide by 0 which in R leads to Inf 
# we clean those values here if they occurred with the baseline value
plotting_frame[plotting_frame$KMS == Inf, "KMS"] <- BASELINE_KMS

# KMS progression plot
p1 <- ggplot(plotting_frame, aes(x = g, y = KMS, color = sim_id, 
                                 linetype = selection,
                                 group = interaction(sim_id, selection, iteration))) +
  geom_line(linewidth = 0.7) +
  scale_linetype_manual(values = c("yes" = "solid", "no" = "dashed")) +
  geom_hline(yintercept = BASELINE_KMS, linetype = "dashed", color = "grey") +
  labs(title = "KMS progression", x = "Generation (cycle)", y = "Score")

# CnFS progression plot
p2 <- ggplot(plotting_frame, aes(x = g, y = CnFS, color = sim_id, 
                                 linetype = selection,
                                 group = interaction(sim_id, selection, iteration))) +
  geom_line(linewidth = 0.7) +
  scale_linetype_manual(values = c("yes" = "solid", "no" = "dashed")) +
  geom_hline(yintercept = BASELINE_CNFS, linetype = "dashed", color = "grey") +
  labs(title = "CnFS progression", x = "Generation (cycle)", y = "Score")

# patchworks the CnFS and KMS plot together with some titles and other dressing
patch_p <- p1 + p2 +
  plot_layout(ncol = 2, nrow = 1, guides = "collect", axes = "collect")
ggsave("plots/figure_S2/fig_s2c.pdf", plot = patch_p, 
       width = 14, height = 7, units = "in")

# and then the representative karyo plots for fig_s2d
p1 <- cnvHeatmap(selection_sims[[highest_KMS_sim]])
p1 <- p1 + 
  labs(title = paste0("Top KMS", " - Sim ", high_KMS_number))
p2 <- cnvHeatmap(selection_sims[[highest_CnFS_sim]])
p2 <- p2 + 
  labs(title = paste0("Top CnFS", " - Sim ", high_CnFS_number))

patch_p <- p1 + p2 +
  plot_layout(ncol = 2, nrow = 1, guides = "collect", axes = "collect")
patch_p
ggsave("plots/figure_S2/fig_s2d.pdf", plot = patch_p, 
       width = 14, height = 7, units = "in")
