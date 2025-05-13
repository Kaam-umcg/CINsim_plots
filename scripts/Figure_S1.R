library(CINsim)
library(tidyverse)

setwd(Sys.getenv("TMPDIR"))

if (!dir.exists("plots/figure_S1")){
  dir.create("plots/figure_S1", recursive = TRUE)
}

# sets seed for reproducibility 
set.seed(42)

# plotting theme and scale colours for copy numbers
# copy number colors
copy_num_cols <- c("gray90", "darkorchid3", "springgreen2", "red3", "gold2", "navy",
                   "lemonchiffon", "dodgerblue", "chartreuse4", "lightcoral", "aquamarine2")
names(copy_num_cols) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

# plotting theme
cinsim_theme <- function() {
  theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          aspect.ratio = 1)
}

# Figure S1a
# mis-segregation probability per copy number
pMissegs <- 10^seq(-6, 0, length.out = 50)
copy_nums <- 1:8
probs <- map_df(pMissegs, function(x) {
  
  p_weighted <- 1 - (1 - x)^(copy_nums*2)
  tibble(copy_num = copy_nums,
         p_misseg = x,
         prob = p_weighted)
  
})

probs %>%
  ggplot(aes(x = p_misseg, y = prob)) +
  geom_line(aes(group = copy_num, col = factor(copy_num, levels = 8:1)), linewidth = 1) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                breaks = 10^(-6:0)) +
  scale_color_manual(values = copy_num_cols) +
  labs(x = "pMisseg", y = "Probability", col = "Copy",
       title = "Mis-segregation probability per copy number") +
  cinsim_theme()
ggsave("plots/figure_S1/figure_s1a.pdf")

# Figure S1b
# generate series of downsampling fractions and population sizes
down_sample_frac <- 2^(-1:-8)
down_sample <- as.vector(c(1, 5) %o% 10^(3:5))

# make an arguments data frame, making all possible combinations
args <- expand.grid(down_sample_frac = down_sample_frac, 
                    down_sample = down_sample)

# generate 200 euploid karyotypes
karyotypes <- makeKaryotypes(numCell = 200, species = "mouse", copies = 2)

# simulate 50 rounds of cell divisions, without CIN
sim_data <- pmap(args, parallelCinsim_simple, 
                 karyotypes = karyotypes,
                 g = 50, 
                 pMisseg = 0,
                 max_num_cells = NULL,
                 cores = 6)

# extract the number of surviving clones per iteration per run
sim_data_compiled <- map_df(sim_data, function(x) {
  
  # get relevant info from the sim_info element for each karyoSim object
  sim_info <- map(x, "sim_info")[[1]]
  ds <- sim_info[["down_sample"]]
  dsf <- sim_info[["down_sample_frac"]]
  
  # extract the number of surviving clones
  num_clones <- x %>%
    map("karyotypes") %>%
    map_dbl(function(y) {
      y %>%
        row.names() %>%
        unique() %>%
        length()
    })
  
  # compile results into a tibble for plotting and return
  result <- tibble(down_sample = as.numeric(ds),
                   down_sample_frac = as.numeric(dsf),
                   group = paste(down_sample, down_sample_frac, sep = "_"),
                   num_clones = num_clones)
  return(result)
  
})

# plot mean clonal frequency as heatmap
sim_data_compiled %>%
  mutate(freq = num_clones/nrow(karyotypes),
         down_sample_frac = factor(down_sample_frac^-1),
         down_sample = factor(down_sample)) %>%
  group_by(group, down_sample, down_sample_frac) %>%
  summarize(mean_freq = mean(freq, na.rm = TRUE),
            sd_freq = sd(freq, na.rm = TRUE)) %>%
  ggplot(aes(x = down_sample, y = down_sample_frac)) +
  geom_tile(aes(fill = mean_freq), col = "grey") +
  geom_text(aes(label = signif(mean_freq, 2)),
            col = "white") +
  geom_rect(aes(xmin = 3.5, xmax = 4.5, ymin = 1.5, ymax = 2.5),
            color ="red", fill = NA, linewidth = 1) +
  scale_fill_viridis_c(begin = 0, end = 1) +
  scale_y_discrete(labels = scales::math_format(frac(1, .x))) +
  
  labs(x = "Down sample threshold", y = "Down sample fraction", 
       fill = "Fraction of clones",
       title = "Clonal survival to 50 generation (200 clones)") 
ggsave("plots/figure_S1/fig_s1b.pdf", width = 7, height = 7, unit = "in")

# Figure S1c
pMissegs <- 10^seq(-6, 0, length.out = 30)
iterations <- 1:6

sim_data <- vector(mode = "list", length = max(iterations))
# simulate 6 iterations, race to 10 billion cells
for(i in iterations) {
  sim_data[[i]] <- parallelCinsim_simple(g = 200,
                                         pMisseg = pMissegs,
                                         max_num_cells = 10e+9,
                                         cores = 6)
}

# extract the number of mis-segregation events per simulation
sim_data_compiled <- map2_df(sim_data, 1:length(sim_data), function(x, y) {
  
  sim_info <- tibble(sim_label = paste0("sim_", 1:length(x)),
                     p_misseg = x %>% map_dbl(function(x) {
                       p_misseg <- x$sim_info[["pMisseg"]]
                       return(as.numeric(p_misseg))
                     }))
  
  result <- x %>%
    map("gen_measures") %>%
    map2_df(names(x), function(z, w) {
      z %>%
        mutate(sim_label = w)}) %>%
    inner_join(sim_info, by = "sim_label") %>%
    select(-sim_label) %>%
    mutate(replicate = y)
  return(result)
  
})

# plot the net growth rate (i.e. the mean number of population doublings)
sim_data_compiled %>%
  mutate(sim_label = paste0(p_misseg, replicate, sep = "_")) %>%
  arrange(replicate, p_misseg, g) %>%
  group_by(sim_label) %>%
  mutate(true_cell_count_prev = lag(true_cell_count, n = 1, default = 1),
         growth = true_cell_count / true_cell_count_prev) %>%
  filter(true_cell_count >= 20) %>%
  ungroup() %>%
  group_by(p_misseg) %>%
  summarize(mean_g = mean(growth, na.rm = TRUE),
            sd_g = sd(growth, na.rm = TRUE)) %>%
  ggplot(aes(x = p_misseg, y = mean_g)) +
  geom_hline(yintercept = 1, linetype = 2, col = "grey") +
  geom_linerange(aes(ymin = mean_g - sd_g,
                     ymax = mean_g + sd_g)) +
  geom_point() +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                breaks = 10^(-6:0), limits = c(10^-6, 1)) +
  expand_limits(y = 1) +
  labs(x = "pMisseg", y = "Growth",
       title = "Net growth") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"),
        aspect.ratio = 1)
ggsave(file = "plots/figure_S1/figure_s1c.pdf")

# Figure S1d
# get the mean + sd of the number of generation required to get to 10 billions cells
plotting_data <- sim_data_compiled %>%
  filter(true_cell_count > 10e+9) %>%
  group_by(p_misseg) %>%
  summarize(mean_g = mean(g, na.rm = TRUE),
            sd_g = sd(g, na.rm = TRUE)) %>%
  ungroup()

# plot the results
plotting_data %>%
  ggplot(aes(x = p_misseg, y = mean_g)) +
  geom_hline(yintercept = min(plotting_data$mean_g), col = "lightgrey", linetype = 2) +
  geom_linerange(aes(ymin = mean_g - sd_g,
                     ymax = mean_g + sd_g)) +
  geom_point(size = 2) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                breaks = 10^(-6:0), limits = c(10^-6, 1)) +
  scale_y_continuous(breaks = seq(0, 80, 10)) +
  expand_limits(y = 0) +
  labs(x = "pMisseg", y = "G") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"),
        aspect.ratio = 1)
ggsave(file = "plots/figure_S1/figure_s1d.pdf")

# Figure S1e is handmade, so we skip it here

# Figure S1f 
# reconstruct a most-fit karyotype based on Mps1 data
best_karyo <- as.integer(apply(Mps1, 2, function(x) rownames(Mps1)[which.max(x)]))
best_karyo <- matrix(best_karyo, nrow = 1)
colnames(best_karyo) <- 1:20
colnames(best_karyo)[20] <- "X"
rownames(best_karyo) <- "cell_1"

# generate a million random karyotypes, append the optimal and euploid karyotypes
num_karyotypes <- 1e+06
karyotypes <- rbind(makeKaryotypes(num_karyotypes, random = TRUE),
                    matrix(sample(1:8, 20*num_karyotypes, replace = TRUE), 
                           ncol = 20, nrow = num_karyotypes),
                    makeKaryotypes(1),
                    best_karyo)
rownames(karyotypes) <- paste(c(rep("Aneuploid-semirandom", num_karyotypes),
                                rep("Aneuploid-random", num_karyotypes),
                                "Euploid", "Optimal"), 
                              1:(num_karyotypes + 2), sep = "_")

# calculate the karyotype scores
# triple colon (:::) calls the internal variable in the package - the manual says to contact the maintainer
# if you're doing this as it's probably a mistake. I, as the maintainer, say the maintainer is stupid so go ahead
scores <- bind_rows(tibble(method = "cn_based",
                           cell = rownames(karyotypes),
                           score = CINsim:::get_score(karyotypes,
                                             selection_mode = "cn_based",
                                             selection_metric = Mps1_X)),
                    tibble(method = "rel_copy",
                           cell = rownames(karyotypes),
                           score = CINsim:::get_score(karyotypes,
                                             selection_mode = "rel_copy",
                                             selection_metric = Mps1_rel_copy))) %>%
  separate(cell, into = c("group", "id"), sep = "_", remove = TRUE) %>%
  mutate(group = factor(group, levels = c("Euploid", "Aneuploid-random", "Aneuploid-semirandom", "Optimal")))

# get the euploid and optimal scores
scores_opt <- scores %>%
  filter(group %in% c("Euploid", "Optimal"),
         method == "cn_based")

# plot the distribution of scores, highlighting the euploid and maximum score
p1 <- scores %>%
  filter(group %in% c("Aneuploid-random", "Aneuploid-semirandom"),
         method == "cn_based") %>%
  ggplot(aes(x = score)) +
  geom_density(aes(group = group, fill = group), col = "grey", alpha = 0.5) +
  geom_hline(yintercept = 0, col = "black") +
  geom_vline(data = scores_opt, aes(group = group, col = group, xintercept = score),
             linewidth = 1) +
  scale_color_manual(values = c("darkgreen", "darkred")) +
  labs(x = "Score", y = "Density", fill = "Karyotypes",
       col = "Reference") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"),
        aspect.ratio = 1)

# score to probability through linear fit
score2prob_lm <- function(x, a, b) {(a*x + b)}

# generate a series of values between the euploid probability and the maximum
sur_0.5 <- seq(0.5, 1, length.out = 1000)
sur_0.9 <- seq(0.9, 1, length.out = 1000)
scores2fit <- seq(min(scores_opt$score), max(scores_opt$score), length.out = 1000)

# make a fit for 0.5 and 0.9 euploid survival
score_tibble <- tibble(scores = scores2fit,
                       sur_0.5 = sur_0.5,
                       sur_0.9 = sur_0.9)

# make the fit and get coefficients
fit_0.5 <- lm(sur_0.5 ~ scores2fit)
fit_0.9 <- lm(sur_0.9 ~ scores2fit)
coef_tibble <- rbind(coef(fit_0.5), coef(fit_0.9)) %>% as.data.frame()
colnames(coef_tibble) <- c("b", "a")
coef_tibble$euploid_sur <- c(0.5, 0.9)

# add the psurvival lines to the plot, making a secondary axis
# divide coeffcients by three to fit in in the frame
p1 <- p1 + 
  geom_abline(slope = coef_tibble$a[1]/3, intercept = coef_tibble$b[1]/3, col = "orange", size = 1) +
  geom_abline(slope = coef_tibble$a[2]/3, intercept = coef_tibble$b[2]/3, col = "red", size = 1) +
  scale_y_continuous(sec.axis = sec_axis(transform = ~ .*3, name = "pSurvival", breaks = seq(0, 1, 0.25)))
ggsave(file = "plots/figure_S1/figure_s1d.pdf", plot = p1)  

# Figure S1g is handmade, so we skip it here
