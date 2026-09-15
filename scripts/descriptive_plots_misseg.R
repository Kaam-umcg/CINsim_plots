library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# MAGIC VARIABLES
N_CHROM <-  23 # human chromosome count
REPEATS <- 10000

# sets the wd to the node directory, otherwise plot saving doesn't work
TMP_DIR <- Sys.getenv("TMPDIR")
setwd(TMP_DIR)

# creates dir to store our figures in later
plot_dir <- file.path(
  TMP_DIR,
  "plots"
)
output_dir <- file.path(
  TMP_DIR,
  "results"
)

if (!dir.exists(plot_dir)){
  dir.create(plot_dir, recursive = TRUE)
}

if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

# function to test how many chromosomes missegregated
test_misseg_count <- function(p_misseg, N_CHROM = 23){
  misseg_count <- sum(runif(N_CHROM, 0, 1) < p_misseg)
  return(misseg_count)
}

# function to test whether a chromosome missegregated
test_misseg <- function(p_misseg, N_CHROM = 23){
  misseg_count <- sum(runif(N_CHROM, 0, 1) < p_misseg)
  if (misseg_count > 0){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# Figure S4c
# generate series of pMisseg values
pMissegs <- 10^seq(-4, 0, length.out = 50)

# dataframe for later plotting
plot_frame <- data.frame(
  p_misseg = pMissegs,
  prob_1_misseg = NA
)

# gets the misseg probability for all p_misseg values
for (i in 1:nrow(plot_frame)){
  row <- plot_frame[i, ]
  # nasty 1 liner, checks what probability there is of at least 1 misseg
  row$prob_1_misseg <- sum(unlist(lapply(rep(row$p_misseg, REPEATS), test_misseg))) / REPEATS
  plot_frame[i, ] <- row
}

# makes the actual plot
p1 <- ggplot(plot_frame, aes(x = p_misseg, y = prob_1_misseg)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans = "log10", breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1)) +
  coord_cartesian(xlim = c(min(pMissegs), max(pMissegs))) +
  labs(x = "p_misseg", y = "Frequency", title ="Mitoses with at least 1 mis-segregating chromosome",
       subtitle = "(human genome of 23 chromosomes)") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(
  filename = file.path(
    plot_dir,
    "mitosis_with_misseg.pdf"
  ),
  plot = p1
)

# define the probability values for missegregation
p_misseg_vals <- c(0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.0015, 0.001)

# create the plot frame
plot_frame <- data.frame(
  p_missegs = rep(p_misseg_vals, each = REPEATS)
)

# map the misseg count for each probability value
plot_frame$n_missegs <- map_int(plot_frame$p_missegs, ~ test_misseg_count(.x))

# summarize the frequency of missegregation counts
frequency_data <- plot_frame %>%
  group_by(p_missegs, n_missegs) %>%
  summarise(frequency = n(), .groups = 'drop')

# ensure p_missegs is treated as a factor
frequency_data$p_missegs <- factor(frequency_data$p_missegs, levels = unique(frequency_data$p_missegs))

# plot the data
p2 <- ggplot(frequency_data, aes(x = n_missegs, y = frequency / REPEATS, color = p_missegs, group = p_missegs)) +
  geom_line() +
  labs(
    title = "Number of mis-segregation events per mitosis",
    subtitle = "(human genome)",
    x = "Number of mis-segregation events",
    y = "Probability",
    color = "p_missegs") +
    scale_x_continuous(breaks = 0:8, limits = c(0, 9)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(
  filename = file.path(
    plot_dir,
    "missegs_per_mitosis.pdf"
  ),
  plot = p2
)