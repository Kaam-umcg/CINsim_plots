library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# MAGIC VARIABLES
N_CHROM <-  23 # human chromosome count
REPEATS <- 10000

# sets the wd to the node directory, otherwise plot saving doesn't work
setwd(Sys.getenv("TMPDIR"))

if (!dir.exists("plots/figure_S4c")){
  dir.create("plots/figure_S4c", recursive = TRUE)
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
pMissegs <- 10^seq(-6, 0, length.out = 50)

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
ggplot(plot_frame, aes(x = p_misseg, y = prob_1_misseg)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans = "log10", breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)) +
  coord_cartesian(xlim = c(min(pMissegs), max(pMissegs))) +
  labs(x = "p_misseg", y = "Frequency", title ="Mitoses with at least 1 mis-segregating chromosome",
       subtitle = "(human genome of 23 chromosomes)") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
ggsave("plots/figure_S4c/figure_s4c.pdf")

# Figure S4d
# define the probability values for missegregation
p_misseg_vals <- c(0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.0015, 0.001)

# create the plot frame
plot_frame <- data.frame(
  p_missegs = rep(p_misseg_vals, each = REPEATS)
)

# map the misseg count for each probability value
plot_frame$n_missegs <- map_int(plot_frame$p_missegs, ~test_misseg_count(.x))

# summarize the frequency of missegregation counts
frequency_data <- plot_frame %>%
  group_by(p_missegs, n_missegs) %>%
  summarise(frequency = n(), .groups = 'drop')

# ensure p_missegs is treated as a factor
frequency_data$p_missegs <- factor(frequency_data$p_missegs, levels = unique(frequency_data$p_missegs))

# plot the data
ggplot(frequency_data, aes(x = n_missegs, y = frequency / REPEATS, color = p_missegs, group = p_missegs)) +
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
ggsave("plots/figure_S4c/figure_s4d.pdf")