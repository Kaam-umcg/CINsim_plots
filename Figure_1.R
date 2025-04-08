library(CINsim)
library(tidyverse)

# sets wd to location of the script - this is dirty, and I do not 
# recommend it for anyone else following in my footsteps.
# Only works in Rstudio
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
0
if (!dir.exists("../plots/figure_1")){
  dir.create("../plots/figure_1")
}

# plotting theme
cinsim_theme <- function() {
  theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          aspect.ratio = 1)
}

# set the seed for the file, as we do simulations.
# this makes future reruns reproducible (I think)!
set.seed(42)

# Figure 1a is created manually using assets created from CINsim
# Figure 1b is manually created
# Figure 1c is manually created

# Figure 1d can be reproduced using the code below:
# define logistic growth function
logistic_growth <- function(pt0, r, M) {
  pt1 <- pt0 + r*(pt0*(1 - (pt0/M)))
}

# set logistic growth parameters
# r = growth rate per generation (t)
r <- seq(0.5, 1, by = 0.05)
# carrying capacity
M <- 10e+09
# max generations
g <- 100

# get data for given parameters
growth_df <- map_df(r, function(x) {
  logistic <- vector(length = g+1)
  logistic[1] <- 1
  exponential_growth <- vector(length = g+1)
  exponential_growth[1] <- 1
  for(i in 1:g) {
    logistic[i+1] <- logistic_growth(logistic[i], x, M)
    exponential_growth[i+1] <- exponential_growth[i] + (exponential_growth[i]*x)
  }
  data <- tibble(r = x,
                 g = 0:g,
                 Logistic = logistic,
                 Exponential = exponential_growth)
  return(data)
}) %>%
  gather(Logistic, Exponential, key = "method", value = "size") %>%
  mutate(r = factor(r))

# plot the population size over time
growth_df %>%
  ggplot(aes(x = g, y = size)) +
  geom_hline(yintercept = c(M/5, M), col = "grey", linetype = 2) +
  geom_line(aes(group = r, col = r), linewidth = 1) +
  scale_color_viridis_d() +
  coord_cartesian(ylim = c(0, M),
                  xlim = c(0, g*0.7)) +
  facet_wrap(~method) +
  labs(x = "Generation", y = "Population size") +
  cinsim_theme()
ggsave(file = "plots/figure_1/figure_1d.pdf")

# Figure 1e can be reproduced using the code below:
# plot the maximum number of generations before 99.9% of M is reached
growth_df %>%
  filter(size <= 0.999*M) %>%
  group_by(r, method) %>%
  summarize(g = max(g),
            size = max(size)) %>%
  ungroup() %>%
  mutate(r = factor(r, levels = rev(unique(r)))) %>%
  ggplot(aes(x = r, y = g)) +
  geom_bar(aes(group = method, fill = method),
           col = "grey", stat = "identity", position = "dodge") +
  geom_text(aes(group = method, label = g), 
            position=position_dodge(width=0.9), vjust=-0.25) +
  labs(x = "R", y = "Generations", fill = "Method",
       title = "Generations to 99.9% of M") +
  cinsim_theme()
ggsave(file = "plots/figure_1/figure_1e.pdf")
