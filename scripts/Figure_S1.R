library(CINsim)
library(tidyverse)

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
ggsave()
