library(CINsim)
library(tidyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(dplyr)

# sets wd to location of the script - this is dirty, and I do not 
# recommend it for anyone else following in my footsteps.
TMP_DIR <- Sys.getenv("TMPDIR")
setwd(TMP_DIR)

plot_dir <- file.path(
  TMP_DIR,
  "plots"
)

if (!dir.exists(plot_dir)){
  dir.create(plot_dir, recursive = TRUE)
}

CORES_AVAIL <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

# plotting theme and scale colours for copy numbers
# copy number colors
copy_num_cols <- c("gray90", "darkorchid3", "springgreen2", "red3", "gold2", "navy",
                   "lemonchiffon", "dodgerblue", "chartreuse4")
names(copy_num_cols) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")

# makes a quick plot
t302_p3 <- CINsim::t302_p3
t302_p3_melt <- melt(t302_p3)
t302_p3_melt$Copy_number <- factor(rep(0:9, times = 20))
colnames(t302_p3_melt) <- c("Chromosome", "Fraction", "Copy_number")

# plots the observed CN states per chromosome in the CINsim style
p1 <- ggplot(t302_p3_melt, aes(x = Chromosome, y = Fraction, 
                      fill = factor(Copy_number, levels = rev(levels(Copy_number))))) + # have to reverse the levels
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = copy_num_cols) +  # Unique colors for each row
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE, n.dodge = 2)) +
  labs(x = "Chromosome", y = "Frequency", fill = "Copies", title = "T302 copy number frequencies") +
  coord_cartesian(ylim = c(0, 1)) +
  cinsim_theme() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15), aspect.ratio = 1)

ggplot2::ggsave(
  plot = p1,
  file = file.path(
    plot_dir,
    "t302_p3_selection_criteria.pdf"
  )
)

# add the code for makingg plot 4C
