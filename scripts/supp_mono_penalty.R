# clean the workspace
rm(list = ls())

setwd(Sys.getenv("TMPDIR"))

`%>%` <- magrittr::`%>%`

plot_dir <- file.path(
    Sys.getenv("TMPDIR"),
    "plots"
)
output_dir <- file.path(
  Sys.getenv("TMPDIR"),
  "output"
)

if (!dir.exists(plot_dir)){
  dir.create(plot_dir, recursive = TRUE)
}

if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

# makes a visualisation of the monosomy penalty count
p_modal <- 0.7
p_mono <- 0.2

# this plot needs:
#' 1. colours for the lines
#' 2. Clearer annotation on the x & y -axis, maybe use custom ticks?
#' 3. Better title
#' 4. General CINsim dressing? use shared theme
p1 <- ggplot2::ggplot() +
    ggplot2::geom_abline(
        intercept = 0,
        slope = 1 / p_modal,
        colour = "black",
        linetype = "solid"
    ) +
    ggplot2::geom_vline(
        xintercept = p_modal - p_mono,
        colour = "red",
        linetype = "dashed") +
    ggplot2::geom_hline(
        yintercept = (p_modal - p_mono) / p_modal,
        colour = "red",
        linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
        x = "\U0394 p_modal - p_mono",
        y = "n_mono",
        title = "Monosomy penalty is defined by \U0394 (p_modal - p_mono)"
    ) +
    ggplot2::lims(
        x = c(0, p_modal),
        y = c(0, 1)
    )

ggplot2::ggsave("~/mono_penalty_plot.png")