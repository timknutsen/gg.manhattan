# Run from the repository root: Rscript examples/ocean.R
suppressPackageStartupMessages(library(ggplot2))
source("gg.manhattan.R")
source("simulate_gwas_data.R")

set.seed(20261021)
example_data <- simulate_gwas_data(
  n_snps = 40000, n_traits = 3, n_chromosomes = 22,
  n_peaks = 3, peak_width = 3e6, peak_effect = 7
)
example_data$trait <- sub("_", " ", example_data$trait)
ocean_plot <- gg.manhattan(example_data, facet_by_trait = TRUE,
                          plot_title = "Manhattan plots by trait") +
  labs(subtitle = "Regional association peaks with gold significance highlights",
       caption = paste0("SIMULATED DATA  |  40,000 loci per trait; correlated neighboring statistics\n",
                        "Dotted: p = 10^-5    Dashed: p = 5 x 10^-8    Gold: p <= 5 x 10^-8"))
ggsave("Rplot.png", ocean_plot, width = 12, height = 8, dpi = 160, bg = "white")
