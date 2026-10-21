# Run from the repository root: Rscript examples/ocean.R
suppressPackageStartupMessages(library(ggplot2))
source("gg.manhattan.R")
source("simulate_gwas_data.R")

set.seed(20261021)
example_data <- simulate_gwas_data(18000, 3, 22)
example_data$P <- runif(nrow(example_data), min = 1e-5, max = 1)

# Synthetic peaks for a visual example; these are not observed GWAS results.
peaks <- list(
  Trait_1 = data.frame(chr = c(3, 11, 18), centre = c(95, 155, 80), height = c(12, 9, 7)),
  Trait_2 = data.frame(chr = c(6, 11, 20), centre = c(130, 155, 175), height = c(10, 11, 8)),
  Trait_3 = data.frame(chr = c(3, 15, 20), centre = c(95, 100, 175), height = c(8, 12, 9))
)
for (trait in names(peaks)) {
  for (i in seq_len(nrow(peaks[[trait]]))) {
    peak <- peaks[[trait]][i, ]
    idx <- which(example_data$trait == trait & example_data$CHR == peak$chr)
    bump <- peak$height * exp(-0.5 * ((example_data$BP[idx] / 1e6 - peak$centre) / 4.5)^2)
    example_data$P[idx] <- 10^(-(-log10(example_data$P[idx]) + bump))
  }
}
example_data$trait <- sub("_", " ", example_data$trait)
ocean_plot <- gg.manhattan(example_data, facet_by_trait = TRUE,
                          plot_title = "Manhattan plots by trait") +
  labs(subtitle = "Ocean palette with gold significance highlights",
       caption = paste0("SIMULATED DATA  |  18,000 loci per trait\n",
                        "Dotted: p = 10^-5    Dashed: p = 5 x 10^-8    Gold: p <= 5 x 10^-8"))
ggsave("Rplot.png", ocean_plot, width = 12, height = 8, dpi = 160, bg = "white")
