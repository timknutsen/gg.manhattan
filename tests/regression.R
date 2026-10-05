# Run from the repository root with Rscript tests/regression.R.
suppressPackageStartupMessages(library(ggplot2))
source("gg.manhattan.R")

# Load the simulator function without running the interactive examples.
simulator_code <- parse("simulate_gwas_data.R")
for (expr in simulator_code) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("simulate_gwas_data"))) {
    eval(expr)
  }
}
stopifnot(exists("simulate_gwas_data", mode = "function"))

failures <- character()
check <- function(name, test) {
  tryCatch({
    test()
    cat("PASS:", name, "\n")
  }, error = function(e) {
    failures <<- c(failures, name)
    cat("FAIL:", name, "-", conditionMessage(e), "\n")
  })
}

regional_data <- data.frame(
  CHR = rep(1, 4), BP = rep(c(10e6, 20e6), 2),
  P = c(1e-3, 1e-6, 1e-4, 1e-7),
  SNP = rep(c("rs1", "rs2"), 2), trait = rep(c("A", "B"), each = 2)
)

check("single-chromosome facets retain readable Mb positions", function() {
  plot <- gg.manhattan(regional_data, facet_by_trait = TRUE)
  built <- ggplot_build(plot)
  stopifnot(nrow(built$layout$layout) == 2L,
            identical(sort(unique(built$data[[1]]$x)), c(10, 20)))
  for (panel in built$layout$panel_params) {
    ticks <- panel$x$get_breaks()
    ticks <- ticks[is.finite(ticks)]
    stopifnot(length(ticks) >= 2L, any(ticks >= 10 & ticks <= 20))
  }
  stopifnot(identical(plot$scales$get_scales("x")$name,
                      "Chromosome 1 position(Mb)"))
})

check("multi-chromosome facets retain custom chromosome labels", function() {
  input <- regional_data
  input$CHR <- rep(c(1, 2), 2)
  plot <- gg.manhattan(input, facet_by_trait = TRUE, chrlabs = c("Chr1", "Chr2"))
  built <- ggplot_build(plot)
  stopifnot(nrow(built$layout$layout) == 2L,
            identical(plot$scales$get_scales("x")$labels, c("Chr1", "Chr2")),
            all(is.finite(built$layout$panel_params[[1]]$x$get_breaks())))
})

check("SNP coordinates are shared across traits with independent p-values", function() {
  set.seed(42)
  data <- simulate_gwas_data(n_snps = 100, n_traits = 3, n_chromosomes = 2)
  stopifnot(nrow(data) == 300L)
  coords <- split(data$BP, data$SNP)
  chromosomes <- split(data$CHR, data$SNP)
  stopifnot(all(vapply(coords, function(x) length(unique(x)) == 1L, logical(1))),
            all(vapply(chromosomes, function(x) length(unique(x)) == 1L, logical(1))))
  a <- subset(data, trait == "Trait_1")
  b <- subset(data, trait == "Trait_2")
  b <- b[match(a$SNP, b$SNP), ]
  stopifnot(identical(a$BP, b$BP), any(a$P != b$P),
            all(is.finite(data$P) & data$P > 0 & data$P <= 1))
  set.seed(42)
  stopifnot(identical(data, simulate_gwas_data(100, 3, 2)))
})

check("small datasets support one through five SNPs per chromosome", function() {
  for (snps_per_chr in 1:5) {
    for (seed in c(1, 42, 99)) {
      set.seed(seed)
      data <- simulate_gwas_data(n_snps = 22 * snps_per_chr,
                                 n_traits = 3, n_chromosomes = 22)
      stopifnot(nrow(data) == 22L * snps_per_chr * 3L,
                all(table(data$CHR, data$trait) == snps_per_chr),
                all(is.finite(data$P) & data$P > 0 & data$P <= 1))
    }
  }
})

check("default-size simulated data builds faceted and single-trait plots", function() {
  set.seed(42)
  data <- simulate_gwas_data()
  faceted <- ggplot_build(gg.manhattan(data, facet_by_trait = TRUE))
  single <- ggplot_build(gg.manhattan(subset(data, trait == "Trait_1")))
  stopifnot(nrow(faceted$layout$layout) == 3L,
            nrow(faceted$data[[1]]) == nrow(data),
            nrow(single$data[[1]]) == sum(data$trait == "Trait_1"),
            all(is.finite(faceted$data[[1]]$x)), all(is.finite(faceted$data[[1]]$y)))
})

if (length(failures)) stop(paste(length(failures), "regression checks failed"))
cat("All regression checks passed.\n")
