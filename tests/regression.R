# Run from the repository root with Rscript tests/regression.R.
suppressPackageStartupMessages(library(ggplot2))
source("gg.manhattan.R")

# Sourcing loads the function without running the example workflow.
source("simulate_gwas_data.R")
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

check("positive SNP counts are exact even with empty chromosomes or remainders", function() {
  for (n_snps in c(1, 10, 11, 22, 23, 43, 45, 101)) {
    set.seed(42)
    data <- simulate_gwas_data(n_snps = n_snps, n_traits = 2)
    stopifnot(nrow(data) == 2L * n_snps,
              length(unique(data$SNP)) == n_snps,
              all(table(data$trait) == n_snps),
              all(data$CHR >= 1 & data$CHR <= 22),
              all(data$BP >= 1 & data$BP <= 250000000),
              !anyDuplicated(data[c("SNP", "trait")]))
  }
})

check("simulator rejects invalid counts with an argument-specific error", function() {
  for (arg in c("n_snps", "n_traits", "n_chromosomes")) {
    for (value in list(0, -1, 1.5, NA_real_, Inf, "10", TRUE, numeric(), c(1, 2))) {
      args <- list(n_snps = 10, n_traits = 2, n_chromosomes = 2)
      args[[arg]] <- value
      message <- tryCatch({do.call(simulate_gwas_data, args); ""},
                          error = function(e) conditionMessage(e))
      stopifnot(grepl(paste(arg, "must be a single positive whole number"),
                      message, fixed = TRUE))
    }
  }
})

check("sourcing the simulator preserves RNG state and runs no examples", function() {
  set.seed(101)
  state <- .Random.seed
  env <- new.env(parent = globalenv())
  source("simulate_gwas_data.R", local = env)
  stopifnot(identical(state, .Random.seed),
            exists("simulate_gwas_data", envir = env, inherits = FALSE),
            !exists("gwas_data", envir = env, inherits = FALSE))
})

check("threshold arguments use plotted y-coordinates", function() {
  default <- ggplot_build(gg.manhattan(regional_data, facet_by_trait = TRUE))
  custom <- ggplot_build(gg.manhattan(regional_data, facet_by_trait = TRUE,
                                     suggestiveline = -log10(1e-3),
                                     genomewideline = -log10(1e-6)))
  off <- ggplot_build(gg.manhattan(regional_data, facet_by_trait = TRUE,
                                  suggestiveline = FALSE, genomewideline = FALSE))
  intercepts <- function(built) {
    sort(unique(unlist(lapply(built$data, function(layer) layer$yintercept))))
  }
  stopifnot(isTRUE(all.equal(intercepts(default), c(5, -log10(5e-8)))),
            identical(intercepts(custom), c(3, 6)),
            length(intercepts(off)) == 0L)
})

check("gold highlights follow the selected threshold in every plot layout", function() {
  input <- regional_data
  input$P <- c(1e-3, 5e-8, 1e-10, 1e-7)
  for (multi_chr in c(FALSE, TRUE)) {
    if (multi_chr) input$CHR <- rep(c(3, 11), 2)
    for (facets in c(FALSE, TRUE)) {
      plot <- gg.manhattan(input, facet_by_trait = facets)
      gold <- Filter(function(layer) identical(layer$aes_params$colour, "#C77B16"), plot$layers)
      stopifnot(length(gold) == 1L, nrow(gold[[1]]$data) == 2L,
                all(gold[[1]]$data$logp >= -log10(5e-8)),
                nrow(ggplot_build(plot)$data[[1]]) == nrow(input))
    }
  }
  custom <- gg.manhattan(input, genomewideline = 6, highlight_colour = "purple")
  purple <- Filter(function(layer) identical(layer$aes_params$colour, "purple"), custom$layers)
  stopifnot(length(purple) == 1L, nrow(purple[[1]]$data) == 3L)
})

check("highlighting can be disabled and is absent for raw-value plots", function() {
  input <- regional_data
  input$P <- c(1e-3, 1e-8, 1e-10, 1e-7)
  for (args in list(list(highlight = FALSE), list(genomewideline = FALSE),
                    list(logp = FALSE))) {
    plot <- do.call(gg.manhattan, c(list(x = input), args))
    stopifnot(!any(vapply(plot$layers, function(layer)
      identical(layer$aes_params$colour, "#C77B16"), logical(1))))
  }
  no_hits <- gg.manhattan(regional_data)
  stopifnot(!any(vapply(no_hits$layers, function(layer)
    identical(layer$aes_params$colour, "#C77B16"), logical(1))))
})

check("Ocean palette and custom colours retain positions on nonconsecutive chromosomes", function() {
  input <- regional_data
  input$CHR <- rep(c(3, 11), 2)
  ocean <- ggplot_build(gg.manhattan(input, facet_by_trait = TRUE))
  custom <- ggplot_build(gg.manhattan(input, cols = c("black", "grey60")))
  stopifnot(setequal(ocean$data[[1]]$colour, c("#294F79", "#278D91")),
            setequal(custom$data[[1]]$colour, c("black", "grey60")),
            all(is.finite(ocean$layout$panel_params[[1]]$x$get_breaks())),
            identical(ocean$data[[1]]$x, custom$data[[1]]$x))
  raw <- ggplot_build(gg.manhattan(input, logp = FALSE))
  stopifnot(identical(sort(raw$data[[1]]$y), sort(input$P)))
  regional <- ggplot_build(gg.manhattan(regional_data))
  stopifnot(length(unique(regional$data[[1]]$colour)) == 2L)
})

check("thresholds remain inside the y-axis without dropping points or lines", function() {
  input <- regional_data
  input$P <- rep(0.5, nrow(input))
  plot <- gg.manhattan(input, facet_by_trait = TRUE)
  built <- withCallingHandlers(ggplot_build(plot), warning = function(w) stop(w))
  for (panel in built$layout$panel_params) {
    stopifnot(max(panel$y.range) >= -log10(5e-8))
  }
  stopifnot(all(is.finite(built$data[[1]]$y)),
            length(Filter(function(layer) !is.null(layer$yintercept), built$data)) == 2L)
})

if (length(failures)) stop(paste(length(failures), "regression checks failed"))
cat("All regression checks passed.\n")
