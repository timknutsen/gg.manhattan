# Simulate GWAS-like summary statistics with regional association peaks.
# Set the seed before calling for reproducible results. Loci are shared across
# traits; local noise and signal strengths vary between traits.
simulate_gwas_data <- function(n_snps = 40000, n_traits = 3, n_chromosomes = 22,
                               n_peaks = 3, peak_width = 3e6, peak_effect = 7) {
  check_count <- function(value, name, minimum = 1) {
    if (!(length(value) == 1L && is.numeric(value) && !is.complex(value) &&
          is.finite(value) && value >= minimum && value == floor(value) &&
          value <= .Machine$integer.max)) {
      kind <- if (minimum == 0) "nonnegative" else "positive"
      stop(paste(name, "must be a single", kind, "whole number no larger than",
                 .Machine$integer.max), call. = FALSE)
    }
  }
  check_count(n_snps, "n_snps")
  check_count(n_traits, "n_traits")
  check_count(n_chromosomes, "n_chromosomes")
  check_count(n_peaks, "n_peaks", minimum = 0)
  for (name in c("peak_width", "peak_effect")) {
    value <- get(name)
    if (!(length(value) == 1L && is.numeric(value) && !is.complex(value) &&
          is.finite(value) && value > 0)) {
      stop(paste(name, "must be a single finite positive number"), call. = FALSE)
    }
  }
  if (ceiling(n_snps / n_chromosomes) > 250000000) {
    stop("Cannot place more than 250000000 distinct SNPs on a chromosome",
         call. = FALSE)
  }

  # Construct shared loci once and preserve the exact requested total.
  per_chr <- n_snps %/% n_chromosomes
  remainder <- n_snps %% n_chromosomes
  chromosomes <- seq_len(min(n_chromosomes, n_snps))
  loci <- do.call(rbind, lapply(chromosomes, function(chr) {
    n <- per_chr + as.integer(chr <= remainder)
    data.frame(
      CHR = chr,
      BP = sort(sample.int(250000000, n, replace = FALSE)),
      SNP = paste0("rs", chr, "_", seq_len(n))
    )
  }))
  rownames(loci) <- NULL
  chr_rows <- split(seq_len(n_snps), loci$CHR)
  peak_n <- min(n_peaks, n_snps)
  shared_centre <- if (peak_n > 0) sample.int(n_snps, 1L) else integer()
  results <- vector("list", n_traits)
  peak_info <- vector("list", n_traits)

  for (trait in seq_len(n_traits)) {
    trait_name <- paste0("Trait_", trait)
    centres <- shared_centre
    while (length(centres) < peak_n) {
      remaining <- setdiff(seq_len(n_snps), centres)
      # Prefer different chromosomes so separate peaks remain distinguishable.
      preferred <- remaining[!loci$CHR[remaining] %in% loci$CHR[centres]]
      if (!length(preferred)) {
        separated <- vapply(remaining, function(idx) {
          all(loci$CHR[centres] != loci$CHR[idx] |
                abs(loci$BP[centres] - loci$BP[idx]) >= 4 * peak_width)
        }, logical(1))
        preferred <- remaining[separated]
      }
      candidates <- if (length(preferred)) preferred else remaining
      centres <- c(centres, candidates[sample.int(length(candidates), 1L)])
    }

    # Distance-correlated standard-normal noise gives neighboring markers
    # related statistics while retaining a uniform marginal null p-value.
    z <- numeric(n_snps)
    for (idx in chr_rows) {
      noise <- rnorm(length(idx))
      if (length(idx) > 1L) {
        distance <- diff(loci$BP[idx])
        rho <- exp(-distance / (peak_width / 6))
        innovation_sd <- sqrt(-expm1(-2 * distance / (peak_width / 6)))
        for (j in seq.int(2L, length(idx))) {
          noise[j] <- rho[j - 1L] * noise[j - 1L] + innovation_sd[j - 1L] * noise[j]
        }
      }
      z[idx] <- noise
    }

    effects <- peak_effect * runif(peak_n, 0.9, 1.1) *
      sample(c(-1, 1), peak_n, replace = TRUE)
    for (k in seq_len(peak_n)) {
      centre <- centres[k]
      idx <- chr_rows[[as.character(loci$CHR[centre])]]
      distance <- (loci$BP[idx] - loci$BP[centre]) / peak_width
      z[idx] <- z[idx] + effects[k] * exp(-0.5 * distance^2)
    }
    p_values <- pmax(.Machine$double.xmin, 2 * pnorm(-abs(z)))
    results[[trait]] <- data.frame(
      CHR = loci$CHR, BP = loci$BP, P = p_values,
      SNP = loci$SNP, trait = trait_name
    )
    peak_info[[trait]] <- data.frame(
      trait = rep(trait_name, peak_n),
      CHR = loci$CHR[centres], BP = loci$BP[centres], SNP = loci$SNP[centres],
      effect = effects, width = rep(peak_width, peak_n)
    )
  }

  data <- do.call(rbind, results)
  rownames(data) <- NULL
  attr(data, "peaks") <- do.call(rbind, peak_info)
  data
}

# Run examples only with Rscript, so source() just loads the function.
if (sys.nframe() == 0L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 before running the examples: install.packages('ggplot2')",
         call. = FALSE)
  }
  source("gg.manhattan.R")
  set.seed(42)
  gwas_data <- simulate_gwas_data()

  print(head(gwas_data, 10))
  str(gwas_data)
  cat("Distinct SNPs:", length(unique(gwas_data$SNP)), "\n")
  cat("Rows across traits:", nrow(gwas_data), "\n")
  cat("Number of traits:", length(unique(gwas_data$trait)), "\n")
  cat("Number of chromosomes:", length(unique(gwas_data$CHR)), "\n")
  cat("P-value range:", range(gwas_data$P), "\n")

  single_trait <- subset(gwas_data, trait == "Trait_1")
  print(gg.manhattan(single_trait, plot_title = "Manhattan Plot - Trait 1"))
  print(gg.manhattan(gwas_data, plot_title = "Manhattan Plot - Multiple Traits",
                    facet_by_trait = TRUE))
}
