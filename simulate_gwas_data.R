# Simulate exactly n_snps distinct loci, repeated once per trait.
# Counts must be positive whole numbers. Set the seed before calling for
# reproducible results. Positions are shared across traits; p-values vary.
simulate_gwas_data <- function(n_snps = 40000, n_traits = 3, n_chromosomes = 22) {
  check_count <- function(value, name) {
    if (!(length(value) == 1L && is.numeric(value) && !is.complex(value) &&
          is.finite(value) && value >= 1 && value == floor(value) &&
          value <= .Machine$integer.max)) {
      stop(paste(name, "must be a single positive whole number no larger than",
                 .Machine$integer.max), call. = FALSE)
    }
  }
  check_count(n_snps, "n_snps")
  check_count(n_traits, "n_traits")
  check_count(n_chromosomes, "n_chromosomes")
  
  data <- data.frame()
  snps_per_chr <- n_snps %/% n_chromosomes
  remainder <- n_snps %% n_chromosomes
  
  # Distribute the remainder over the first chromosomes to preserve the total.
  for (chr in seq_len(n_chromosomes)) {
    n_chr_snps <- snps_per_chr + as.integer(chr <= remainder)
    if (n_chr_snps == 0L) next
    # Each SNP has the same genomic position across traits.
    bp_positions <- sort(sample.int(250000000, n_chr_snps, replace = FALSE))
    snp_ids <- paste0("rs", chr, "_", seq_along(bp_positions))

    for (trait in seq_len(n_traits)) {
      # Generate p-values with some realistic structure
      # Most SNPs have non-significant p-values
      # A few have strong signals (smaller p-values)
      p_values <- runif(n_chr_snps, min = 0.001, max = 1)
      
      # Add 1-5 hits per chromosome and trait, bounded by the SNP count.
      n_hits <- sample.int(min(5L, n_chr_snps), 1L)
      if (n_hits > 0) {
        hit_indices <- sample.int(n_chr_snps, n_hits, replace = FALSE)
        p_values[hit_indices] <- runif(n_hits, min = 1e-8, max = 1e-4)
      }
      
      # Combine into dataframe
      trait_data <- data.frame(
        CHR = chr,
        BP = bp_positions,
        P = p_values,
        SNP = snp_ids,
        trait = paste0("Trait_", trait)
      )
      
      data <- rbind(data, trait_data)
    }
  }
  
  return(data)
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
