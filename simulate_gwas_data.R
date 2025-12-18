# Install ggplot2 if not already installed
if (!require(ggplot2)) install.packages("ggplot2")

# Set seed for reproducibility
set.seed(42)

source("gg.manhattan.R")

# Function to simulate GWAS-like data
simulate_gwas_data <- function(n_snps = 40000, n_traits = 3, n_chromosomes = 22) {
  
  data <- data.frame()
  
  # Create SNPs distributed across chromosomes
  for (chr in 1:n_chromosomes) {
    # Number of SNPs per chromosome (roughly proportional to chromosome size)
    snps_per_chr <- round(n_snps / n_chromosomes)
    
    for (trait in 1:n_traits) {
      # Generate base pair positions (realistic range per chromosome)
      bp_positions <- sort(sample(1:250000000, snps_per_chr, replace = FALSE))
      
      # Generate p-values with some realistic structure
      # Most SNPs have non-significant p-values
      # A few have strong signals (smaller p-values)
      p_values <- runif(snps_per_chr, min = 0.001, max = 1)
      
      # Add some "hits" (significant SNPs) - typically 2-5 per chromosome per trait with more SNPs
      n_hits <- sample(1:5, 1)
      if (n_hits > 0) {
        hit_indices <- sample(1:snps_per_chr, n_hits, replace = FALSE)
        p_values[hit_indices] <- runif(n_hits, min = 1e-8, max = 1e-4)
      }
      
      # Create SNP identifiers
      snp_ids <- paste0("rs", chr, "_", seq_along(bp_positions))
      
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

# Generate example data for 3 traits with 40k SNPs
gwas_data <- simulate_gwas_data(n_snps = 40000, n_traits = 3, n_chromosomes = 22)

# View the structure
head(gwas_data, 10)
str(gwas_data)

# Summary statistics
print(paste("Total SNPs:", nrow(gwas_data)))
print(paste("Number of traits:", length(unique(gwas_data$trait))))
print(paste("Number of chromosomes:", length(unique(gwas_data$CHR))))
print(paste("P-value range:", min(gwas_data$P), "to", max(gwas_data$P)))

# Example: Plot a single trait
single_trait <- subset(gwas_data, trait == "Trait_1")
gg.manhattan(single_trait, plot_title = "Manhattan Plot - Trait 1")

# Example: Plot all traits together with faceting
gg.manhattan(gwas_data, plot_title = "Manhattan Plot - Multiple Traits", facet_by_trait = TRUE)
