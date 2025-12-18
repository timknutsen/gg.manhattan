# Manhattan Plot Visualization

A ggplot2-based function for creating publication-quality Manhattan plots from GWAS (Genome-Wide Association Study) data.

## Example Output

![Manhattan Plot Example](Rplot.png)

## Files

- **gg.manhattan.R** - Main plotting function
- **simulate_gwas_data.R** - Script to generate example GWAS data

## Features

- Single or multi-chromosome visualization
- Support for multiple traits with faceted plots
- Customizable significance threshold lines
- Automatic chromosome labeling and positioning
- Clean, publication-ready styling

## Quick Start

```r
source("gg.manhattan.R")
source("simulate_gwas_data.R")

# Generate example data (40k SNPs, 3 traits)
gwas_data <- simulate_gwas_data(n_snps = 40000, n_traits = 3, n_chromosomes = 22)

# Create faceted plot
gg.manhattan(gwas_data, plot_title = "Manhattan Plot - Multiple Traits",
             facet_by_trait = TRUE)
```

## Input Data Requirements

Your data frame must contain:
- **CHR** - Chromosome number (numeric)
- **BP** - Base pair position (numeric)
- **P** - P-value (numeric, 0-1)

Optional:
- **SNP** - SNP identifier
- **trait** - Trait name (for multi-trait plots)

## Parameters

- `facet_by_trait` - Separate plots by trait (default: FALSE)
- `suggestiveline` - Suggestive p-value threshold (default: 1e-5)
- `genomewideline` - Genome-wide significance threshold (default: 5e-8)
- `logp` - Use -log10(p) transformation (default: TRUE)

See function documentation for additional options.
