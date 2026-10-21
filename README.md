# Manhattan plots with ggplot2

Plot GWAS results for one or more chromosomes, with optional panels for each trait.
The repository also includes a reproducible example-data generator.

![Example Manhattan plot](Rplot.png)

The default Ocean style uses alternating blue and teal chromosomes, pale trait
headers, subtle horizontal gridlines, and gold points at or above the selected
genome-wide significance threshold. The image uses simulated data.

## Setup

Run the examples and tests from the repository root. Install ggplot2 once:

```r
install.packages("ggplot2")
```

The simulator itself uses base R. Plotting and the regression checks need ggplot2.
The current checks have been run with R 4.5.2 and ggplot2 4.0.2.

## Quick start

```r
library(ggplot2)
source("gg.manhattan.R")
source("simulate_gwas_data.R")

set.seed(42)
gwas_data <- simulate_gwas_data(
  n_snps = 40000, n_traits = 3, n_chromosomes = 22
)
stopifnot(nrow(gwas_data) == 120000)

plot <- gg.manhattan(
  gwas_data,
  facet_by_trait = TRUE,
  plot_title = "Manhattan plot by trait"
)
print(plot)
ggsave("manhattan.png", plot = plot, width = 10, height = 8, dpi = 150)
```

Sourcing the simulator only defines the function: it does not generate data,
reset the random seed, or run the plotting examples. To run the complete example
workflow instead, use:

```sh
Rscript simulate_gwas_data.R
```

The script prints a data summary and two plots; in a noninteractive R session,
the default graphics device normally writes them to `Rplots.pdf`.

## Input data

Pass a data frame to `gg.manhattan()`. Column names can be overridden with
`chr`, `bp`, `p`, `snp`, and `trait`.

| Default column | Contents |
| --- | --- |
| `CHR` | Numeric chromosome identifiers; encode X, Y, or MT as numbers before plotting. |
| `BP` | Numeric base-pair positions within each chromosome. |
| `P` | Finite p-values with `0 < P <= 1` when `logp = TRUE`. |
| `SNP` | Optional SNP identifiers. |
| `trait` | Trait labels required for panels by trait; numeric or character labels are accepted. |

Supply valid, nonmissing chromosome identifiers, positions, and p-values.
With `facet_by_trait = TRUE` and at least two distinct traits, each trait gets
a separate panel. Otherwise the function uses its ordinary plot layout.

For one chromosome, the x-axis shows positions in Mb. For multiple chromosomes,
it uses cumulative positions with chromosome labels.

## Significance thresholds

Threshold arguments are **plotted y-coordinates**, not raw p-values. With the
default `logp = TRUE`, convert a p-value threshold using `-log10(p)`:

```r
threshold_plot <- gg.manhattan(
  gwas_data,
  facet_by_trait = TRUE,
  suggestiveline = -log10(1e-4),   # y = 4
  genomewideline = -log10(1e-6)   # y = 6
)
print(threshold_plot)
```

The defaults are `suggestiveline = 5` and
`genomewideline = -log10(5e-8)` (approximately 7.30103).
Pass `FALSE` to either argument to hide that line.

The suggestive line is dotted and the genome-wide line is dashed. Gold
highlights follow `genomewideline`: by default they mark `P <= 5e-8`, including
points exactly on the threshold. Set `highlight = FALSE` to turn off the gold
overlay, or set `highlight_colour` to choose another colour. Hiding the
genome-wide line also disables the automatic highlights. In log-p plots the
y-axis includes enabled thresholds even if all observations fall below them.

With `logp = FALSE`, `P` is plotted directly and thresholds still use
y-coordinates. For a score plot, disable the p-value threshold lines and set an
appropriate y-axis label with `+ labs(y = "Score")`.
Automatic significance highlights are disabled when `logp = FALSE`.

## Ocean styling

```r
ocean_plot <- gg.manhattan(
  gwas_data,
  facet_by_trait = TRUE,
  highlight = TRUE,
  highlight_colour = "#C77B16",
  point_size = 0.6,
  point_alpha = 0.75
)
print(ocean_plot)
```

Use `cols` to change chromosome colours. For an un-faceted regional plot with
multiple traits, the points instead use trait colours drawn from that palette.
The highlight overlay uses slightly larger opaque points so significant
associations remain visible in dense plots. You can extend the returned ggplot
with `labs()` or `theme()` as usual.

To reproduce the README figure from simulated data:

```sh
Rscript examples/ocean.R
```

## Plotting arguments

| Argument | Default | Meaning |
| --- | --- | --- |
| `x` | Required | Input data frame. |
| `chr`, `bp`, `p` | `"CHR"`, `"BP"`, `"P"` | Columns containing chromosomes, positions, and p-values. |
| `snp`, `trait` | `"SNP"`, `"trait"` | Optional SNP and trait columns. |
| `cols` | `c("#294F79", "#278D91")` | Alternating Ocean blue and teal chromosome colours. |
| `plot_title` | `NULL` | Optional plot title. |
| `chrlabs` | `NULL` | Labels for multi-chromosome plots, in sorted chromosome order; provide one per chromosome. |
| `suggestiveline` | `-log10(1e-5)` | Suggestive threshold at y = 5; `FALSE` hides it. |
| `genomewideline` | `-log10(5e-8)` | Genome-wide threshold at y approximately 7.30103; `FALSE` hides it. |
| `logp` | `TRUE` | Transform `P` to `-log10(P)`; `FALSE` plots the input values directly. |
| `facet_by_trait` | `FALSE` | Stack panels by trait when two or more traits are present. |
| `highlight` | `TRUE` | Highlight points at or above `genomewideline` in log-p plots. |
| `highlight_colour` | `"#C77B16"` | Gold colour for significant points. |
| `point_size` | `0.6` | Base point size; highlighted points are 1.5 times this size. |
| `point_alpha` | `0.75` | Base point opacity; highlighted points remain opaque. |

The return value is a ggplot object, so it can be printed, extended with
ggplot2 layers, or saved with `ggsave()`.

## Simulated data

`simulate_gwas_data(n_snps = 40000, n_traits = 3, n_chromosomes = 22,
n_peaks = 3, peak_width = 3e6, peak_effect = 7)`
returns columns `CHR`, `BP`, `P`, `SNP`, and `trait`.

- All three counts must be single positive whole numbers.
- `n_snps` is the exact number of distinct loci **per trait**. The output has
  `n_snps * n_traits` rows, including when the SNP count is not divisible by
  the chromosome count.
- SNPs are distributed as evenly as possible. The first chromosomes receive
  any remainder; chromosomes with no SNPs are omitted. If there are fewer
  SNPs than chromosomes, fewer chromosomes will appear in the output.
- Each locus keeps the same SNP ID, chromosome, and position across traits.
  Positions are sampled without replacement from 1 to 250,000,000 bp within
  each chromosome.
- P-values are derived from two-sided standard-normal tests, with correlated
  background statistics for nearby SNPs. Away from peaks, the marginal null
  p-values are uniform rather than truncated.
- By default each trait has three regional signals. Effects taper with distance
  from each centre, producing groups of associated neighboring SNPs. One centre
  is shared across traits; the other centres and effect strengths vary. Centres
  use different chromosomes when possible, otherwise they prefer separated loci.
- Call `set.seed()` before generating data to reproduce an example.

| Simulator argument | Default | Meaning |
| --- | --- | --- |
| `n_snps` | `40000` | Exact number of distinct loci per trait. |
| `n_traits` | `3` | Number of traits. |
| `n_chromosomes` | `22` | Number of chromosomes over which loci are distributed. |
| `n_peaks` | `3` | Regional signal centres per trait, capped at the available SNP count; `0` generates background only. |
| `peak_width` | `3e6` | Gaussian effect-profile standard deviation in base pairs; larger values produce broader peaks. |
| `peak_effect` | `7` | Approximate absolute effect on the standard-normal statistic at a peak centre, varied by 10% between peaks; larger values tend to give stronger associations. |

The count arguments must be whole numbers; `n_peaks` may be zero, and the other
counts must be positive. Width and effect must be finite positive numbers.
Background correlation decays exponentially with physical distance, with a
decay length of `peak_width / 6`. Changing width therefore changes both the
signal footprint and local noise correlation. The actual peak height varies
with the simulated noise and marker density; sparse datasets may not contain
enough neighboring markers to show a smooth regional peak.

`attr(gwas_data, "peaks")` records each simulated centre's trait, chromosome,
position, SNP ID, signed effect, and width. It describes the injected signals,
not detected lead SNPs. The example image is generated directly by this
simulator without adding separate artificial p-value peaks in the plotting code.

```r
set.seed(42)
small_data <- simulate_gwas_data(n_snps = 10, n_traits = 2)
stopifnot(
  nrow(small_data) == 20,
  length(unique(small_data$SNP)) == 10
)

single_trait <- subset(gwas_data, trait == "Trait_1")
print(gg.manhattan(single_trait, plot_title = "Trait 1"))

regional_data <- subset(gwas_data, CHR == 1)
print(gg.manhattan(regional_data, facet_by_trait = TRUE))
```

These are illustrative summary statistics with LD-like local correlation,
not a genotype-based GWAS or a calibrated model of linkage disequilibrium,
population structure, heritability, or genetic correlations between traits.

## Regression checks

```sh
Rscript tests/regression.R
```

The checks cover regional facet axes, chromosome labels, shared SNP coordinates,
independent trait p-values, small and nondivisible SNP counts, invalid arguments,
source-time behaviour, threshold coordinates, and default-size plots. Style
checks also cover exact-threshold highlights, custom and disabled highlights,
raw-value plots, nonconsecutive chromosome labels, and visible threshold lines.
Simulation checks cover regional signal clusters, null p-value calibration,
peak controls, and shared centre metadata.

## Files

| File | Purpose |
| --- | --- |
| `gg.manhattan.R` | Plotting function. |
| `simulate_gwas_data.R` | Simulator and examples run with Rscript. |
| `tests/regression.R` | Standalone regression checks. |
| `examples/ocean.R` | Reproducible Ocean figure using synthetic association peaks. |
| `Rplot.png` | Example output. |
| `LICENSE.md` | GPL-3.0 licence. |
