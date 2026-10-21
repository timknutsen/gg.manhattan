gg.manhattan <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", trait = "trait",
                         cols=c("#294F79", "#278D91"), plot_title=NULL, chrlabs=NULL,
                         suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8),
                         logp=TRUE, facet_by_trait=FALSE, highlight=TRUE,
                         highlight_colour="#C77B16", point_size=0.6, point_alpha=0.75) {
  
  # Require ggplot2, install if not present in labrary.
  if (!require(ggplot2)) install.packages("ggplot2")
  
  # Not sure why, but package check will warn without this.
  CHR=BP=P=index=NULL
  
  # Check for sensible dataset
  ## Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  ## warn if you don't have a snp column
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!(trait %in% names(x))) warning(paste("No trait column found. OK unless you're trying to plot multiple traits."))
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # Create a new data.frame with columns called CHR, BP, and P.
  d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])
  
  # If the input data frame has a SNP column, add it to the new data frame you're creating.
  #	if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
  # If the input data frame has a trait column, add it to the new data frame you're creating.
  if (!is.null(x[[trait]])) d=transform(d, TRAIT=x[[trait]])
  
  # Set positions, ticks, and labels for plotting
  ## Sort and keep only values where is numeric.
  #d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }
  d$pos=NA
  
  
  # Fixes the bug where one chromosome is missing by adding a sequential index column.
  d$index=NA
  ind = 0
  for (i in unique(d$CHR)){
    ind = ind + 1
    d[d$CHR==i,]$index = ind
  }
  
  # This section sets up positions and ticks. Ticks should be placed in the
  # middle of a chromosome. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromsome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    ## Uncomment the next line to plot single chr results in Mb
    d$pos=d$BP/1e6
    ticks=floor(length(d$pos))/2+1
    xlabel = paste('Chromosome',unique(d$CHR),'position(Mb)')
    labs = ticks
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
        d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
      }
      # Old way: assumes SNPs evenly distributed
      # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
      # New way: doesn't make that assumption
      ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2)
    }
    xlabel = 'Chromosome'
    ## Set x labels to Chromsomes, or provided custom list from chrlabs argument.
    if (!is.null(chrlabs)) {
      labs <- chrlabs
    } else {
      labs <- unique(d$CHR)
    }
  }
  
  # Build a common plot so regional, faceted, and genome-wide views share a style.
  overlay_traits <- nchr == 1 && !facet_by_trait &&
    !is.null(d$TRAIT) && length(unique(d$TRAIT)) >= 2
  colour_mapping <- if (overlay_traits) aes(colour = factor(TRAIT)) else aes(colour = factor(CHR))
  palette <- if (overlay_traits) {
    grDevices::colorRampPalette(cols)(length(unique(d$TRAIT)))
  } else rep(cols, length.out = nchr)
  plot <- ggplot(d, aes(x = pos, y = logp)) +
    geom_point(colour_mapping, size = point_size, alpha = point_alpha, stroke = 0) +
    scale_colour_manual(values = palette, guide = "none") +
    labs(x = xlabel, y = if (logp) expression(-log[10](italic(p))) else p,
         title = plot_title)

  if (nchr == 1) {
    plot <- plot + scale_x_continuous(name = xlabel,
                                      expand = expansion(mult = c(0.008, 0.008)))
  } else {
    plot <- plot + scale_x_continuous(name = xlabel, breaks = ticks, labels = labs,
                                      expand = expansion(mult = c(0.008, 0.008)))
  }

  # Include the threshold lines in the limits even when no points reach them.
  thresholds <- numeric()
  if (suggestiveline) thresholds <- c(thresholds, suggestiveline)
  if (genomewideline) thresholds <- c(thresholds, genomewideline)
  if (logp) {
    ymax <- max(2, ceiling(max(c(d$logp, thresholds)) / 2) * 2)
    plot <- plot + scale_y_continuous(
      breaks = seq(0, ymax, 2), limits = c(0, ymax),
      expand = expansion(mult = c(0, 0.025)))
  }

  if (suggestiveline) {
    plot <- plot + geom_hline(yintercept = suggestiveline, colour = "#617887",
                              linewidth = 0.35, linetype = "dotted")
  }
  if (genomewideline) {
    plot <- plot + geom_hline(yintercept = genomewideline, colour = "#617887",
                              linewidth = 0.45, linetype = "dashed")
  }
  # Reuse the selected genome-wide threshold. Raw-value plots are not
  # automatically interpreted as p-values, so they do not receive highlights.
  if (highlight && logp && genomewideline) {
    hits <- d[d$logp >= genomewideline, , drop = FALSE]
    if (nrow(hits)) {
      plot <- plot + geom_point(data = hits, colour = highlight_colour,
                                size = point_size * 1.5, alpha = 1, stroke = 0)
    }
  }

  if (facet_by_trait && !is.null(d$TRAIT) && length(unique(d$TRAIT)) >= 2) {
    plot <- plot + facet_wrap(~TRAIT, ncol = 1)
  }

  plot + theme_minimal(base_size = 12, base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      text = element_text(colour = "#123348"),
      axis.text = element_text(colour = "#526B79", size = 10),
      axis.title = element_text(colour = "#123348", size = 12),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.line.x = element_line(colour = "#E4EBEE", linewidth = 0.35),
      axis.ticks.x = element_line(colour = "#E4EBEE", linewidth = 0.35),
      panel.grid.major.y = element_line(colour = "#E4EBEE", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.y = grid::unit(10, "pt"),
      strip.background = element_rect(fill = "#EDF4F6", colour = NA),
      strip.text = element_text(colour = "#123348", face = "bold", hjust = 0,
                                size = 11, margin = margin(5, 8, 5, 8)),
      plot.title = element_text(size = 20, face = "bold", margin = margin(b = 12)),
      plot.subtitle = element_text(colour = "#526B79", margin = margin(b = 12)),
      plot.caption = element_text(colour = "#526B79", hjust = 0, size = 9, lineheight = 1.4,
                                   margin = margin(t = 12)),
      plot.title.position = "plot", plot.caption.position = "plot",
      plot.margin = margin(16, 20, 12, 14),
      legend.position = "none"
    )
}
