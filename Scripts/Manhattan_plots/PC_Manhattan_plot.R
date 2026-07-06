# This script generates a Manhattan plot from GWAS summary statistics.
# It visualises -log10(p-values) across genomic positions for each chromosome.
# Genome-wide significant SNPs are highlighted.


#load the packages 

library(dplyr)
library(ggplot2)

#prepare GWAS results to clean the data, convert chromosome position to numeric values and identifies genome wode significant SNPs.
prepare_gwas <- function(df) {

  df <- df %>%
    mutate(
      CHR = as.numeric(CHR),              # chromosome as numeric
      POS = as.numeric(POS),              # base pair position
      P = pmax(P, 1e-300),                
      logP = -log10(P),                   # convert p-values to -log10 scale
      Significant = P < 5e-8              # genome wide significant threshold
    ) %>%
    arrange(CHR, POS)                     # sort by chromosome and position


  # compute chromosome lengths and offsets to place chromosomes end to end on x-axis
 
  chrom_sizes <- df %>%
    group_by(CHR) %>%
    summarise(chr_len = max(POS, na.rm = TRUE), .groups = "drop") %>%
    mutate(offset = cumsum(chr_len) - chr_len)

  # Add cumulative base-pair position across chromosomes
  
  df <- df %>%
    left_join(chrom_sizes, by = "CHR") %>%
    mutate(BPcum = POS + offset)

  # Add cumulative base-pair position across chromosomes
  
  axis_df <- df %>%
    group_by(CHR) %>%
    summarise(center = (max(BPcum) + min(BPcum)) / 2,
              .groups = "drop")

  list(data = df, axis = axis_df)
}

# add function to PC GWAS data - must be GWAS summary statistics dataframe
pc_prep <- prepare_gwas(pc_all)

pc <- pc_prep$data
axis_pc <- pc_prep$axis

# Subsample points for faster plotting (reduces overplotting)

pc_plot <- pc %>% slice(seq(1, n(), by = 5))

# create the manhattan plot
 p_pc <- ggplot(pc_plot, aes(x = BPcum, y = logP)) +

  # Alternate chromosome colours for visual clarity
  geom_point(aes(color = factor(CHR %% 2)),
             size = 0.6, alpha = 0.8) +

  # Highlight genome-wide significant SNPs in red
  geom_point(data = subset(pc_plot, Significant),
             colour = "red3", size = 0.9) +
  # Genome-wide significance threshold line
  geom_hline(yintercept = -log10(5e-8),
             linetype = "dashed") +

  # Chromosome labels on x-axis
  scale_x_continuous(breaks = axis_pc$center,
                     labels = axis_pc$CHR) +

  scale_color_manual(values = c("grey20", "grey70")) +

  labs(title = "Prostate Cancer GWAS",
       x = "Chromosome",
       y = expression(-log[10](italic(P)))) +

  theme_classic(base_size = 14) +

  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))


pc_plot <- pc %>% slice(seq(1, n(), by = 10))

#save the plot
ggsave(
  filename = "/scratch/prj/dtr_incom/sanna/manhattan_pc.pdf",
  plot = p_pc,
  width = 12,
  height = 6
)
