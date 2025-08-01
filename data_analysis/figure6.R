# import packages
library(tidyverse)

### figure 6A, isoelectric point, up- and down-regulated O-GlcNAcylation site
# remove Tn antigen OG glycoprotein from sequence parameter database
Tn_antigen_database <- read_csv(
  'data_source/Tn_antigen_database/Tn_antigen_database.csv'
)

library(readxl)

OG_site_Mer31_sequence_parameter_filtered <- read_xlsx(
  'data_source/OGlcNAc_site/OG_site_Mer31_regulated_combined_sequence_parameter.xlsx'
) |> 
  filter(
    ! UniprotID %in% Tn_antigen_database$Uniprot.Entry
  )

# wilcoxon test
library(rstatix)
library(ggpubr)

isoelectric_point_wilcoxon_test <- OG_site_Mer31_sequence_parameter_filtered |> 
  filter(category %in% c('up', 'down')) |> 
  group_by(cell) |> 
  wilcox_test(Isoelectric_point ~ category, p.adjust.method = 'BH') |> 
  add_significance('p')

# violin boxplot
figure6A <- OG_site_Mer31_sequence_parameter_filtered |> 
  filter(category %in% c('up', 'down')) |> 
  ggplot() +
  geom_violin(
    aes(x = factor(category, levels = c("up", "down")), y = Isoelectric_point, fill = cell), 
    color = "transparent"
  ) +
  geom_boxplot(
    aes(x = factor(category, levels = c("up", "down")), y = Isoelectric_point), color = "black",
    outliers = FALSE, width = 0.2
  ) +
  labs(x = "", y = "isoelectric point") +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ cell, scales = "free_x") +
  stat_pvalue_manual(data = isoelectric_point_wilcoxon_test |> filter(), label = "p.signif", tip.length = 0, 
                     size = 7,
                     hide.ns = "p", y.position = c(14)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(
  filename = 'figures/figure6/figure6A.eps',
  plot = figure6A,
  height = 2.5, width = 2.5, units = 'in'
)

### figure 6B, Buried vs. Exposed, up- and down-regulated O-GlcNAcylation site
# remove Tn antigen OG glycoprotein from sequence parameter database
Tn_antigen_database <- read_csv(
  'data_source/Tn_antigen_database/Tn_antigen_database.csv'
)

library(readxl)

OG_site_BuriedExposed_filtered <- read_xlsx(
  'data_source/OGlcNAc_site/OG_glycopeptide_localized_singlesite_secondary_structure_combined.xlsx'
) |> 
  filter(
    ! UniprotID %in% Tn_antigen_database$Uniprot.Entry
  )

# wilcoxon test
library(rstatix)
library(ggpubr)

OG_site_BuriedExposed_wilcoxon_test <- OG_site_BuriedExposed_filtered |> 
  filter(!is.na(Class.Assignment)) |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ Class.Assignment, p.adjust.method = "BH") |> 
  add_significance("p")

# violin boxplot
figure6B <- OG_site_BuriedExposed_filtered |> 
  filter(!is.na(Class.Assignment)) |> 
  ggplot() +
  geom_violin(
    aes(x = factor(Class.Assignment, levels = c("E", "B")), y = logFC, fill = cell), 
    color = "transparent"
  ) +
  geom_boxplot(
    aes(x = factor(Class.Assignment, levels = c("E", "B")), y = logFC), 
    color = "black", width = 0.2, outliers = FALSE
  ) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ cell, scales = "free_x") +
  scale_x_discrete(labels = c("E" = "Exposed",
                              "B" = "Buried")) +
  stat_pvalue_manual(data = OG_site_BuriedExposed_wilcoxon_test, 
                     label = "p.signif", tip.length = 0, size = 7,
                     hide.ns = "p", y.position = c(1.8)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1, lineheight = 0.1),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(
  filename = 'figures/figure6/figure6B.eps',
  plot = figure6B,
  height = 2.8, width = 2.5, units = 'in'
)



