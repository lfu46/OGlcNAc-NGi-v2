# import packages
library(tidyverse)

## figure 2A, UpSet plot for identified O-GlcNAcylated proteins 
library(UpSetR)

# generate total glycoprotein
OG_glycoprotein_total <- bind_rows(
  OG_glycoprotein_Top_tb_HepG2 |> distinct(UniprotID),
  OG_glycoprotein_Top_tb_HEK293T |> distinct(UniprotID),
  OG_glycoprotein_Top_tb_Jurkat |> distinct(UniprotID)
) |> 
  distinct()

OG_glycoprotein_upset_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> select(UniprotID) |> mutate(HepG2 = 1)
OG_glycoprotein_upset_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> select(UniprotID) |> mutate(HEK293T = 1)
OG_glycoprotein_upset_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> select(UniprotID) |> mutate(Jurkat = 1)

OG_glycoprotein_total_upset <- OG_glycoprotein_total |> 
  left_join(OG_glycoprotein_upset_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_upset_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_upset_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  mutate(
    HepG2 = ifelse(is.na(HepG2), 0, 1),
    HEK293T = ifelse(is.na(HEK293T), 0, 1),
    Jurkat = ifelse(is.na(Jurkat), 0, 1)
  )

OG_glycoprotein_total_upset_dataframe <- as.data.frame(OG_glycoprotein_total_upset)

upset(OG_glycoprotein_total_upset_dataframe, sets = c("HepG2", "HEK293T", "Jurkat"),
      main.bar.color = Color_2,
      sets.bar.color = Color_2,
      order.by = "freq", point.size = 5,
      text.scale = c(2, 1.5, 1.5, 1.5, 2, 3),
      mainbar.y.label = "# of quantified \nO-GlcNAcylated proteins")

## figure 2B, Venn Diagram of identified O-GlcNAcylation sites
library(eulerr)

# extract site information
OG_site_HepG2 <- OG_glycopeptide_Top_tb_HepG2 |> 
  filter(!is.na(combined_site)) |> 
  pull(Index)

OG_site_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  pull(Index)

OG_site_Jurkat <- OG_glycopeptide_Top_tb_Jurkat |> 
  filter(!is.na(combined_site)) |> 
  pull(Index)

# euler diagram
mat <- c(
  "HEK293T (1921)" = 1262,
  "HepG2 (382)" = 101,
  "Jurkat (1092)" = 471,
  "HEK293T (1921)&HepG2 (382)" = 68,
  "HEK293T (1921)&Jurkat (1092)" = 408,
  "HepG2 (382)&Jurkat (1092)" = 30,
  "HepG2 (382)&HEK293T (1921)&Jurkat (1092)" = 183
)

fit <- euler(mat)

tiff(filename = "figures/figure2/figure2B.tiff", 
     width = 2, height = 2, units = "in", compression = "lzw",
     res = 1200)

plot(fit, 
     fill = list(fill = c(Color_2, Color_3, Color_4), alpha = 0.8),
     edges = list(col = "white", lwd = 3),
     labels = list(alpha = c(0)),
     legend = list(fontsize = 0))

dev.off()

## figure 2C, gene ontology analysis for common and unique O-GlcNAcylated proteins
library(org.Hs.eg.db)
library(clusterProfiler)

# total OG protein
OG_glycoprotein_total <- bind_rows(
  OG_glycoprotein_Top_tb_HepG2 |> select(UniprotID),
  OG_glycoprotein_Top_tb_HEK293T |> select(UniprotID),
  OG_glycoprotein_Top_tb_Jurkat |> select(UniprotID)
) |> 
  distinct() |> 
  pull()

# unique OG protein
OG_glycoprotein_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> distinct(UniprotID)
OG_glycoprotein_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> distinct(UniprotID)
OG_glycoprotein_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> distinct(UniprotID)

OG_glycoprotein_unique_HepG2 <- OG_glycoprotein_HepG2 |> 
  filter(!(UniprotID %in% OG_glycoprotein_HEK293T$UniprotID) & !(UniprotID %in% OG_glycoprotein_Jurkat$UniprotID)) |> 
  pull()

OG_glycoprotein_unique_HEK293T <- OG_glycoprotein_HEK293T |> 
  filter(!(UniprotID %in% OG_glycoprotein_HepG2$UniprotID) & !(UniprotID %in% OG_glycoprotein_Jurkat$UniprotID)) |> 
  pull()

OG_glycoprotein_unique_Jurkat <- OG_glycoprotein_Jurkat |> 
  filter(!(UniprotID %in% OG_glycoprotein_HepG2$UniprotID) & !(UniprotID %in% OG_glycoprotein_HEK293T$UniprotID)) |> 
  pull()

# common OG protein
OG_glycoprotein_common <- OG_glycoprotein_HepG2 |> 
  filter((UniprotID %in% OG_glycoprotein_HEK293T$UniprotID) & (UniprotID %in% OG_glycoprotein_Jurkat$UniprotID)) |> 
  pull()

# gene ontology analysis
OG_glycoprotein_unique_HepG2_GO <- enrichGO(
  gene = OG_glycoprotein_unique_HepG2,
  OrgDb = org.Hs.eg.db,
  universe = OG_glycoprotein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_unique_HepG2_GO@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_unique_HepG2_GO.csv'
)

OG_glycoprotein_unique_HEK293T_GO <- enrichGO(
  gene = OG_glycoprotein_unique_HEK293T,
  OrgDb = org.Hs.eg.db,
  universe = OG_glycoprotein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_unique_HEK293T_GO@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_unique_HEK293T_GO.csv'
)

OG_glycoprotein_unique_Jurkat_GO <- enrichGO(
  gene = OG_glycoprotein_unique_Jurkat,
  OrgDb = org.Hs.eg.db,
  universe = OG_glycoprotein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_unique_Jurkat_GO@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_unique_Jurkat_GO.csv'
)

OG_glycoprotein_common_GO <- enrichGO(
  gene = OG_glycoprotein_common,
  OrgDb = org.Hs.eg.db,
  universe = OG_glycoprotein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_common_GO@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_common_GO.csv'
)

# import GO result
library(tidyverse)

OG_glycoprotein_common_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_common_GO.csv'
)
OG_glycoprotein_unique_HEK293T_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_unique_HEK293T_GO.csv'
)
OG_glycoprotein_unique_HepG2_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_unique_HepG2_GO.csv'
)
OG_glycoprotein_unique_Jurkat_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_unique_Jurkat_GO.csv'
)

# combine common and unique GO terms
gene_ontology_OG_glycoprotein_common_unique <- bind_rows(
  OG_glycoprotein_common_GO |> filter(
    Description %in% c(
      'RNA binding', 'cytoplasmic ribonucleoprotein granule', 'P-body', 'nucleocytoplasmic transport', 'DNA binding'
    )
  ) |> 
    mutate(group = 'Common'),
  
  OG_glycoprotein_unique_HEK293T_GO |> filter(
    Description %in% c(
      'positive regulation of canonical Wnt signaling pathway', 'double-stranded DNA binding', 'humoral immune response'
    )
  ) |> 
    mutate(group = 'HEK293T only'),
  
  OG_glycoprotein_unique_HepG2_GO |> filter(
    Description %in% c(
      'lipid metabolic process', 'sterol metabolic process', 'cholesterol homeostasis', 'cell-matrix adhesion'
    )
  ) |> 
    mutate(group = 'HepG2 only'),
  
  OG_glycoprotein_unique_Jurkat_GO |> filter(
    Description %in% c(
      'immune response', 'positive regulation of T cell proliferation', 'cell surface', 'plasma membrane'
    )
  ) |> 
    mutate(group = 'Jurkat only')
) |> 
  select(Description, pvalue, group) |> 
  mutate(
    Description = factor(Description, levels = c(
      'RNA binding', 'cytoplasmic ribonucleoprotein granule', 'P-body', 'nucleocytoplasmic transport', 'DNA binding',
      'positive regulation of canonical Wnt signaling pathway', 'double-stranded DNA binding', 'humoral immune response',
      'lipid metabolic process', 'sterol metabolic process', 'cholesterol homeostasis', 'cell-matrix adhesion',
      'immune response', 'positive regulation of T cell proliferation', 'cell surface', 'plasma membrane'
    ))
  )

# heatmap
figure2C <- gene_ontology_OG_glycoprotein_common_unique |> 
  ggplot() +
  geom_tile(
    aes(x = group, y = Description, fill = pvalue)
  ) +
  labs(fill = 'P.Value') +
  scale_fill_stepsn(limits = c(0, 1),
                    breaks = c(0, 1E-3, 0.01, 0.05, 1),
                    labels = c("0", "1E-3", "0.01", "0.05", "1"),
                    n.breaks = 4,
                    values = scales::rescale(c(0, 1E-3, 0.01 ,0.05, 1)),
                    colours = c(Color_9, NA)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0),
    panel.grid.minor = element_line(color = "gray", linewidth = 0),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black", lineheight = 0.6),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2),
    legend.key.size = unit(0.15, "in")
  )

ggsave(
  filename = 'figures/figure2/figure2C.eps',
  device = cairo_ps,
  plot = figure2C, 
  height = 3.5, width = 6, 
  units = "in", fallback_resolution = 1200
)

