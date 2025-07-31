# import packages
library(tidyverse)
library(ggpubr)
library(rstatix)

## figure 3A, OG glycoprotein fold change
# combine results for HEK293T, HepG2 and Jurkat cells
OG_glycoprotein_logFC_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  mutate(cell = "HepG2") |> 
  select(UniprotID, logFC, cell)

OG_glycoprotein_logFC_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  mutate(cell = "HEK293T") |> 
  select(UniprotID, logFC, cell)

OG_glycoprotein_logFC_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  mutate(cell = "Jurkat") |> 
  select(UniprotID, logFC, cell)

OG_glycoprotein_logFC_combined <- bind_rows(
  OG_glycoprotein_logFC_HepG2, 
  OG_glycoprotein_logFC_HEK293T, 
  OG_glycoprotein_logFC_Jurkat
)

# ks test
ks.test(OG_glycoprotein_Top_tb_HepG2$logFC, OG_glycoprotein_Top_tb_HEK293T$logFC)
ks.test(OG_glycoprotein_Top_tb_HEK293T$logFC, OG_glycoprotein_Top_tb_Jurkat$logFC)
ks.test(OG_glycoprotein_Top_tb_Jurkat$logFC, OG_glycoprotein_Top_tb_HepG2$logFC)

OG_glycoprotein_ks_test <- tribble(
  ~ .y., ~ group1, ~ group2, ~ p,
  "logFC", "HEK293T", "HepG2", 1.289e-11,
  "logFC", "Jurkat", "HEK293T", 0.001891,
  "logFC", "HepG2", "Jurkat", 1.651e-13
) |> 
  add_significance('p')

# violin boxplot
figure3A <- OG_glycoprotein_logFC_combined |>
  ggplot() +
  geom_violin(aes(x = cell, y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  stat_pvalue_manual(data = OG_glycoprotein_ks_test, label = "p.signif", tip.length = 0, 
                     size = 7, y.position = c(1.5, 1.9, 1.7)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "none"
  )

ggsave(
  filename = "figures/figure3/figure3A.eps",
  device = "eps",
  plot = figure3A, 
  height = 2, width = 2, units = "in"
)

## figure 3B, OG glycoprotein vs. WP protein
# extract information for HepG2
OG_WP_distribution_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  select(UniprotID, logFC_OG = logFC) |> 
  left_join(WP_protein_Top_tb_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_OG, logFC_WP = logFC) |> 
  filter(!is.na(logFC_WP)) |> 
  pivot_longer(cols = logFC_OG:logFC_WP, names_to = "Exp", values_to = "logFC") |> 
  mutate(Cell = "HepG2")

# extract information for HEK293T
OG_WP_distribution_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  select(UniprotID, logFC_OG = logFC) |> 
  left_join(WP_protein_Top_tb_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_OG, logFC_WP = logFC) |> 
  filter(!is.na(logFC_WP)) |> 
  pivot_longer(cols = logFC_OG:logFC_WP, names_to = "Exp", values_to = "logFC") |> 
  mutate(Cell = "HEK293T")

# extract information for Jurkat
OG_WP_distribution_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  select(UniprotID, logFC_OG = logFC) |> 
  left_join(WP_protein_Top_tb_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_OG, logFC_WP = logFC) |> 
  filter(!is.na(logFC_WP)) |> 
  pivot_longer(cols = logFC_OG:logFC_WP, names_to = "Exp", values_to = "logFC") |> 
  mutate(Cell = "Jurkat")

# combine results for HEK293T, HepG2 and Jurkat
OG_WP_distribution_combined <- bind_rows(
  OG_WP_distribution_HepG2,
  OG_WP_distribution_HEK293T,
  OG_WP_distribution_Jurkat
)

# ks test
ks.test(
  OG_WP_distribution_HepG2 |> filter(Exp == "logFC_OG") |> pull(logFC),
  OG_WP_distribution_HepG2 |> filter(Exp == "logFC_WP") |> pull(logFC)
)

ks.test(
  OG_WP_distribution_HEK293T |> filter(Exp == "logFC_OG") |> pull(logFC),
  OG_WP_distribution_HEK293T |> filter(Exp == "logFC_WP") |> pull(logFC)
)

ks.test(
  OG_WP_distribution_Jurkat |> filter(Exp == "logFC_OG") |> pull(logFC),
  OG_WP_distribution_Jurkat |> filter(Exp == "logFC_WP") |> pull(logFC)
)

ks_test_label <- tibble(
  Cell = c("HEK293T", "HEK293T", "HepG2", "HepG2", "Jurkat", "Jurkat"),
  p_signif = c("****", "****", "****", "****", "", ""),
  logFC = 1.5,
  Name = c("glycoprotein_HEK293T", "glycoprotein_HEK293T", "glycoprotein_HepG2", "glycoprotein_HepG2", "glycoprotein_Jurkat", "glycoprotein_Jurkat"),
  Exp = c("OG_HEK293T", "WP", "OG_HepG2", "WP", "OG_Jurkat", "WP")
)

# split violin plot
library(introdataviz)

figure3B <- OG_WP_distribution_combined |> 
  mutate(Name = paste("glycoprotein", Cell, sep = "_"), Exp = ifelse(Exp == "logFC_OG", paste("OG", Cell, sep = "_"), "WP")) |> 
  ggplot(aes(x = Name, y = logFC, fill = factor(Exp, c("OG_HEK293T", "OG_HepG2", "OG_Jurkat", "WP")))) +
  geom_split_violin(color = "transparent") +
  geom_text(data = ks_test_label, aes(y = logFC, label = p_signif), size = 7) +
  scale_fill_manual(values = c(
    "OG_HEK293T" = Color_2,
    "OG_HepG2" = Color_3,
    "OG_Jurkat" = Color_4,
    "WP" = "gray"
  )) +
  facet_grid(~ Cell, scale = "free_x") +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)"), fill = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.length.x = unit(0, "in"),
    legend.position = "none",
    strip.text = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(
  filename = 'figures/figure3/figure3B.eps',
  device = "eps",
  plot = figure3B, 
  height = 2, width = 2.5, units = 'in'
)

## figure 3C, GO common terms, OG glycoprotein abundance change
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

# gene ontology analysis
# HepG2
OG_glycoprotein_HepG2_GO <- enrichGO(
  gene = OG_glycoprotein_Top_tb_HepG2 |> pull(UniprotID),
  OrgDb = org.Hs.eg.db,
  keyType = 'UNIPROT',
  universe = OG_glycoprotein_total,
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_HepG2_GO@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_HepG2_GO.csv'
)

# HEK293T
OG_glycoprotein_HEK293T_GO <- enrichGO(
  gene = OG_glycoprotein_Top_tb_HEK293T |> pull(UniprotID),
  OrgDb = org.Hs.eg.db,
  keyType = 'UNIPROT',
  universe = OG_glycoprotein_total,
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_HEK293T_GO@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_HEK293T_GO.csv'
)

# Jurkat
OG_glycoprotein_Jurkat_GO <- enrichGO(
  gene = OG_glycoprotein_Top_tb_Jurkat |> pull(UniprotID),
  OrgDb = org.Hs.eg.db,
  keyType = 'UNIPROT',
  universe = OG_glycoprotein_total,
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_Jurkat_GO@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_Jurkat_GO.csv'
)

### common terms
# import gene ontology result for OG glycoprotein
OG_glycoprotein_HepG2_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_HepG2_GO.csv'
)

OG_glycoprotein_HEK293T_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_HEK293T_GO.csv'
)

OG_glycoprotein_Jurkat_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_Jurkat_GO.csv'
)

## RNA binding
# HepG2
RNA_binding_OG_glycoprotein_HepG2 <- OG_glycoprotein_HepG2_GO |> 
  filter(Description == 'RNA binding') |> 
  select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  pull()

OG_glycoprotein_HepG2_RNA_binding <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(UniprotID %in% RNA_binding_OG_glycoprotein_HepG2) |> 
  select(UniprotID, logFC) |> 
  mutate(cell = 'HepG2')

# HEK293T
RNA_binding_OG_glycoprotein_HEK293T <- OG_glycoprotein_HEK293T_GO |> 
  filter(Description == 'RNA binding') |> 
  select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  pull()

OG_glycoprotein_HEK293T_RNA_binding <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID %in% RNA_binding_OG_glycoprotein_HEK293T) |> 
  select(UniprotID, logFC) |> 
  mutate(cell = 'HEK293T')

# Jurkat
RNA_binding_OG_glycoprotein_Jurkat <- OG_glycoprotein_Jurkat_GO |> 
  filter(Description == 'RNA binding') |> 
  select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  pull()

OG_glycoprotein_Jurkat_RNA_binding <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(UniprotID %in% RNA_binding_OG_glycoprotein_Jurkat) |> 
  select(UniprotID, logFC) |> 
  mutate(cell = 'Jurkat')

## DNA binding
# HepG2
DNA_binding_OG_glycoprotein_HepG2 <- OG_glycoprotein_HepG2_GO |> 
  filter(Description == 'DNA binding') |> 
  select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  pull()

OG_glycoprotein_HepG2_DNA_binding <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(UniprotID %in% DNA_binding_OG_glycoprotein_HepG2) |> 
  select(UniprotID, logFC) |> 
  mutate(cell = 'HepG2')

# HEK293T
DNA_binding_OG_glycoprotein_HEK293T <- OG_glycoprotein_HEK293T_GO |> 
  filter(Description == 'DNA binding') |> 
  select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  pull()

OG_glycoprotein_HEK293T_DNA_binding <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID %in% DNA_binding_OG_glycoprotein_HEK293T) |> 
  select(UniprotID, logFC) |> 
  mutate(cell = 'HEK293T')

# Jurkat
DNA_binding_OG_glycoprotein_Jurkat <- OG_glycoprotein_Jurkat_GO |> 
  filter(Description == 'DNA binding') |> 
  select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  pull()

OG_glycoprotein_Jurkat_DNA_binding <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(UniprotID %in% DNA_binding_OG_glycoprotein_Jurkat) |> 
  select(UniprotID, logFC) |> 
  mutate(cell = 'Jurkat')

## nucleocytoplasmic transport
# HepG2
Nuc_Cyto_Transport_OG_glycoprotein_HepG2 <- OG_glycoprotein_HepG2_GO |> 
  filter(Description == 'nucleocytoplasmic transport') |> 
  select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  pull()

OG_glycoprotein_HepG2_Nuc_Cyto_Transport <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(UniprotID %in% Nuc_Cyto_Transport_OG_glycoprotein_HepG2) |> 
  select(UniprotID, logFC) |> 
  mutate(cell = 'HepG2')

# HEK293T
Nuc_Cyto_Transport_OG_glycoprotein_HEK293T <- OG_glycoprotein_HEK293T_GO |> 
  filter(Description == 'nucleocytoplasmic transport') |> 
  select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  pull()

OG_glycoprotein_HEK293T_Nuc_Cyto_Transport <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID %in% Nuc_Cyto_Transport_OG_glycoprotein_HEK293T) |> 
  select(UniprotID, logFC) |> 
  mutate(cell = 'HEK293T')

# Jurkat
Nuc_Cyto_Transport_OG_glycoprotein_Jurkat <- OG_glycoprotein_Jurkat_GO |> 
  filter(Description == 'nucleocytoplasmic transport') |> 
  select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  pull()

OG_glycoprotein_Jurkat_Nuc_Cyto_Transport <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(UniprotID %in% Nuc_Cyto_Transport_OG_glycoprotein_Jurkat) |> 
  select(UniprotID, logFC) |> 
  mutate(cell = 'Jurkat')

# combine results for RNA binding, DNA binding and nucleocytoplasmic transport
OG_glycoprotein_common_terms_combined <- bind_rows(
  # RNA binding
  bind_rows(
    OG_glycoprotein_HepG2_RNA_binding,
    OG_glycoprotein_HEK293T_RNA_binding,
    OG_glycoprotein_Jurkat_RNA_binding
  ) |> 
    mutate(
      Term = 'RNA binding'
    ),
  
  # DNA binding
  bind_rows(
    OG_glycoprotein_HepG2_DNA_binding,
    OG_glycoprotein_HEK293T_DNA_binding,
    OG_glycoprotein_Jurkat_DNA_binding
  ) |> 
    mutate(
      Term = 'DNA binding'
    ),
  
  # nucleoplasmic transport
  bind_rows(
    OG_glycoprotein_HepG2_Nuc_Cyto_Transport,
    OG_glycoprotein_HEK293T_Nuc_Cyto_Transport,
    OG_glycoprotein_Jurkat_Nuc_Cyto_Transport
  ) |> 
    mutate(
      Term = 'Nuc Cyto Transport'
    )
)

# ks test
# RNA binding
ks.test(OG_glycoprotein_HepG2_RNA_binding$logFC, OG_glycoprotein_HEK293T_RNA_binding$logFC)
ks.test(OG_glycoprotein_HEK293T_RNA_binding$logFC, OG_glycoprotein_Jurkat_RNA_binding$logFC)
ks.test(OG_glycoprotein_Jurkat_RNA_binding$logFC, OG_glycoprotein_HepG2_RNA_binding$logFC)

# DNA binding
ks.test(OG_glycoprotein_HepG2_DNA_binding$logFC, OG_glycoprotein_HEK293T_DNA_binding$logFC)
ks.test(OG_glycoprotein_HEK293T_DNA_binding$logFC, OG_glycoprotein_Jurkat_DNA_binding$logFC)
ks.test(OG_glycoprotein_Jurkat_DNA_binding$logFC, OG_glycoprotein_HepG2_DNA_binding$logFC)

# nucleocytoplasmic transport
ks.test(OG_glycoprotein_HepG2_Nuc_Cyto_Transport$logFC, OG_glycoprotein_HEK293T_Nuc_Cyto_Transport$logFC)
ks.test(OG_glycoprotein_HEK293T_Nuc_Cyto_Transport$logFC, OG_glycoprotein_Jurkat_Nuc_Cyto_Transport$logFC)
ks.test(OG_glycoprotein_Jurkat_Nuc_Cyto_Transport$logFC, OG_glycoprotein_HepG2_Nuc_Cyto_Transport$logFC)

# coming ks test results
library(ggpubr)
library(rstatix)

OG_glycoprotein_ks_test <- tribble(
  ~ .y., ~ group1, ~ group2, ~ p, ~ Term,
  "logFC", "HEK293T", "HepG2", 0.0007829, "RNA binding",
  "logFC", "Jurkat", "HEK293T", 0.0153, "RNA binding",
  "logFC", "HepG2", "Jurkat", 5.327e-05, "RNA binding",
  "logFC", "HEK293T", "HepG2", 0.0001767, "DNA binding",
  "logFC", "Jurkat", "HEK293T", 0.0304, "DNA binding",
  "logFC", "HepG2", "Jurkat", 0.0008344, "DNA binding",
  "logFC", "HEK293T", "HepG2", 0.0002324, "Nuc Cyto Transport",
  "logFC", "HepG2", "Jurkat", 0.0001177, "Nuc Cyto Transport"
) |> 
  add_significance('p')

# violin boxplot
figure3C <- OG_glycoprotein_common_terms_combined |> 
  ggplot() +
  geom_violin(aes(x = cell, y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  stat_pvalue_manual(data = OG_glycoprotein_ks_test, label = "p.signif", tip.length = 0, 
                     size = 7, hide.ns = "p.signif", y.position = c(1.5, 1.9, 1.7, 1.5)) +
  facet_grid(~ Term, scale = "free_x") +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    legend.title = element_text(size = 9, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(
  filename = 'figures/figure3/figure3C.eps',
  plot = figure3C,
  height = 2.5, width = 5, units = 'in'
)

### figure 3D, gene ontology analysis, generally up- and down- regulated O-GlcNAcylated proteins
## extract regulated O-GlcNAcylated proteins
# HepG2 up
OG_glycoprotein_up_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

# HepG2 down
OG_glycoprotein_down_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

# HepG2 median
OG_glycoprotein_median_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_HepG2$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_HepG2$UniprotID)) |> 
  select(UniprotID)

# HEK293T up
OG_glycoprotein_up_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

# HEK293T down
OG_glycoprotein_down_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

# HEK293T median
OG_glycoprotein_median_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_HEK293T$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_HEK293T$UniprotID)) |> 
  select(UniprotID)

# Jurkat up
OG_glycoprotein_up_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

# Jurkat down
OG_glycoprotein_down_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

# Jurkat median
OG_glycoprotein_median_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_Jurkat$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_Jurkat$UniprotID)) |> 
  select(UniprotID)

# generally upregulated O-GlcNAcylated proteins
# significantly upregulated in at least two types of cells and not downregulated in the other one
# downregulated OG glycoprotein list in all three cell types
OG_glycoprotein_lower_zero <- bind_rows(
  OG_glycoprotein_Top_tb_HepG2 |> 
    filter(logFC < 0) |> 
    select(UniprotID),
  
  OG_glycoprotein_Top_tb_HEK293T |> 
    filter(logFC < 0) |> 
    select(UniprotID),
  
  OG_glycoprotein_Top_tb_Jurkat |> 
    filter(logFC < 0) |> 
    select(UniprotID)
) |> 
  distinct()

# general upregulated O-GlcNAcylated protein list
general_up_OG_glycoprotein_list <- bind_rows(
  OG_glycoprotein_up_HepG2,
  OG_glycoprotein_up_HEK293T,
  OG_glycoprotein_up_Jurkat
) |> 
  count(UniprotID) |> 
  filter(n > 1) |> 
  filter(! UniprotID %in% OG_glycoprotein_lower_zero$UniprotID) |> 
  pull(UniprotID)

# generally downregulated O-GlcNAcylated proteins
# significantly downregulated in at least two types of cells and not upregulated in the other one
# upregulated OG glycoprotein list in all three cell types
OG_glycoprotein_higher_zero <- bind_rows(
  OG_glycoprotein_Top_tb_HepG2 |> 
    filter(logFC > 0) |> 
    select(UniprotID),
  
  OG_glycoprotein_Top_tb_HEK293T |> 
    filter(logFC > 0) |> 
    select(UniprotID),
  
  OG_glycoprotein_Top_tb_Jurkat |> 
    filter(logFC > 0) |> 
    select(UniprotID)
) |> 
  distinct()

# general downregulated O-GlcNAcylated protein list
general_down_OG_glycoprotein_list <- bind_rows(
  OG_glycoprotein_down_HepG2,
  OG_glycoprotein_down_HEK293T,
  OG_glycoprotein_down_Jurkat
) |> 
  count(UniprotID) |> 
  filter(n > 1) |> 
  filter(! UniprotID %in% OG_glycoprotein_higher_zero$UniprotID) |> 
  pull(UniprotID)

# gene ontology analysis, general up- or down-regulated O-GlcNAcylated protein
library(org.Hs.eg.db)
library(clusterProfiler)

# total OG glycoprotein
OG_glycoprotein_total <- bind_rows(
  OG_glycoprotein_Top_tb_HepG2 |> select(UniprotID),
  OG_glycoprotein_Top_tb_HEK293T |> select(UniprotID),
  OG_glycoprotein_Top_tb_Jurkat |> select(UniprotID)
) |> 
  distinct() |> 
  pull()

# gene ontology analysis
# general up
OG_glycoprotein_general_up_GO <- enrichGO(
  gene = general_up_OG_glycoprotein_list,
  OrgDb = org.Hs.eg.db,
  universe = OG_glycoprotein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_general_up_GO@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_general_up_GO.csv'
)

# general down
OG_glycoprotein_general_down_GO <- enrichGO(
  gene = general_down_OG_glycoprotein_list,
  OrgDb = org.Hs.eg.db,
  universe = OG_glycoprotein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_general_down_GO@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_general_down_GO.csv'
)

# import gene ontology results for general up- and down-regulated OG glycoprotein
OG_glycoprotein_general_up_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_general_up_GO.csv'
)

OG_glycoprotein_general_down_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_general_down_GO.csv'
)

# combine result with selected terms
OG_glycoprotein_general_up_down_GO <- bind_rows(
  OG_glycoprotein_general_up_GO |> 
    filter(
      Description %in% c(
        'catalytic activity, acting on RNA', 'nucleotide binding', 'ligase activity', 'amino acid metabolic process'
      )
    ) |> 
    select(Description, pvalue, Count) |> 
    mutate(Group = "up"),
  
    OG_glycoprotein_general_down_GO |> 
    filter(
      Description %in% c(
        'regulation of response to external stimulus', 'P-body', 'cytoplasmic ribonucleoprotein granule', 'post-transcriptional gene silencing'
      )
    ) |> 
    select(Description, pvalue, Count) |> 
    mutate(Group = "down")
)

# dot plot
figure3D <- OG_glycoprotein_general_up_down_GO |> 
  ggplot() +
  geom_point(
    aes(
      x = factor(Group, levels = c("up", "down")), 
      y = factor(
        Description, levels = c(
          'catalytic activity, acting on RNA', 'nucleotide binding', 'ligase activity', 'amino acid metabolic process',
          'regulation of response to external stimulus', 'P-body', 'cytoplasmic ribonucleoprotein granule', 'post-transcriptional gene silencing'
          )
      ), 
      size = Count, fill = pvalue
    ),
    shape = 21
  ) +
  scale_size(range = c(5, 11), breaks = c(4, 12)) +
  scale_fill_stepsn(limits = c(0, 1), 
                    breaks = c(0, 5E-3, 7.5E-3, 0.01, 0.05, 1), 
                    labels = c("0", "5E-3", "7.5E-3", "0.01", "0.05", "1"), 
                    n.breaks = 6, 
                    values = scales::rescale(c(0, 5E-3, 7.5E-3, 0.01, 0.05, 1)),
                    colours = c(Color_9, "transparent")) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.2),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black", lineheight = 0.1),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2)
  )

ggsave(
  filename = 'figures/figure3/figure3D.eps',
  device = cairo_ps,
  plot = figure3D,
  height = 2.5, width = 4.7, units = 'in',
  fallback_resolution = 1200
)




