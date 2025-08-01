# import packages
library(tidyverse)

### figure 5A, extracellular exosoem, nucleus and mitochondrion, in HepG2
library(org.Hs.eg.db)
library(clusterProfiler)
library(Biostrings)

# import human fasta file
human_fasta <- readAAStringSet(
  'data_source/fasta_file/uniprotkb_reviewed_true_AND_model_organ_2025_08_01.fasta'
)

# build tibble using human fasta
human_fasta_tibble <- tibble(
  Name = names(human_fasta),
  Sequence = as.character(human_fasta),
  Length = width(human_fasta)
) |> 
  mutate(
    Name = sub(' .*', '', Name),
    Full_Protein_Length = as.numeric(Length)
  ) |> 
  separate(Name, into = c('sp', 'UniProt_Accession', 'name'), sep = '\\|') |> 
  select(UniProt_Accession, Sequence, Full_Protein_Length)

# gene ontology analysis, OG glycoprotein HepG2, background: human proteome
OG_glycoprotein_HepG2_GO_bghp <- enrichGO(
  gene = OG_glycoprotein_Top_tb_HepG2 |> pull(UniprotID),
  OrgDb = org.Hs.eg.db,
  keyType = 'UNIPROT',
  universe = human_fasta_tibble$UniProt_Accession,
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_HepG2_GO_bghp@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_HepG2_GO_bghp.csv'
)

# extract extracellular exosome, nucleus and mitochondrion, OG glycoprotein list
# use quantified total OG glycoprotein in this experiment as background
# different background, different enrichment result
OG_glycoprotein_HepG2_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_HepG2_GO.csv'
)

# extracellular exosome
exo_OG_list <- OG_glycoprotein_HepG2_GO |> 
  filter(
    Description %in% 'extracellular exosome'
  ) |> 
  select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  pull()

# nucleus
nuc_OG_list <- OG_glycoprotein_HepG2_GO |> 
  filter(
    str_detect(Description, 'nucleus')
  ) |> 
  select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  pull()

# mitochondrion list
mito_OG_list <- OG_glycoprotein_HepG2_GO |> 
  filter(
    Description %in% 'mitochondrion'
  ) |> 
  select(geneID) |> 
  separate_rows(geneID, sep = '/') |> 
  pull()

# combine exo, nuc and mito OG glycoprotein in HepG2
OG_glycoprotein_HepG2_exo_nuc_mito <- bind_rows(
  OG_glycoprotein_Top_tb_HepG2 |> 
    filter(
      UniprotID %in% exo_OG_list
    ) |> 
    mutate(
      group = 'Exosome'
    ),
  
  OG_glycoprotein_Top_tb_HepG2 |> 
    filter(
      UniprotID %in% nuc_OG_list
    ) |> 
    mutate(
      group = 'Nucleus'
    ),
  
  OG_glycoprotein_Top_tb_HepG2 |> 
    filter(
      UniprotID %in% mito_OG_list
    ) |> 
    mutate(
      group = 'Mitochondrion'
    )
)

# wilcoxon test
library(rstatix)
library(ggpubr)

Exo_Nuc_Mito_wilcoxon_test <- OG_glycoprotein_HepG2_exo_nuc_mito |> 
  wilcox_test(logFC ~ group, p.adjust.method = 'BH') |> 
  add_significance('p')

# violin boxplot
figure5A <- OG_glycoprotein_HepG2_exo_nuc_mito |> 
  ggplot() +
  geom_violin(
    aes(
      x = factor(group, levels = c("Exosome", "Nucleus", "Mitochondrion")), 
      y = logFC, fill = factor(group, levels = c("Exosome", "Nucleus", "Mitochondrion"))
    ), 
    color = "transparent"
  ) +
  geom_boxplot(
    aes(
      x = factor(group, levels = c("Exosome", "Nucleus", "Mitochondrion")), 
      y = logFC
    ), 
    color = "black", outliers = FALSE, width = 0.2
  ) +
  labs(x = "", y = "") +
  scale_fill_manual(
    values = c(
      "Exosome" = Color_2,
      "Nucleus" = Color_3,
      "Mitochondrion" = Color_4
    )
  ) +
  stat_pvalue_manual(data = Exo_Nuc_Mito_wilcoxon_test, 
                     label = "p.signif", tip.length = 0, size = 7,
                     hide.ns = "p", y.position = c(2.0, 2.2)) +
  coord_cartesian(ylim = c(-1.8, 2.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "none"
  )

ggsave(
  filename = 'figures/figure5/figure5A.eps', 
  plot = figure5A, 
  height = 2, width = 2, units = "in"
)

### figure 5B, GSEA, OG glycoprotein HepG2
# CORUM database
CORUM_database <- read_delim(
  "data_source/CORUM/corum_allComplexes.txt",
  col_names = TRUE,
  name_repair = "universal"
) |> 
  filter(organism == "Human") |> 
  select(complex_name, UniprotID = subunits_uniprot_id) |> 
  separate_rows(UniprotID, sep = ";")

# OG glycoprotein HepG2, protein complex, GSEA
library(clusterProfiler)

OG_glycoprotein_List_HepG2 <- OG_glycoprotein_Top_tb_HepG2$logFC
names(OG_glycoprotein_List_HepG2) <- OG_glycoprotein_Top_tb_HepG2$UniprotID
OG_glycoprotein_List_HepG2 <- sort(OG_glycoprotein_List_HepG2, decreasing = TRUE)

protein_complex_GSEA_HepG2 <- GSEA(OG_glycoprotein_List_HepG2, TERM2GENE = CORUM_database)

write_csv(
  protein_complex_GSEA_HepG2@result,
  file = 'data_source/CORUM/protein_complex_GSEA_HepG2.csv'
)

# GSEA plot
figure5B <- gseaplot(
  protein_complex_GSEA_HepG2,
  geneSetID = 1, color.line = Color_3, by = 'runningScore'
) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black")
  )

ggsave(
  filename = 'figures/figure5/figure5B.eps',
  plot = figure5B,
  height = 2, width = 2, units = 'in'
)

### figure5C, mitochondrial central dogma, OG glycoprotein HepG2
library(readxl)

# import Human.MitoCarta3.xlsx
# https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways
HumanMitoCarta_3 <- read_xlsx(
  "data_source/Mito_OG/Human.MitoCarta3.xlsx",
  sheet = 'A Human MitoCarta3.0',
  col_names = TRUE,
  .name_repair = "universal"
) 

Mito_pathway <- HumanMitoCarta_3 |> 
  select(Symbol, UniprotID = UniProt, MitoCarta3.0_MitoPathways) |> 
  separate_rows(MitoCarta3.0_MitoPathways, sep = ' \\| ') |> 
  distinct()

# extract OG glycoprotein information from HEK293T, HepG2 and Jurkat
# HEK293T
OG_glycoprotein_mito_pathway_HEK293T <- Mito_pathway |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HEK293T")

# HepG2
OG_glycoprotein_mito_pathway_HepG2 <- Mito_pathway |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HepG2")

# Jurkat
OG_glycoprotein_mito_pathway_Jurkat <- Mito_pathway |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "Jurkat")

# combine result from HEK293T, HepG2 and Jurkat
OG_glycoprotein_mito_pathway_combined <- bind_rows(
  OG_glycoprotein_mito_pathway_HEK293T,
  OG_glycoprotein_mito_pathway_HepG2,
  OG_glycoprotein_mito_pathway_Jurkat
)

# wilcoxon test
library(rstatix)
library(ggpubr)

OG_glycoprotein_central_dogma_wilcoxon_test <- OG_glycoprotein_mito_pathway_combined |> 
  filter(str_detect(MitoCarta3.0_MitoPathways, "Mitochondrial central dogma")) |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

# violin boxplot
figure5C <- OG_glycoprotein_mito_pathway_combined |> 
  ggplot() +
  geom_jitter(
    aes(
      x = cell, y = logFC, color = cell
    ), 
    position = position_jitter(width = 0.2)
  ) +
  geom_boxplot(
    aes(
      x = cell, y = logFC
    ),
    color = "black", fill = "transparent", outliers = FALSE
  ) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_color_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  stat_pvalue_manual(data = OG_glycoprotein_central_dogma_wilcoxon_test, label = "p.signif", 
                     tip.length = 0, size = 7,
                     hide.ns = "p", y.position = c(1.7, 2.2, 1.8)) +
  coord_cartesian(ylim = c(-1.8, 2.4)) +
  ggtitle("Mito. Central Dogma") +
  theme_bw() +
  theme(
    title = element_text(color = "black", size = 9, lineheight = 0.15),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(
  filename = 'figures/figure5/figure5C.eps',
  plot = figure5C,
  height = 2.5, width = 2, units = 'in'
)

### figure5D, mitochondrial disease, OG glycoprotein HepG2
# import MitoCop.xlsx
# https://www.sciencedirect.com/science/article/pii/S1550413121005295
# Supplemental information, Table S1
library(readxl)

MitoCop_mito <- read_xlsx(
  "data_source/Mito_OG/MitoCoP.xlsx",
  sheet = '(B) MitoCoP (1,134 genes)',
  skip = 1,
  col_names = TRUE,
  .name_repair = "universal"
)

MitoCop_mito_disease <- MitoCop_mito |> 
  select(UniprotID = Simplified.protein.IDs, MitoCoP.disease.gene, 
         Central.nervous.system:Tumors) |> 
  filter(MitoCoP.disease.gene == 1) |> 
  separate_rows(UniprotID, sep = ";") |> 
  pivot_longer(Central.nervous.system:Tumors, 
               names_to = 'Function', values_to = 'value') |> 
  filter(!is.na(value)) |> 
  distinct()

# extract OG glycoprotein information from HEK293T, HepG2 and Jurkat
# HEK293T
OG_glycoprotein_mito_3_HEK293T <- MitoCop_mito_disease |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HEK293T")

# HepG2
OG_glycoprotein_mito_3_HepG2 <- MitoCop_mito_disease |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HepG2")

# Jurkat
OG_glycoprotein_mito_3_Jurkat <- MitoCop_mito_disease |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "Jurkat")

# combine result from HEK293T, HepG2 and Jurkat
OG_glycoprotein_mitocop_mito_3_combined <- bind_rows(
  OG_glycoprotein_mito_3_HEK293T,
  OG_glycoprotein_mito_3_HepG2,
  OG_glycoprotein_mito_3_Jurkat
)

# wilcoxon test
library(rstatix)
library(ggpubr)

OG_glycoprotein_mitocop_disease_wilcoxon_test <- OG_glycoprotein_mitocop_mito_3_combined |> 
  filter(Function == "Metabolism" | Function == "Liver") |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

# violin boxplot
figure5D <- OG_glycoprotein_mitocop_mito_3_combined |> 
  ggplot(
    aes(x = cell, y = logFC)
  ) +
  geom_jitter(
    aes(color = cell), position = position_jitter(width = 0.2)
  ) +
  geom_boxplot(
    color = "black", fill = "transparent", outliers = FALSE
  ) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_color_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  stat_pvalue_manual(data = OG_glycoprotein_mitocop_disease_wilcoxon_test, label = "p.signif", 
                     tip.length = 0, size = 7,
                     hide.ns = "p", y.position = c(1.4, 1.7, 2.0)) +
  coord_cartesian(ylim = c(-1.8, 2.2)) +
  ggtitle("Liver & Metabolism") +
  theme_bw() +
  theme(
    title = element_text(color = "black", size = 9, lineheight = 0.15),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(
  filename = 'figures/figure5/figure5D.eps',
  plot = figure5D,
  height = 2.5, width = 2, units = 'in'
)

