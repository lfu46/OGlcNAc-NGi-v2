# import packages
library(tidyverse)
library(circlize)
library(scales)
library(ComplexHeatmap)
library(gridBase)

# figure 4A, circular heatmap
# extract OG glycoprotein, up, down and median information from HepG2, HEK293T and Jurkat
# HepG2
OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "up", cell = "HepG2") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "down", cell = "HepG2") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

OG_glycoprotein_median_sl_tmm_logFC_adjp_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2$UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "median", cell = "HepG2") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

# HEK293T
OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "up", cell = "HEK293T") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "down", cell = "HEK293T") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

OG_glycoprotein_median_sl_tmm_logFC_adjp_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T$UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "median", cell = "HEK293T") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

# Jurkat
OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "up", cell = "Jurkat") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "down", cell = "Jurkat") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

OG_glycoprotein_median_sl_tmm_logFC_adjp_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat$UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "median", cell = "Jurkat") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

# combine differential analysis result for HepG2, HEK293T and Jurkat cell
OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total <- bind_rows(
  OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2,
  OG_glycoprotein_median_sl_tmm_logFC_adjp_HepG2,
  OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2,
  OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T,
  OG_glycoprotein_median_sl_tmm_logFC_adjp_HEK293T,
  OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T,
  OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat,
  OG_glycoprotein_median_sl_tmm_logFC_adjp_Jurkat,
  OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat
) |> 
  rowwise() |> 
  mutate(
    max_value = max(across(ends_with("sl_tmm"))),
    min_value = min(across(ends_with("sl_tmm")))
  ) |> 
  ungroup() |> 
  mutate(
    scaled_Tuni_1_sl_tmm = (Tuni_1_sl_tmm - min_value) * 2/(max_value - min_value) - 1,
    scaled_Tuni_2_sl_tmm = (Tuni_2_sl_tmm - min_value) * 2/(max_value - min_value) - 1,
    scaled_Tuni_3_sl_tmm = (Tuni_3_sl_tmm - min_value) * 2/(max_value - min_value) - 1,
    scaled_Ctrl_4_sl_tmm = (Ctrl_4_sl_tmm - min_value) * 2/(max_value - min_value) - 1,
    scaled_Ctrl_5_sl_tmm = (Ctrl_5_sl_tmm - min_value) * 2/(max_value - min_value) - 1,
    scaled_Ctrl_6_sl_tmm = (Ctrl_6_sl_tmm - min_value) * 2/(max_value - min_value) - 1
  ) |> 
  mutate(
    cell_line = ifelse(cell == "HEK293T", 1, ifelse(cell == "HepG2", 2, 3))
  )

write_csv(
  OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total,
  file = 'data_source/circular_heatmap/OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total.csv'
)

# import combined differential analysis data from HepG2, HEK923T and Jurkat cells
library(tidyverse)

OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total <- read_csv(
  'data_source/circular_heatmap/OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total.csv'
)

## key parameters for circular heatmap
# normalized TMT abundance
mat_sl_tmm <- data.matrix(OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total |> select(scaled_Tuni_1_sl_tmm:scaled_Ctrl_6_sl_tmm))

# cell type
cell <- OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total$cell

# cell type index
cell_line <- OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total$cell_line

# adjusted p value
adjpVal <- OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total$adj.P.Val

# regulated: up, down or median
category <- OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total$category

## data for circular heatmap intralink, overlap OG glycoprotein from up or down category in different cell
# from down Jurkat to up HEK293T
downJurkat_upHEK293T <- semi_join(
  OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat, OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T,
  by = "UniprotID"
)

indices_from_downJurkat_to_upHEK293T_downJurkat <- which(OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat$UniprotID %in% donwJurkat_upHEK293T$UniprotID)
indices_from_downJurkat_to_upHEK293T_upHEK293T <- which(OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T$UniprotID %in% donwJurkat_upHEK293T$UniprotID)

link_from_downJurkat_to_upHEK293T <- tribble(
  ~ group1, ~ group1_index, ~ group2, ~ group2_index,
  "Jurkat", indices_from_downJurkat_to_upHEK293T_downJurkat, "HKE293T", indices_from_downJurkat_to_upHEK293T_upHEK293T
) |> unnest(cols = c(group1_index, group2_index))

# from down HEK293T to up Jurkat
downHEK293T_upJurkat <- semi_join(
  OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T, OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat,
  by = "UniprotID"
)

indices_from_downHEK293T_to_upJurkat_downHEK293T <- which(OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T$UniprotID %in% downHEK293T_upJurkat$UniprotID)
indices_from_downHEK293T_to_upJurkat_upJurkat <- which(OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat$UniprotID %in% downHEK293T_upJurkat$UniprotID)

link_from_downHEK293T_to_upJurkat <- tribble(
  ~ group1, ~ group1_index, ~ group2, ~ group2_index,
  "HEK293T", indices_from_downHEK293T_to_upJurkat_downHEK293T, "Jurkat", indices_from_downHEK293T_to_upJurkat_upJurkat
) |> unnest(cols = c(group1_index, group2_index))

# from down HepG2 to up HEK293T
downHepG2_upHEK293T <- semi_join(
  OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2, OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T,
  by = "UniprotID"
)

indices_from_downHepG2_to_upHEK293T_downHepG2 <- which(OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2$UniprotID %in% downHepG2_upHEK293T$UniprotID)
indices_from_downHepG2_to_upHEK293T_upHEK293T <- which(OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T$UniprotID %in% downHepG2_upHEK293T$UniprotID)

link_from_downHepG2_to_upHEK293T <- tribble(
  ~ group1, ~ group1_index, ~ group2, ~ group2_index,
  "HepG2", indices_from_downHepG2_to_upHEK293T_downHepG2, "HEK293T", indices_from_downHepG2_to_upHEK293T_upHEK293T
) |> unnest(cols = c(group1_index, group2_index))

# from down HEK293T to up HepG2
downHEK293T_upHepG2 <- semi_join(
  OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T, OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2,
  by = "UniprotID"
)

indices_from_downHEK293T_to_upHepG2_downHEK293T <- which(OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T$UniprotID %in% downHEK293T_upHepG2$UniprotID)
indices_from_donwHEK293T_to_upHepG2_upHepG2 <- which(OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2$UniprotID %in% downHEK293T_upHepG2$UniprotID)

link_from_downHEK293T_to_upHepG2 <- tribble(
  ~ group1, ~ group1_index, ~ group2, ~ group2_index,
  "HEK293T", indices_from_downHEK293T_to_upHepG2_downHEK293T, "HepG2", indices_from_donwHEK293T_to_upHepG2_upHepG2
) |> unnest(cols = c(group1_index, group2_index))

# from down HepG2 to up Jurkat
downHepG2_upJurkat <- semi_join(
  OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2, OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat,
  by = "UniprotID"
)

indices_from_downHepG2_to_upJurkat_downHepG2 <- which(OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2$UniprotID %in% downHepG2_upJurkat$UniprotID)
indices_from_downHepG2_to_upJurkat_upJurkat <- which(OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat$UniprotID %in% downHepG2_upJurkat$UniprotID)

link_from_downHepG2_to_upJurkat <- tribble(
  ~ group1, ~ group1_index, ~ group2, ~ group2_index,
  "HepG2", indices_from_downHepG2_to_upJurkat_downHepG2, "Jurkat", indices_from_downHepG2_to_upJurkat_upJurkat
) |> unnest(cols = c(group1_index, group2_index))

# from down Jurkat to up HepG2
downJurkat_upHepG2 <- semi_join(
  OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat, OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2,
  by = "UniprotID"
)

indices_from_downJurkat_to_upHepG2_downJurkat <- which(OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat$UniprotID %in% downJurkat_upHepG2$UniprotID)
indices_from_downJurkat_to_upHepG2_upHepG2 <- which(OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2$UniprotID %in% downJurkat_upHepG2$UniprotID)

link_from_downJurkat_to_upHepG2 <- tribble(
  ~ group1, ~ group1_index, ~ group2, ~ group2_index,
  "Jurkat", indices_from_downJurkat_to_upHepG2_downJurkat, "HepG2", indices_from_downJurkat_to_upHepG2_upHepG2
) |> unnest(cols = c(group1_index, group2_index))

# circular heatmap
circlize_plot = function() {
  circos.heatmap(category, split = cell_line, col = col_category, track.height = 0.03)
  circos.heatmap(mat_sl_tmm, col = col_mat, track.height = 0.1)
  circos.heatmap(cell, col = col_cell, track.height = 0.03)
  circos.heatmap(adjpVal, col = col_adjpval_2, track.height = 0.03)
  
  # heatmap link 
  # from down jurkat to up HEK293T
  for(i in seq_len(nrow(link_from_downJurkat_to_upHEK293T))) {
    group1 = 3
    group2 = 1
    
    x1 = link_from_downJurkat_to_upHEK293T$group1_index[i] + 790
    x2 = link_from_downJurkat_to_upHEK293T$group2_index[i]
    
    circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_4, lwd = 3)
  }
  
  # from down HEK293T to up Jurkat
  for(i in seq_len(nrow(link_from_downHEK293T_to_upJurkat))) {
    group1 = 1
    group2 = 3
    
    x1 = link_from_downHEK293T_to_upJurkat$group1_index[i] + 1075
    x2 = link_from_downHEK293T_to_upJurkat$group2_index[i]
    
    circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_2, lwd = 3)
  }
  
  # from down HepG2 to up HEK293T
  for(i in seq_len(nrow(link_from_downHepG2_to_upHEK293T))) {
    group1 = 2
    group2 = 1
    
    x1 = link_from_downHepG2_to_upHEK293T$group1_index[i] + 327
    x2 = link_from_downHepG2_to_upHEK293T$group2_index[i]
    
    circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_3, lwd =3)
  }
  
  # from down HEK293T to up HepG2
  for(i in seq_len(nrow(link_from_downHEK293T_to_upHepG2))) {
    group1 = 1
    group2 = 2
    
    x1 = link_from_downHEK293T_to_upHepG2$group1_index[i] + 1075
    x2 = link_from_downHEK293T_to_upHepG2$group2_index[i]
    
    circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_2, lwd = 3)
  }
  
  # from down HepG2 to up Jurkat
  for(i in seq_len(nrow(link_from_downHepG2_to_upJurkat))) {
    group1 = 2
    group2 = 3
    
    x1 = link_from_downHepG2_to_upJurkat$group1_index[i] + 327
    x2 = link_from_downHepG2_to_upJurkat$group2_index[i]
    
    circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_3, lwd = 3)
  }
  
  # from down Jurkat to up HepG2
  for(i in seq_len(nrow(link_from_downJurkat_to_upHepG2))) {
    group1 = 3
    group2 = 2
    
    x1 = link_from_downJurkat_to_upHepG2$group1_index[i] + 790
    x2 = link_from_downJurkat_to_upHepG2$group2_index[i]
    
    circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_4, lwd =3)
  }
  
  circos.clear()
}

# circular heatmap color
col_category <- c("up" = Color_5, "median" = "gray", "down" = Color_6)

col_mat <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

col_cell = c("HEK293T" = Color_2, "HepG2" = Color_3, "Jurkat" = Color_4)

col_adjpval_2 = colorRamp2(c(0, 0.001, 0.01, 0.01+1e-6, 0.05, 1), 
                           c(Color_1, alpha(Color_1, 0.8), alpha(Color_1, 0.6), 
                             alpha(Color_1, 0.6), alpha(Color_1, 0.4), "white"))

# generate legend for circular heatmap
lgd_category <- Legend(title = "Category", at = names(col_category), legend_gp = gpar(fill = col_category))
lgd_mat <- Legend(title = "Z-score", col_fun = col_mat)
lgd_cell <- Legend(title = "Cell", at = names(col_cell), legend_gp = gpar(fill = col_cell))
lgd_adjpval <- Legend(title = "adj.P.Value", col_fun = col_adjpval_2, at = c(0, 0.001, 0.01, 0.05, 1),
                      break_dist = c(1, 1, 3, 3))

# combine circular heatmap and legends
plot.new()
circle_size = unit(1, "snpc")

pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circlize_plot()
upViewport()

h = dev.size()[2]
lgd_list <- packLegend(lgd_mat, lgd_cell, lgd_category, lgd_adjpval, 
                       max_height = unit(0.9*h, "inch"))
draw(lgd_list, x = circle_size, just = "left")

### figure 4B, OG glycoprotein differentailly regulated, HEK293T vs. Jurkat
# extract OG glycoprotein HEK293T vs. Jurkat
OG_glycoprotein_diff_HEK293TvsJurkat <- bind_rows(downHEK293T_upJurkat, downJurkat_upHEK293T) |> 
  select(UniprotID)

write_csv(
  OG_glycoprotein_diff_HEK293TvsJurkat,
  file = 'data_source/HEK293TvsJurkat/OG_glycoprotein_diff_HEK293TvsJurkat.csv'
)

# import OG glycoprotein HEK923T vs. Jurkat
library(tidyverse)

OG_glycoprotein_diff_HEK293TvsJurkat <- read_csv(
  'data_source/HEK293TvsJurkat/OG_glycoprotein_diff_HEK293TvsJurkat.csv'
)

OG_glycoprotein_diff_HEK293TvsJurkat_logFC <- OG_glycoprotein_diff_HEK293TvsJurkat |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  select(UniprotID, logFC_HEK293T = logFC) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  select(UniprotID, logFC_HEK293T, logFC_Jurkat = logFC)

OG_glycoprotein_diff_HEK293TvsJurkat_sl_tmm_logFC <- OG_glycoprotein_diff_HEK293TvsJurkat_logFC |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_HEK293T, by = "UniprotID") |> 
  select(UniprotID, logFC_HEK293T, logFC_Jurkat, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_HEK293T = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    log2_ratio2_HEK293T = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    log2_ratio3_HEK293T = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, logFC_HEK293T, log2_ratio1_HEK293T:log2_ratio3_HEK293T, logFC_Jurkat) |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_Jurkat, by = "UniprotID") |> 
  mutate(
    log2_ratio1_Jurkat = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    log2_ratio2_Jurkat = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    log2_ratio3_Jurkat = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, logFC_HEK293T, log2_ratio1_HEK293T:log2_ratio3_HEK293T, logFC_Jurkat, log2_ratio1_Jurkat:log2_ratio3_Jurkat)

# complex heatmap
library(ComplexHeatmap)

mat <- data.matrix(OG_glycoprotein_diff_HEK293TvsJurkat_sl_tmm_logFC |> select(starts_with("log2_ratio")))
rownames(mat) <- OG_glycoprotein_diff_HEK293TvsJurkat_sl_tmm_logFC$UniprotID
colnames(mat) <- c("HEK293T 1", "HEK293T 2", "HEK293T 3", "Jurkat 1", "Jurkat 2", "Jurkat 3")

anno_matrix <- data.matrix(OG_glycoprotein_diff_HEK293TvsJurkat_sl_tmm_logFC |> select(starts_with("log2_ratio")))

ha1 = rowAnnotation("log2(Tuni/Ctrl)" = anno_points(anno_matrix,
                                                    pch = c(1, 1, 1, 2, 2, 2), 
                                                    gp = gpar(col = c(Color_2, Color_2, Color_2, Color_4, Color_4, Color_4), lwd = 2),
                                                    ylim = c(-3, 3),
                                                    size = unit(4, "mm"),
                                                    width = unit(4, "cm"),
                                                    axis_param = list(
                                                      at = c(-3, 0, 3)
                                                    )
))

Heatmap(mat, name = "log2(Tuni/Ctrl)", show_row_names = TRUE, show_column_names = TRUE,
              right_annotation = ha1)

### figure 4C, OG up HEK293T, gene ontology
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

# upregulated OG glycoprotein, HEK293T
OG_glycoprotein_up_HEK293T_list <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID) |> 
  pull()

# gene ontology analysis
OG_glycoprotein_up_HEK293T_GO <- enrichGO(
  gene = OG_glycoprotein_up_HEK293T_list,
  OrgDb = org.Hs.eg.db,
  universe = OG_glycoprotein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_up_HEK293T_GO@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_up_HEK293T_GO.csv'
)

# import gene ontology, OG up HEK293T
OG_glycoprotein_up_HEK293T_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_up_HEK293T_GO.csv'
)

# barplot
figure4C <- OG_glycoprotein_up_HEK293T_GO |> 
  filter(
    Description %in% c(
      'preribosome', 'rRNA processing', 'nucleotide binding', 'RNA binding', 'nucleolus'
    )
  ) |> 
  ggplot() +
  geom_bar(
    aes(x = -log10(pvalue), y = fct_reorder(Description, -log10(pvalue))), 
    stat = "identity", fill = Color_2, width = 0.7
  ) +
  labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
  geom_text(aes(label = Description, x = 0, y = Description), hjust = 0, size = 3) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "in")
  )

ggsave(
  filename = 'figures/figure4/figure4C.eps', 
  plot = figure4C, 
  height = 1.5, width = 2, units = 'in'
)

### figure 4D, OG up HepG2, gene ontology
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

# upregulated OG glycoprotein, HepG2
OG_glycoprotein_up_HepG2_list <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID) |> 
  pull()

# gene ontology analysis
OG_glycoprotein_up_HepG2_GO <- enrichGO(
  gene = OG_glycoprotein_up_HepG2_list,
  OrgDb = org.Hs.eg.db,
  universe = OG_glycoprotein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_up_HepG2_GO@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_up_HepG2_GO.csv'
)

# import gene ontology, OG up HepG2
OG_glycoprotein_up_HepG2_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_up_HepG2_GO.csv'
)

# barplot
figure4D <- OG_glycoprotein_up_HepG2_GO |> 
  filter(
    Description %in% c(
      'extracellular exosome', 'extracellular vesicle', 'mitochondrial matrix', 'RNA binding', 'GTP binding'
    )
  ) |> 
  ggplot() +
  geom_bar(
    aes(x = -log10(pvalue), y = fct_reorder(Description, -log10(pvalue))), 
    stat = "identity", fill = Color_3, width = 0.7
  ) +
  labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
  geom_text(aes(label = Description, x = 0, y = Description), hjust = 0, size = 3) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "in")
  )

ggsave(
  filename = 'figures/figure4/figure4D.eps', 
  plot = figure4D, 
  height = 1.5, width = 2, units = 'in'
)

### figure 4E, OG down Jurkat, gene ontology
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

# downregulated OG glycoprotein, Jurkat
OG_glycoprotein_down_Jurkat_list <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID) |> 
  pull()

# gene ontology analysis
OG_glycoprotein_down_Jurkat_GO <- enrichGO(
  gene = OG_glycoprotein_down_Jurkat_list,
  OrgDb = org.Hs.eg.db,
  universe = OG_glycoprotein_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OG_glycoprotein_down_Jurkat_GO@result,
  file = 'data_source/gene_ontology/OG_glycoprotein_down_Jurkat_GO.csv'
)

# import gene ontology, OG down Jurkat
OG_glycoprotein_down_Jurkat_GO <- read_csv(
  'data_source/gene_ontology/OG_glycoprotein_down_Jurkat_GO.csv'
)

# barplot
figure4E <- OG_glycoprotein_down_Jurkat_GO |> 
  filter(
    Description %in% c(
      'glycoprotein biosynthetic process', 'carbohydrate derivative metabolic process', 'Golgi membrane', 'protein glycosylation', 'condensed nuclear chromosome'
    )
  ) |> 
  ggplot() +
  geom_bar(
    aes(x = -log10(pvalue), y = fct_reorder(Description, -log10(pvalue))), 
    stat = "identity", fill = Color_4, width = 0.7
  ) +
  labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
  # geom_text(aes(label = Description, x = 0, y = Description), hjust = 0, size = 3) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "in")
  )

ggsave(
  filename = 'figures/figure4/figure4E.eps', 
  plot = figure4E, 
  height = 1.5, width = 2, units = 'in'
)
