# import packages
library(tidyverse)
library(corrplot)

# common OG glycoprotein
OG_glycoprotein_common <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter((UniprotID %in% OG_glycoprotein_Top_tb_HEK293T$UniprotID) & (UniprotID %in% OG_glycoprotein_Top_tb_Jurkat$UniprotID)) |> 
  pull(UniprotID)

# extract information for the common OG glycoprotein in HEK293T, HepG2 and Jurkat
OG_glycoprotein_common_ratio <- OG_glycoprotein_noTn_sl_tmm_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_common) |> 
  mutate(
    ratio_1_HepG2 = Tuni_1_sl_tmm/Ctrl_4_sl_tmm,
    ratio_2_HepG2 = Tuni_2_sl_tmm/Ctrl_5_sl_tmm,
    ratio_3_HepG2 = Tuni_3_sl_tmm/Ctrl_6_sl_tmm
  ) |> 
  select(UniprotID, ends_with("HepG2")) |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  mutate(
    ratio_1_HEK293T = Tuni_1_sl_tmm/Ctrl_4_sl_tmm,
    ratio_2_HEK293T = Tuni_2_sl_tmm/Ctrl_5_sl_tmm,
    ratio_3_HEK293T = Tuni_3_sl_tmm/Ctrl_6_sl_tmm
  ) |> 
  select(UniprotID, starts_with("ratio")) |> 
  left_join(OG_glycoprotein_noTn_sl_tmm_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  mutate(
    ratio_1_Jurkat = Tuni_1_sl_tmm/Ctrl_4_sl_tmm,
    ratio_2_Jurkat = Tuni_2_sl_tmm/Ctrl_5_sl_tmm,
    ratio_3_Jurkat = Tuni_3_sl_tmm/Ctrl_6_sl_tmm
  ) |> 
  select(UniprotID, starts_with("ratio"))

# correlation matrix
OG_glycorprotein_common_cor_matrix <- cor(OG_glycoprotein_common_ratio |> select(!UniprotID), method = "pearson")

rownames(OG_glycorprotein_common_cor_matrix) <- c(
  "HepG2 rep1", 
  "HepG2 rep2",
  "HepG2 rep3",
  "HEK293T rep1",
  "HEK293T rep2",
  "HEK293T rep3",
  "Jurkat rep1",
  "Jurkat rep2",
  "Jurkat rep3")

colnames(OG_glycorprotein_common_cor_matrix) <- c(
  "HepG2 rep1", 
  "HepG2 rep2",
  "HepG2 rep3",
  "HEK293T rep1",
  "HEK293T rep2",
  "HEK293T rep3",
  "Jurkat rep1",
  "Jurkat rep2",
  "Jurkat rep3")

# correlation matrix plot
corrplot(OG_glycorprotein_common_cor_matrix, type = "lower", method = "square",
         tl.cex = 1.2, tl.col = "black", tl.srt = 90,
         col = COL1('YlGn', n = 10),
         col.lim = c(0, 1),
         cl.cex = 0.8,
         is.corr = FALSE)
