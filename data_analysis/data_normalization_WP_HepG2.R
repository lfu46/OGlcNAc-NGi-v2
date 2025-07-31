# import packages
library(tidyverse)
library(limma)
library(edgeR)

## protein level normalization
# check the column total
WP_protein_HepG2 |> select(Tuni_1:Ctrl_6) |> colSums()

# sample loading normalization
target_mean_protein_HepG2 <- mean(colSums(WP_protein_HepG2 |> select(Tuni_1:Ctrl_6)))

norm_facs_protein_HepG2 <- target_mean_protein_HepG2/colSums(WP_protein_HepG2 |> select(Tuni_1:Ctrl_6))

WP_protein_Tuni1_Ctrl6_sl_HepG2 <- sweep(WP_protein_HepG2 |> select(Tuni_1:Ctrl_6), 2, norm_facs_protein_HepG2, FUN = "*")

WP_protein_Tuni1_Ctrl6_sl_tb_HepG2 <- tibble(WP_protein_Tuni1_Ctrl6_sl_HepG2)
colnames(WP_protein_Tuni1_Ctrl6_sl_tb_HepG2) <- c("Tuni_1_sl", "Tuni_2_sl", "Tuni_3_sl", "Ctrl_4_sl", "Ctrl_5_sl", "Ctrl_6_sl")

WP_protein_HepG2_sl <- bind_cols(WP_protein_HepG2, WP_protein_Tuni1_Ctrl6_sl_tb_HepG2)

# check the sl column total
WP_protein_HepG2_sl |> select(Tuni_1_sl:Ctrl_6_sl) |> colSums()

# TMM normalization
sl_tmm_WP_protein_HepG2 <- calcNormFactors(WP_protein_HepG2_sl |> select(Tuni_1_sl:Ctrl_6_sl))

WP_protein_Tuni1sl_Ctrl6sl_tmm_HepG2 <- sweep(WP_protein_HepG2_sl |> select(Tuni_1_sl:Ctrl_6_sl), 2, sl_tmm_WP_protein_HepG2, FUN = "/")

WP_protein_Tuni1sl_Ctrl6sl_tmm_tb_HepG2 <- as_tibble(WP_protein_Tuni1sl_Ctrl6sl_tmm_HepG2)
colnames(WP_protein_Tuni1sl_Ctrl6sl_tmm_tb_HepG2) <- c("Tuni_1_sl_tmm", "Tuni_2_sl_tmm", "Tuni_3_sl_tmm", "Ctrl_4_sl_tmm", "Ctrl_5_sl_tmm", "Ctrl_6_sl_tmm")

WP_protein_sl_tmm_HepG2 <- bind_cols(WP_protein_HepG2_sl, WP_protein_Tuni1sl_Ctrl6sl_tmm_tb_HepG2)

# export protein level normalization result
write_csv(
  WP_protein_sl_tmm_HepG2,
  file = 'data_source/WP_normalization/WP_protein_sl_tmm_HepG2.csv'
)

