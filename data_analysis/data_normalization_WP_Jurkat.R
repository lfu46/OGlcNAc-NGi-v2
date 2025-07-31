# import packages
library(tidyverse)
library(limma)
library(edgeR)

## protein level normalization
# check the column total
WP_protein_Jurkat |> select(Tuni_1:Ctrl_6) |> colSums()

# sample loading normalization
target_mean_protein_Jurkat <- mean(colSums(WP_protein_Jurkat |> select(Tuni_1:Ctrl_6)))

norm_facs_protein_Jurkat <- target_mean_protein_Jurkat/colSums(WP_protein_Jurkat |> select(Tuni_1:Ctrl_6))

WP_protein_Tuni1_Ctrl6_sl_Jurkat <- sweep(WP_protein_Jurkat |> select(Tuni_1:Ctrl_6), 2, norm_facs_protein_Jurkat, FUN = "*")

WP_protein_Tuni1_Ctrl6_sl_tb_Jurkat <- tibble(WP_protein_Tuni1_Ctrl6_sl_Jurkat)
colnames(WP_protein_Tuni1_Ctrl6_sl_tb_Jurkat) <- c("Tuni_1_sl", "Tuni_2_sl", "Tuni_3_sl", "Ctrl_4_sl", "Ctrl_5_sl", "Ctrl_6_sl")

WP_protein_Jurkat_sl <- bind_cols(WP_protein_Jurkat, WP_protein_Tuni1_Ctrl6_sl_tb_Jurkat)

# check the sl column total
WP_protein_Jurkat_sl |> select(Tuni_1_sl:Ctrl_6_sl) |> colSums()

# TMM normalization
sl_tmm_WP_protein_Jurkat <- calcNormFactors(WP_protein_Jurkat_sl |> select(Tuni_1_sl:Ctrl_6_sl))

WP_protein_Tuni1sl_Ctrl6sl_tmm_Jurkat <- sweep(WP_protein_Jurkat_sl |> select(Tuni_1_sl:Ctrl_6_sl), 2, sl_tmm_WP_protein_Jurkat, FUN = "/")

WP_protein_Tuni1sl_Ctrl6sl_tmm_tb_Jurkat <- as_tibble(WP_protein_Tuni1sl_Ctrl6sl_tmm_Jurkat)
colnames(WP_protein_Tuni1sl_Ctrl6sl_tmm_tb_Jurkat) <- c("Tuni_1_sl_tmm", "Tuni_2_sl_tmm", "Tuni_3_sl_tmm", "Ctrl_4_sl_tmm", "Ctrl_5_sl_tmm", "Ctrl_6_sl_tmm")

WP_protein_sl_tmm_Jurkat <- bind_cols(WP_protein_Jurkat_sl, WP_protein_Tuni1sl_Ctrl6sl_tmm_tb_Jurkat)

# export protein level normalization result
write_csv(
  WP_protein_sl_tmm_Jurkat,
  file = 'data_source/WP_normalization/WP_protein_sl_tmm_Jurkat.csv'
)

