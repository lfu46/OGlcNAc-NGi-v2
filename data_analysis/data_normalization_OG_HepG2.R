# import packages
library(tidyverse)
library(limma)
library(edgeR)

# check the column total
OG_psm_HepG2_noTn |> select(Tuni_1:Ctrl_6) |> colSums()

# sample loading normalization
target_mean_glycopeptide_HepG2 <- mean(colSums(OG_glycopeptide_raw_HepG2 |> select(Tuni_1:Ctrl_6)))

norm_facs_glycopeptide_HepG2 <- target_mean_glycopeptide_HepG2/colSums(OG_glycopeptide_raw_HepG2 |> select(Tuni_1:Ctrl_6))

OG_glycopeptide_raw_Tuni1_Ctrl6_sl_HepG2 <- sweep(OG_glycopeptide_raw_HepG2 |> select(Tuni_1:Ctrl_6), 2, norm_facs_glycopeptide_HepG2, FUN = "*")

OG_glycopeptide_raw_Tuni1_Ctrl6_sl_tb_HepG2 <- tibble(OG_glycopeptide_raw_Tuni1_Ctrl6_sl_HepG2)
colnames(OG_glycopeptide_raw_Tuni1_Ctrl6_sl_tb_HepG2) <- c("Tuni_1_sl", "Tuni_2_sl", "Tuni_3_sl", "Ctrl_4_sl", "Ctrl_5_sl", "Ctrl_6_sl")

OG_glycopeptide_raw_sl_HepG2 <- bind_cols(OG_glycopeptide_raw_HepG2, OG_glycopeptide_raw_Tuni1_Ctrl6_sl_tb_HepG2)
