# import packages
library(tidyverse)
library(limma)
library(edgeR)

## glycopeptide level normalization
# check the column total
OG_glycopeptide_Jurkat_noTn |> select(Tuni_1:Ctrl_6) |> colSums()

# sample loading normalization
target_mean_glycopeptide_Jurkat <- mean(colSums(OG_glycopeptide_Jurkat_noTn |> select(Tuni_1:Ctrl_6)))

norm_facs_glycopeptide_Jurkat <- target_mean_glycopeptide_Jurkat/colSums(OG_glycopeptide_Jurkat_noTn |> select(Tuni_1:Ctrl_6))

OG_glycopeptide_Tuni1_Ctrl6_sl_Jurkat <- sweep(OG_glycopeptide_Jurkat_noTn |> select(Tuni_1:Ctrl_6), 2, norm_facs_glycopeptide_Jurkat, FUN = "*")

OG_glycopeptide_Tuni1_Ctrl6_sl_tb_Jurkat <- tibble(OG_glycopeptide_Tuni1_Ctrl6_sl_Jurkat)
colnames(OG_glycopeptide_Tuni1_Ctrl6_sl_tb_Jurkat) <- c("Tuni_1_sl", "Tuni_2_sl", "Tuni_3_sl", "Ctrl_4_sl", "Ctrl_5_sl", "Ctrl_6_sl")

OG_glycopeptide_Jurkat_noTn_sl <- bind_cols(OG_glycopeptide_Jurkat_noTn, OG_glycopeptide_Tuni1_Ctrl6_sl_tb_Jurkat)

# check the sl column total
OG_glycopeptide_Jurkat_noTn_sl |> select(Tuni_1_sl:Ctrl_6_sl) |> colSums()

# TMM normalization
sl_tmm_OG_glycopeptide_Jurkat <- calcNormFactors(OG_glycopeptide_Jurkat_noTn_sl |> select(Tuni_1_sl:Ctrl_6_sl))

OG_glycopeptide_Tuni1sl_Ctrl6sl_tmm_Jurkat <- sweep(OG_glycopeptide_Jurkat_noTn_sl |> select(Tuni_1_sl:Ctrl_6_sl), 2, sl_tmm_OG_glycopeptide_Jurkat, FUN = "/")

OG_glycopeptide_Tuni1sl_Ctrl6sl_tmm_tb_Jurkat <- as_tibble(OG_glycopeptide_Tuni1sl_Ctrl6sl_tmm_Jurkat)
colnames(OG_glycopeptide_Tuni1sl_Ctrl6sl_tmm_tb_Jurkat) <- c("Tuni_1_sl_tmm", "Tuni_2_sl_tmm", "Tuni_3_sl_tmm", "Ctrl_4_sl_tmm", "Ctrl_5_sl_tmm", "Ctrl_6_sl_tmm")

OG_glycopeptide_noTn_sl_tmm_Jurkat <- bind_cols(OG_glycopeptide_Jurkat_noTn_sl, OG_glycopeptide_Tuni1sl_Ctrl6sl_tmm_tb_Jurkat)

# export glycopeptide level normalization result
write_csv(
  OG_glycopeptide_noTn_sl_tmm_Jurkat,
  file = 'data_source/OG_normalization/OG_glycopeptide_noTn_sl_tmm_Jurkat.csv'
)

## glycoprotein level normalization
# check the column total
OG_glycoprotein_Jurkat_noTn |> select(Tuni_1:Ctrl_6) |> colSums()

# sample loading normalization
target_mean_glycoprotein_Jurkat <- mean(colSums(OG_glycoprotein_Jurkat_noTn |> select(Tuni_1:Ctrl_6)))

norm_facs_glycoprotein_Jurkat <- target_mean_glycoprotein_Jurkat/colSums(OG_glycoprotein_Jurkat_noTn |> select(Tuni_1:Ctrl_6))

OG_glycoprotein_Tuni1_Ctrl6_sl_Jurkat <- sweep(OG_glycoprotein_Jurkat_noTn |> select(Tuni_1:Ctrl_6), 2, norm_facs_glycoprotein_Jurkat, FUN = "*")

OG_glycoprotein_Tuni1_Ctrl6_sl_tb_Jurkat <- tibble(OG_glycoprotein_Tuni1_Ctrl6_sl_Jurkat)
colnames(OG_glycoprotein_Tuni1_Ctrl6_sl_tb_Jurkat) <- c("Tuni_1_sl", "Tuni_2_sl", "Tuni_3_sl", "Ctrl_4_sl", "Ctrl_5_sl", "Ctrl_6_sl")

OG_glycoprotein_Jurkat_noTn_sl <- bind_cols(OG_glycoprotein_Jurkat_noTn, OG_glycoprotein_Tuni1_Ctrl6_sl_tb_Jurkat)

# check the sl column total
OG_glycoprotein_Jurkat_noTn_sl |> select(Tuni_1_sl:Ctrl_6_sl) |> colSums()

# TMM normalization
sl_tmm_OG_glycoprotein_Jurkat <- calcNormFactors(OG_glycoprotein_Jurkat_noTn_sl |> select(Tuni_1_sl:Ctrl_6_sl))

OG_glycoprotein_Tuni1sl_Ctrl6sl_tmm_Jurkat <- sweep(OG_glycoprotein_Jurkat_noTn_sl |> select(Tuni_1_sl:Ctrl_6_sl), 2, sl_tmm_OG_glycoprotein_Jurkat, FUN = "/")

OG_glycoprotein_Tuni1sl_Ctrl6sl_tmm_tb_Jurkat <- as_tibble(OG_glycoprotein_Tuni1sl_Ctrl6sl_tmm_Jurkat)
colnames(OG_glycoprotein_Tuni1sl_Ctrl6sl_tmm_tb_Jurkat) <- c("Tuni_1_sl_tmm", "Tuni_2_sl_tmm", "Tuni_3_sl_tmm", "Ctrl_4_sl_tmm", "Ctrl_5_sl_tmm", "Ctrl_6_sl_tmm")

OG_glycoprotein_noTn_sl_tmm_Jurkat <- bind_cols(OG_glycoprotein_Jurkat_noTn_sl, OG_glycoprotein_Tuni1sl_Ctrl6sl_tmm_tb_Jurkat)

# export glycoprotein level normalization result
write_csv(
  OG_glycoprotein_noTn_sl_tmm_Jurkat,
  file = 'data_source/OG_normalization/OG_glycoprotein_noTn_sl_tmm_Jurkat.csv'
)

