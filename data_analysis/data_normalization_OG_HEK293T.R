# import packages
library(tidyverse)
library(limma)
library(edgeR)

## glycopeptide level normalization
# check the column total
OG_glycopeptide_HEK293T_noTn |> select(Tuni_1:Ctrl_6) |> colSums()

# sample loading normalization
target_mean_glycopeptide_HEK293T <- mean(colSums(OG_glycopeptide_HEK293T_noTn |> select(Tuni_1:Ctrl_6)))

norm_facs_glycopeptide_HEK293T <- target_mean_glycopeptide_HEK293T/colSums(OG_glycopeptide_HEK293T_noTn |> select(Tuni_1:Ctrl_6))

OG_glycopeptide_Tuni1_Ctrl6_sl_HEK293T <- sweep(OG_glycopeptide_HEK293T_noTn |> select(Tuni_1:Ctrl_6), 2, norm_facs_glycopeptide_HEK293T, FUN = "*")

OG_glycopeptide_Tuni1_Ctrl6_sl_tb_HEK293T <- tibble(OG_glycopeptide_Tuni1_Ctrl6_sl_HEK293T)
colnames(OG_glycopeptide_Tuni1_Ctrl6_sl_tb_HEK293T) <- c("Tuni_1_sl", "Tuni_2_sl", "Tuni_3_sl", "Ctrl_4_sl", "Ctrl_5_sl", "Ctrl_6_sl")

OG_glycopeptide_HEK293T_noTn_sl <- bind_cols(OG_glycopeptide_HEK293T_noTn, OG_glycopeptide_Tuni1_Ctrl6_sl_tb_HEK293T)

# check the sl column total
OG_glycopeptide_HEK293T_noTn_sl |> select(Tuni_1_sl:Ctrl_6_sl) |> colSums()

# TMM normalization
sl_tmm_OG_glycopeptide_HEK293T <- calcNormFactors(OG_glycopeptide_HEK293T_noTn_sl |> select(Tuni_1_sl:Ctrl_6_sl))

OG_glycopeptide_Tuni1sl_Ctrl6sl_tmm_HEK293T <- sweep(OG_glycopeptide_HEK293T_noTn_sl |> select(Tuni_1_sl:Ctrl_6_sl), 2, sl_tmm_OG_glycopeptide_HEK293T, FUN = "/")

OG_glycopeptide_Tuni1sl_Ctrl6sl_tmm_tb_HEK293T <- as_tibble(OG_glycopeptide_Tuni1sl_Ctrl6sl_tmm_HEK293T)
colnames(OG_glycopeptide_Tuni1sl_Ctrl6sl_tmm_tb_HEK293T) <- c("Tuni_1_sl_tmm", "Tuni_2_sl_tmm", "Tuni_3_sl_tmm", "Ctrl_4_sl_tmm", "Ctrl_5_sl_tmm", "Ctrl_6_sl_tmm")

OG_glycopeptide_noTn_sl_tmm_HEK293T <- bind_cols(OG_glycopeptide_HEK293T_noTn_sl, OG_glycopeptide_Tuni1sl_Ctrl6sl_tmm_tb_HEK293T)

# export glycopeptide level normalization result
write_csv(
  OG_glycopeptide_noTn_sl_tmm_HEK293T,
  file = 'data_source/OG_normalization/OG_glycopeptide_noTn_sl_tmm_HEK293T.csv'
)

## glycoprotein level normalization
# check the column total
OG_glycoprotein_HEK293T_noTn |> select(Tuni_1:Ctrl_6) |> colSums()

# sample loading normalization
target_mean_glycoprotein_HEK293T <- mean(colSums(OG_glycoprotein_HEK293T_noTn |> select(Tuni_1:Ctrl_6)))

norm_facs_glycoprotein_HEK293T <- target_mean_glycoprotein_HEK293T/colSums(OG_glycoprotein_HEK293T_noTn |> select(Tuni_1:Ctrl_6))

OG_glycoprotein_Tuni1_Ctrl6_sl_HEK293T <- sweep(OG_glycoprotein_HEK293T_noTn |> select(Tuni_1:Ctrl_6), 2, norm_facs_glycoprotein_HEK293T, FUN = "*")

OG_glycoprotein_Tuni1_Ctrl6_sl_tb_HEK293T <- tibble(OG_glycoprotein_Tuni1_Ctrl6_sl_HEK293T)
colnames(OG_glycoprotein_Tuni1_Ctrl6_sl_tb_HEK293T) <- c("Tuni_1_sl", "Tuni_2_sl", "Tuni_3_sl", "Ctrl_4_sl", "Ctrl_5_sl", "Ctrl_6_sl")

OG_glycoprotein_HEK293T_noTn_sl <- bind_cols(OG_glycoprotein_HEK293T_noTn, OG_glycoprotein_Tuni1_Ctrl6_sl_tb_HEK293T)

# check the sl column total
OG_glycoprotein_HEK293T_noTn_sl |> select(Tuni_1_sl:Ctrl_6_sl) |> colSums()

# TMM normalization
sl_tmm_OG_glycoprotein_HEK293T <- calcNormFactors(OG_glycoprotein_HEK293T_noTn_sl |> select(Tuni_1_sl:Ctrl_6_sl))

OG_glycoprotein_Tuni1sl_Ctrl6sl_tmm_HEK293T <- sweep(OG_glycoprotein_HEK293T_noTn_sl |> select(Tuni_1_sl:Ctrl_6_sl), 2, sl_tmm_OG_glycoprotein_HEK293T, FUN = "/")

OG_glycoprotein_Tuni1sl_Ctrl6sl_tmm_tb_HEK293T <- as_tibble(OG_glycoprotein_Tuni1sl_Ctrl6sl_tmm_HEK293T)
colnames(OG_glycoprotein_Tuni1sl_Ctrl6sl_tmm_tb_HEK293T) <- c("Tuni_1_sl_tmm", "Tuni_2_sl_tmm", "Tuni_3_sl_tmm", "Ctrl_4_sl_tmm", "Ctrl_5_sl_tmm", "Ctrl_6_sl_tmm")

OG_glycoprotein_noTn_sl_tmm_HEK293T <- bind_cols(OG_glycoprotein_HEK293T_noTn_sl, OG_glycoprotein_Tuni1sl_Ctrl6sl_tmm_tb_HEK293T)

# export glycoprotein level normalization result
write_csv(
  OG_glycoprotein_noTn_sl_tmm_HEK293T,
  file = 'data_source/OG_normalization/OG_glycoprotein_noTn_sl_tmm_HEK293T.csv'
)

