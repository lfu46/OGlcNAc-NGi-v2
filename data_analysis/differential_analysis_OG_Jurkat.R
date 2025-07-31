# import packages
library(tidyverse)
library(limma)

# experimental model
Experiment_Model <- model.matrix(~ 0 + factor(rep(c("case1", "case2"), each = 3), levels = c("case1", "case2")))
colnames(Experiment_Model) <- c("case1", "case2")
Contrast_matrix <- makeContrasts(case1_case2 = case1 - case2, levels = Experiment_Model)

# glycopeptide level
OG_glycopeptide_noTn_sl_tmm_log2_Jurkat <- OG_glycopeptide_noTn_sl_tmm_Jurkat |> 
  mutate(
    log2_Tuni_1_sl_tmm = log2(Tuni_1_sl_tmm),
    log2_Tuni_2_sl_tmm = log2(Tuni_2_sl_tmm),
    log2_Tuni_3_sl_tmm = log2(Tuni_3_sl_tmm),
    log2_Ctrl_4_sl_tmm = log2(Ctrl_4_sl_tmm),
    log2_Ctrl_5_sl_tmm = log2(Ctrl_5_sl_tmm),
    log2_Ctrl_6_sl_tmm = log2(Ctrl_6_sl_tmm)
  )

OG_glycopeptide_noTn_sl_tmm_log2_Data_Matrix_Jurkat <- data.matrix(OG_glycopeptide_noTn_sl_tmm_log2_Jurkat |> select(starts_with("log2")))
rownames(OG_glycopeptide_noTn_sl_tmm_log2_Data_Matrix_Jurkat) <- OG_glycopeptide_noTn_sl_tmm_log2_Jurkat$Index

OG_glycopeptide_noTn_sl_tmm_log2_Fit_Jurkat <- lmFit(OG_glycopeptide_noTn_sl_tmm_log2_Data_Matrix_Jurkat, Experiment_Model)
OG_glycopeptide_noTn_sl_tmm_log2_Fit_Contrast_Jurkat <- contrasts.fit(OG_glycopeptide_noTn_sl_tmm_log2_Fit_Jurkat, Contrast_matrix)
OG_glycopeptide_noTn_sl_tmm_log2_Fit_Contrast_Jurkat <- eBayes(OG_glycopeptide_noTn_sl_tmm_log2_Fit_Contrast_Jurkat)
OG_glycopeptide_Top_Jurkat <- topTable(OG_glycopeptide_noTn_sl_tmm_log2_Fit_Contrast_Jurkat, number = Inf, adjust.method = "BH")

Rownames_OG_glycopeptide_Top_Jurkat <- rownames(OG_glycopeptide_Top_Jurkat)
OG_glycopeptide_Top_tb_Jurkat <- as_tibble(OG_glycopeptide_Top_Jurkat)
OG_glycopeptide_Top_tb_Jurkat$Index <- Rownames_OG_glycopeptide_Top_Jurkat

write_csv(
  OG_glycopeptide_Top_tb_Jurkat,
  file = 'data_source/OG_differential_analysis/OG_glycopeptide_Top_tb_Jurkat.csv'
)

# glycoprotein level
OG_glycoprotein_noTn_sl_tmm_log2_Jurkat <- OG_glycoprotein_noTn_sl_tmm_Jurkat |> 
  mutate(
    log2_Tuni_1_sl_tmm = log2(Tuni_1_sl_tmm),
    log2_Tuni_2_sl_tmm = log2(Tuni_2_sl_tmm),
    log2_Tuni_3_sl_tmm = log2(Tuni_3_sl_tmm),
    log2_Ctrl_4_sl_tmm = log2(Ctrl_4_sl_tmm),
    log2_Ctrl_5_sl_tmm = log2(Ctrl_5_sl_tmm),
    log2_Ctrl_6_sl_tmm = log2(Ctrl_6_sl_tmm)
  )

OG_glycoprotein_noTn_sl_tmm_log2_Data_Matrix_Jurkat <- data.matrix(OG_glycoprotein_noTn_sl_tmm_log2_Jurkat |> select(starts_with("log2")))
rownames(OG_glycoprotein_noTn_sl_tmm_log2_Data_Matrix_Jurkat) <- OG_glycoprotein_noTn_sl_tmm_log2_Jurkat$UniprotID

OG_glycoprotein_noTn_sl_tmm_log2_Fit_Jurkat <- lmFit(OG_glycoprotein_noTn_sl_tmm_log2_Data_Matrix_Jurkat, Experiment_Model)
OG_glycoprotein_noTn_sl_tmm_log2_Fit_Contrast_Jurkat <- contrasts.fit(OG_glycoprotein_noTn_sl_tmm_log2_Fit_Jurkat, Contrast_matrix)
OG_glycoprotein_noTn_sl_tmm_log2_Fit_Contrast_Jurkat <- eBayes(OG_glycoprotein_noTn_sl_tmm_log2_Fit_Contrast_Jurkat)
OG_glycoprotein_Top_Jurkat <- topTable(OG_glycoprotein_noTn_sl_tmm_log2_Fit_Contrast_Jurkat, number = Inf, adjust.method = "BH")

Rownames_OG_glycoprotein_Top_Jurkat <- rownames(OG_glycoprotein_Top_Jurkat)
OG_glycoprotein_Top_tb_Jurkat <- as_tibble(OG_glycoprotein_Top_Jurkat)
OG_glycoprotein_Top_tb_Jurkat$UniprotID <- Rownames_OG_glycoprotein_Top_Jurkat

write_csv(
  OG_glycoprotein_Top_tb_Jurkat,
  file = 'data_source/OG_differential_analysis/OG_glycoprotein_Top_tb_Jurkat.csv'
)
