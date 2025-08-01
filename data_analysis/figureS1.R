# import packages
library(tidyverse)
library(psych)

# figure S1A, reproducibility for OG glycoprotein HEK293T
OG_glycoprotein_repro_HEK293T <- OG_glycoprotein_noTn_sl_tmm_HEK293T |> 
  mutate(
    Rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    Rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    Rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(starts_with("rep"))

pairs.panels(OG_glycoprotein_repro_HEK293T[1:3], lm = TRUE, main = "HEK293T")

# figure S1B, reproducibility for OG glycoprotein HepG2
OG_glycoprotein_repro_HepG2 <- OG_glycoprotein_noTn_sl_tmm_HepG2 |> 
  mutate(
    Rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    Rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    Rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(starts_with("Rep"))

pairs.panels(OG_glycoprotein_repro_HepG2[1:3], lm = TRUE, main = "HepG2")

# figure S1C, reproducibility for OG glycoprotein Jurkat
OG_glycoprotein_repro_Jurkat <- OG_glycoprotein_noTn_sl_tmm_Jurkat |> 
  mutate(
    Rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    Rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    Rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(starts_with("Rep"))

pairs.panels(OG_glycoprotein_repro_Jurkat[1:3], lm = TRUE, main = "Jurkat")
