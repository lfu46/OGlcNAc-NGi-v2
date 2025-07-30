# import packages
library(tidyverse)

## protein level quantification
# HepG2
WP_protein_HepG2 <- WP_psm_raw_HepG2 |> 
  group_by(UniprotID) |> 
  summarize(
    Tuni_1 = sum(Sn.126),
    Tuni_2 = sum(Sn.127n),
    Tuni_3 = sum(Sn.128c),
    Ctrl_4 = sum(Sn.129n),
    Ctrl_5 = sum(Sn.130c),
    Ctrl_6 = sum(Sn.131)
  )

write_csv(
  WP_protein_HepG2,
  file = 'data_source/WP_quantification/WP_protein_HepG2.csv'
)

# HEK293T
WP_protein_HEK293T <- WP_psm_raw_HEK293T |> 
  group_by(UniprotID) |> 
  summarize(
    Tuni_1 = sum(Sn.126),
    Tuni_2 = sum(Sn.127n),
    Tuni_3 = sum(Sn.128c),
    Ctrl_4 = sum(Sn.129n),
    Ctrl_5 = sum(Sn.130c),
    Ctrl_6 = sum(Sn.131)
  )

write_csv(
  WP_protein_HEK293T,
  file = 'data_source/WP_quantification/WP_protein_HEK293T.csv'
)

# Jurkat
WP_protein_Jurkat <- WP_psm_raw_Jurkat |> 
  group_by(UniprotID) |> 
  summarize(
    Tuni_1 = sum(Sn.126),
    Tuni_2 = sum(Sn.127n),
    Tuni_3 = sum(Sn.128c),
    Ctrl_4 = sum(Sn.129n),
    Ctrl_5 = sum(Sn.130c),
    Ctrl_6 = sum(Sn.131)
  )

write_csv(
  WP_protein_Jurkat,
  file = 'data_source/WP_quantification/WP_protein_Jurkat.csv'
)
