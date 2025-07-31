# import packages
library(tidyverse)

## glycopeptide level quantification
# OG HepG2
OG_glycopeptide_HepG2_noTn <- OG_psm_HepG2_noTn |> 
  group_by(Index, UniprotID, combined_site, Start.Position, End.Position, identified, localized) |> 
  summarize(
    Tuni_1 = sum(Sn.126),
    Tuni_2 = sum(Sn.127n),
    Tuni_3 = sum(Sn.128c),
    Ctrl_4 = sum(Sn.129n),
    Ctrl_5 = sum(Sn.130c),
    Ctrl_6 = sum(Sn.131)
  ) |> 
  ungroup()

write_csv(
  OG_glycopeptide_HepG2_noTn,
  file = 'data_source/OG_quantification/OG_glycopeptide_HepG2_noTn.csv'
)

# OG HEK293T
OG_glycopeptide_HEK293T_noTn <- OG_psm_HEK293T_noTn |> 
  group_by(Index, UniprotID, combined_site, Start.Position, End.Position, identified, localized) |> 
  summarize(
    Tuni_1 = sum(Sn.126),
    Tuni_2 = sum(Sn.127n),
    Tuni_3 = sum(Sn.128c),
    Ctrl_4 = sum(Sn.129n),
    Ctrl_5 = sum(Sn.130c),
    Ctrl_6 = sum(Sn.131)
  ) |> 
  ungroup()

write_csv(
  OG_glycopeptide_HEK293T_noTn,
  file = 'data_source/OG_quantification/OG_glycopeptide_HEK293T_noTn.csv'
)

# OG Jurkat
OG_glycopeptide_Jurkat_noTn <- OG_psm_Jurkat_noTn |> 
  group_by(Index, UniprotID, combined_site, Start.Position, End.Position, identified, localized) |> 
  summarize(
    Tuni_1 = sum(Sn.126),
    Tuni_2 = sum(Sn.127n),
    Tuni_3 = sum(Sn.128c),
    Ctrl_4 = sum(Sn.129n),
    Ctrl_5 = sum(Sn.130c),
    Ctrl_6 = sum(Sn.131)
  ) |> 
  ungroup()

write_csv(
  OG_glycopeptide_Jurkat_noTn,
  file = 'data_source/OG_quantification/OG_glycopeptide_Jurkat_noTn.csv'
)

## glycoprotein level quantification
# OG HepG2
OG_glycoprotein_HepG2_noTn <- OG_glycopeptide_HepG2_noTn |> 
  group_by(UniprotID) |> 
  summarize(
    Tuni_1 = sum(Tuni_1),
    Tuni_2 = sum(Tuni_2),
    Tuni_3 = sum(Tuni_3),
    Ctrl_4 = sum(Ctrl_4),
    Ctrl_5 = sum(Ctrl_5),
    Ctrl_6 = sum(Ctrl_6)
  )

write_csv(
  OG_glycoprotein_HepG2_noTn,
  file = 'data_source/OG_quantification/OG_glycoprotein_HepG2_noTn.csv'
)

# OG HEK293T
OG_glycoprotein_HEK293T_noTn <- OG_glycopeptide_HEK293T_noTn |> 
  group_by(UniprotID) |> 
  summarize(
    Tuni_1 = sum(Tuni_1),
    Tuni_2 = sum(Tuni_2),
    Tuni_3 = sum(Tuni_3),
    Ctrl_4 = sum(Ctrl_4),
    Ctrl_5 = sum(Ctrl_5),
    Ctrl_6 = sum(Ctrl_6)
  )

write_csv(
  OG_glycoprotein_HEK293T_noTn,
  file = 'data_source/OG_quantification/OG_glycoprotein_HEK293T_noTn.csv'
)

# OG Jurkat
OG_glycoprotein_Jurkat_noTn <- OG_glycopeptide_Jurkat_noTn |> 
  group_by(UniprotID) |> 
  summarize(
    Tuni_1 = sum(Tuni_1),
    Tuni_2 = sum(Tuni_2),
    Tuni_3 = sum(Tuni_3),
    Ctrl_4 = sum(Ctrl_4),
    Ctrl_5 = sum(Ctrl_5),
    Ctrl_6 = sum(Ctrl_6)
  )

write_csv(
  OG_glycoprotein_Jurkat_noTn,
  file = 'data_source/OG_quantification/OG_glycoprotein_Jurkat_noTn.csv'
)
