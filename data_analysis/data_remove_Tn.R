# import packages
library(tidyverse)

# remove protein with Tn antigen
# OG HepG2
OG_psm_HepG2_noTn <- OG_psm_raw_HepG2 |> 
  filter(! UniprotID %in% Tn_antigen_database$Uniprot.Entry)

write_csv(
  OG_psm_HepG2_noTn,
  file = 'data_source/OGlcNAc_noTn/OG_psm_HepG2_noTn.csv'
)

# OG HEK293T
OG_psm_HEK293T_noTn <- OG_psm_raw_HEK293T |> 
  filter(! UniprotID %in% Tn_antigen_database$Uniprot.Entry)

write_csv(
  OG_psm_HEK293T_noTn,
  file = 'data_source/OGlcNAc_noTn/OG_psm_HEK293T_noTn.csv'
)

# OG Jurkat
OG_psm_Jurkat_noTn <- OG_psm_raw_Jurkat |> 
  filter(! UniprotID %in% Tn_antigen_database$Uniprot.Entry)

write_csv(
  OG_psm_Jurkat_noTn,
  file = 'data_source/OGlcNAc_noTn/OG_psm_Jurkat_noTn.csv'
)
