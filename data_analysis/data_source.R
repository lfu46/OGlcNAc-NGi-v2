# import packages
library(tidyverse)

# raw file
# OG
OG_psm_raw_HepG2 <- read_csv(
  'data_source/raw_file/OG_psm_raw_HepG2.csv'
)
OG_psm_raw_HEK293T <- read_csv(
  'data_source/raw_file/OG_psm_raw_HEK293T.csv'
)
OG_psm_raw_Jurkat <- read_csv(
  'data_source/raw_file/OG_psm_raw_Jurkat.csv'
)

# WP
WP_psm_raw_HepG2 <- read_csv(
  'data_source/raw_file/WP_psm_raw_HepG2.csv'
)
WP_psm_raw_HEK293T <- read_csv(
  'data_source/raw_file/WP_psm_raw_HEK293T.csv'
)
WP_psm_raw_Jurkat <- read_csv(
  'data_source/raw_file/WP_psm_raw_Jurkat.csv'
)

# Tn antigen database
Tn_antigen_database <- read_csv(
  'data_source/Tn_antigen_database/Tn_antigen_database.csv'
)

# OG remove Tn antigen
OG_psm_HepG2_noTn <- read_csv(
  'data_source/OGlcNAc_noTn/OG_psm_HepG2_noTn.csv'
)
OG_psm_HEK293T_noTn <- read_csv(
  'data_source/OGlcNAc_noTn/OG_psm_HEK293T_noTn.csv'
)
OG_psm_Jurkat_noTn <- read_csv(
  'data_source/OGlcNAc_noTn/OG_psm_Jurkat_noTn.csv'
)

### OG quantification
## glycopeptide level
# HepG2
OG_glycopeptide_HepG2_noTn <- read_csv(
  'data_source/OG_quantification/OG_glycopeptide_HepG2_noTn.csv'
)

# HEK293T
OG_glycopeptide_HEK293T_noTn <- read_csv(
  'data_source/OG_quantification/OG_glycopeptide_HEK293T_noTn.csv'
)

# Jurkat
OG_glycopeptide_Jurkat_noTn <- read_csv(
  'data_source/OG_quantification/OG_glycopeptide_Jurkat_noTn.csv'
)

## glycoprotein level
# HepG2
OG_glycoprotein_HepG2_noTn <- read_csv(
  'data_source/OG_quantification/OG_glycoprotein_HepG2_noTn.csv'
)

# HEK293T
OG_glycoprotein_HEK293T_noTn <- read_csv(
  'data_source/OG_quantification/OG_glycoprotein_HEK293T_noTn.csv'
)

# Jurkat
OG_glycoprotein_Jurkat_noTn <- read_csv(
  'data_source/OG_quantification/OG_glycoprotein_Jurkat_noTn.csv'
)

# WP quantification
# HepG2
WP_protein_HepG2 <- read_csv(
  'data_source/WP_quantification/WP_protein_HepG2.csv'
)

# HEK293T
WP_protein_HEK293T <- read_csv(
  'data_source/WP_quantification/WP_protein_HEK293T.csv'
)

# Jurkat
WP_protein_Jurkat <- read_csv(
  'data_source/WP_quantification/WP_protein_HEK293T.csv'
)

### OG normalization
## glycopeptide level
# HepG2
OG_glycopeptide_noTn_sl_tmm_HepG2 <- read_csv(
  'data_source/OG_normalization/OG_glycopeptide_noTn_sl_tmm_HepG2.csv'
)

# HEK293T
OG_glycopeptide_noTn_sl_tmm_HEK293T <- read_csv(
  'data_source/OG_normalization/OG_glycopeptide_noTn_sl_tmm_HEK293T.csv'
)

# Jurkat
OG_glycopeptide_noTn_sl_tmm_Jurkat <- read_csv(
  'data_source/OG_normalization/OG_glycopeptide_noTn_sl_tmm_Jurkat.csv'
)

## glycoprotein level
# HepG2
OG_glycoprotein_noTn_sl_tmm_HepG2 <- read_csv(
  'data_source/OG_normalization/OG_glycoprotein_noTn_sl_tmm_HepG2.csv'
)

# HEK293T
OG_glycoprotein_noTn_sl_tmm_HEK293T <- read_csv(
  'data_source/OG_normalization/OG_glycoprotein_noTn_sl_tmm_HEK293T.csv'
)

# Jurkat
OG_glycoprotein_noTn_sl_tmm_Jurkat <- read_csv(
  'data_source/OG_normalization/OG_glycoprotein_noTn_sl_tmm_Jurkat.csv'
)

## WP normalization
# HepG2
WP_protein_sl_tmm_HepG2 <- read_csv(
  'data_source/WP_normalization/WP_protein_sl_tmm_HepG2.csv'
)

# HEK293T 
WP_protein_sl_tmm_HEK293T <- read_csv(
  'data_source/WP_normalization/WP_protein_sl_tmm_HEK293T.csv'
)

# Jurkat
WP_protein_sl_tmm_Jurkat <- read_csv(
  'data_source/WP_normalization/WP_protein_sl_tmm_Jurkat.csv'
)
