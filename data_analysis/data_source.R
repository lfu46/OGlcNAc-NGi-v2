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
  'data_source/raw_data/OG_psm_raw_Jurkat.csv'
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

# 