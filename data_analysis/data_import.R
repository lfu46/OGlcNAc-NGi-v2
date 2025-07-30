# import packages
library(tidyverse)

# import O-GlcNAcylation result
OG_psm_raw_HepG2 <- read_csv(
  '/Volumes/Expansion/Glycosylation_Crosstalk/Data/NOCT_HepG2_07042024/OG_sequest/ronghuwulab_1729609351.csv',
  col_names = TRUE,
  name_repair = "universal"
)

OG_psm_raw_HEK293T <- read_csv(
  '/Volumes/Expansion/Glycosylation_Crosstalk/Data/NOCT_HEK293T_08102024/OG_sequest/ronghuwulab_1729615074.csv',
  col_names = TRUE,
  name_repair = "universal"
)

OG_psm_raw_Jurkat <- read_csv(
  '/Volumes/Expansion/Glycosylation_Crosstalk/Data/NOCT_Jurkat_08102024/OG_sequest/ronghuwulab_1729615323.csv',
  col_names = TRUE,
  name_repair = "universal"
)

# import whole proteome result
WP_psm_raw_HepG2 <- read_csv(
  '',
  col_names = TRUE,
  name_repair = "universal"
)

WP_psm_raw_HEK293T <- read_csv(
  '/Volumes/Expansion/LocalDisk_E/LocalDisk_1/Glycosylation_Crosstalk/Data/NOCT_HEK293T_08102024/WP_sequest/ronghuwulab_1725719351.csv',
  col_names = TRUE,
  name_repair = "universal"
)

WP_psm_raw_Jurkat <- read_csv(
  '/Volumes/Expansion/LocalDisk_E/LocalDisk_1/Glycosylation_Crosstalk/Data/NOCT_Jurkat_08102024/WP_sequest/ronghuwulab_1725717952.csv',
  col_names = TRUE,
  name_repair = "universal"
)
