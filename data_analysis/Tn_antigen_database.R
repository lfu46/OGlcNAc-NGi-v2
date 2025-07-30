# import packages
library(tidyverse)
library(readxl)

### generate Tn antigen database
## import Tn data 
# Anal. Chem. 2022, 94, 7, 3343–3351.: https://pubs.acs.org/doi/full/10.1021/acs.analchem.1c05438
# supporting information: https://pubs.acs.org/doi/suppl/10.1021/acs.analchem.1c05438/suppl_file/ac1c05438_si_003.xlsx
Tn_AnalChem2022_Exp1 <- read_xlsx(
  'data_source/Tn_antigen_database/ac1c05438_si_004.xlsx',
  skip = 1,
  col_names = TRUE,
  .name_repair = 'universal',
  sheet = 'Jurkat Exp#1'
) |> 
  separate(Protein...Name, into = c('sp', 'Uniprot.Entry', 'HUMAN'), sep = '\\|')

Tn_AnalChem2022_Exp2 <- read_xlsx(
  'data_source/Tn_antigen_database/ac1c05438_si_004.xlsx',
  skip = 1,
  col_names = TRUE,
  .name_repair = 'universal',
  sheet = 'Jurkat Exp#2'
) |> 
  separate(Protein...Name, into = c('sp', 'Uniprot.Entry', 'HUMAN'), sep = '\\|')

# Angew. Chem. Int. Ed. Engl. 2017;56(25):7107-7111.: https://onlinelibrary.wiley.com/doi/full/10.1002/anie.201702191
# supporting information: https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fanie.201702191&file=anie201702191-sup-0001-misc_information.pdf
Tn_Angew2017_Exp1 <- read_xlsx(
  'data_source/Tn_antigen_database/Table_S1.xlsx',
  skip = 1,
  col_names = TRUE,
  .name_repair = 'universal'
)

Tn_Angew2017_Exp2 <- read_xlsx(
  'data_source/Tn_antigen_database/Table_S2.xlsx',
  skip = 1,
  col_names = TRUE,
  .name_repair = 'universal'
)

Tn_Angew2017_Exp3 <- read_xlsx(
  'data_source/Tn_antigen_database/Table_S3.xlsx',
  skip = 1,
  col_names = TRUE,
  .name_repair = 'universal'
)

# combine all results for protein with Tn antigen
Tn_antigen_database <- bind_rows(
  Tn_AnalChem2022_Exp1 |> select(Uniprot.Entry),
  Tn_AnalChem2022_Exp2 |> select(Uniprot.Entry),
  Tn_Angew2017_Exp1 |> select(Uniprot.Entry),
  Tn_Angew2017_Exp2 |> select(Uniprot.Entry),
  Tn_Angew2017_Exp3 |> select(Uniprot.Entry)
) |> 
  distinct()

write_csv(
  Tn_antigen_database,
  file = 'data_source/Tn_antigen_database/Tn_antigen_database.csv'
)

