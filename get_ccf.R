library(tidyverse)
library(readxl)

setwd("/orfeo/LTS/CDSLab/LT_storage/CLL/scDNA/h5_file/")
tapestri = readxl::read_excel("Tapestry CLL Targeted Panel.xlsx", sheet = 10)

metadata_tapestri = tapestri %>% 
  dplyr::mutate(ids = paste0(Chr, ":", Start, ":", Ref, "/", Var)) %>% 
  dplyr::select(Case, ids,viber_cluster, dplyr::starts_with("CCF"), dplyr::starts_with("VAF"))

x = metadata_tapestri %>% 
  dplyr::mutate(CCF_T3_neg = ifelse(CCF_T3_neg == "nd", NA, as.numeric(CCF_T3_neg))) %>% 
  dplyr::mutate(CCF_T3_pos = ifelse(CCF_T3_pos == "nd", NA, as.numeric(CCF_T3_pos))) %>% 
  dplyr::mutate(VAF_T3_neg = ifelse(VAF_T3_neg == "nd", NA, as.numeric(VAF_T3_neg))) %>% 
  dplyr::mutate(VAF_T3_pos = ifelse(VAF_T3_pos == "nd", NA, as.numeric(VAF_T3_pos))) %>% 
  # tidyr::pivot_longer(cols = c(starts_with("CCF"), starts_with("VAF"), starts_with("Filter")), 
  #                     names_to = c("sample_ccf", "sample_vcf", "sample_filter"), values_to = c("CCF", "VAF", "FILTER")) 
  tidyr::pivot_longer(cols = starts_with("CCF"), names_to = "sample_ccf", values_to = "CCF") %>% 
  dplyr::rename(sample = sample_ccf) %>% 
  dplyr::mutate(sample = paste0(Case, gsub("CCF", "", sample))) %>% 
  dplyr::select(sample, Case, ids, viber_cluster, CCF) %>% 
  dplyr::distinct()

  # tidyr::pivot_longer(cols = starts_with("VAF"), names_to = "sample_vaf", values_to = "VAF") %>% 
  # tidyr::pivot_longer(cols = starts_with("FILTER"), names_to = "sample_filter", values_to = "FILTER")

y = metadata_tapestri %>% 
  dplyr::mutate(CCF_T3_neg = ifelse(CCF_T3_neg == "nd", NA, as.numeric(CCF_T3_neg))) %>% 
  dplyr::mutate(CCF_T3_pos = ifelse(CCF_T3_pos == "nd", NA, as.numeric(CCF_T3_pos))) %>% 
  dplyr::mutate(VAF_T3_neg = ifelse(VAF_T3_neg == "nd", NA, as.numeric(VAF_T3_neg))) %>% 
  dplyr::mutate(VAF_T3_pos = ifelse(VAF_T3_pos == "nd", NA, as.numeric(VAF_T3_pos))) %>% 
  # tidyr::pivot_longer(cols = c(starts_with("CCF"), starts_with("VAF"), starts_with("Filter")), 
  #                     names_to = c("sample_ccf", "sample_vcf", "sample_filter"), values_to = c("CCF", "VAF", "FILTER")) 
  tidyr::pivot_longer(cols = starts_with("VAF"), names_to = "sample_vaf", values_to = "VAF") %>% 
  dplyr::rename(sample = sample_vaf) %>% 
  dplyr::mutate(sample = paste0(Case, gsub("VAF", "", sample))) %>% 
  dplyr::select(Case, ids, viber_cluster, sample, VAF) %>% 
  dplyr::distinct()
  
new_metadata = full_join(x,y, by = join_by("Case" == "Case", "ids" == "ids", "viber_cluster" == "viber_cluster", "sample" == "sample"))

write.table(new_metadata, file = "/orfeo/scratch/cdslab/vgazziero/tapestri_cll/data/mutations_ccfs.csv", sep = ",", quote = F, row.names = F, col.names = T)
