rm(list=ls())

library(tidyverse)
library(readxl)

# read excel file 
excel_file =  read_excel("~/Desktop/LTS_cdslab/Tapestry CLL Targeted Panel.xlsx", sheet = 10)

# read the rds file
rds_mutations = readRDS("~/Desktop/scratch_cdslab/tapestri_cll/data/mutations_clusters_table.rds")


# format them similarly 

excel_file = excel_file %>% 
  mutate(nd...1 = NULL) %>% 
  distinct()

rds_mutations_v2 = rds_mutations %>% 
  mutate(sample = case_when(sample == "1p" ~ "T1_pos", 
                            sample == "1n" ~ "T1_neg", 
                            sample == "2p" ~ "T2_pos", 
                            sample == "2n" ~ "T2_neg", 
                            sample == "3p" ~ "T3_pos", 
                            sample == "3n" ~ "T3_neg", .default = sample)) %>% 
  group_by(patient) %>% 
  group_split()

names(rds_mutations_v2) = lapply(rds_mutations_v2, function(x) { x$patient %>% unique }) %>% unlist

excel_file = excel_file %>% 
  group_by(Case) %>% 
  group_split()

names(excel_file) = lapply(excel_file, function(x) {x$Case %>% unique}) %>% unlist

rds_mutations_v2_noNA = lapply(rds_mutations_v2, function(x) {
  x %>% 
    distinct() %>% 
    # filter(!is.na(CCF)) %>% 
    pivot_wider(names_from = sample, values_from = CCF, names_prefix = "CCF_", values_fill = 0)
}) 

rds_mutations_v2_noNA = rds_mutations_v2_noNA %>% 
  purrr::list_modify("CT287" = NULL)

results = lapply(names(excel_file), function(x) {
  muts_excel = excel_file[[x]]$id
  muts_rds = rds_mutations_v2_noNA[[x]]$id
  
  # print(x)
  
  tibble(sample = x, 
         intersection_excel_rds = intersect(unique(muts_excel), unique(muts_rds)) %>% length(), 
         excel = length(unique(muts_excel)))
}) %>% bind_rows()
results

