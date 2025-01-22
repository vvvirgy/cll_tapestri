library(tidyverse)
library(patchwork)

files = list.files("cll_tapestri/excel_positions", recursive = T, full.names = T)

ccf_comparison = lapply(files, function(x) {
  read.table(x, header = T, sep = ",")[,-1] %>% 
    dplyr::select(everything(), -Tapestri_id) %>% 
    dplyr::relocate(mutation, .before = starts_with("Tapestri_CCF*")) 
})

samples = lapply(strsplit(files, split = "/"), function(x) {
  x %>% 
    last() %>%  
    gsub("_comparison_tapestri_vs_bulk_excel_positions.csv", "", .)
}) %>% unlist()

names(ccf_comparison) = samples



pp = lapply(ccf_comparison %>% names, function(x) {
  colnames(ccf_comparison[[x]]) = c("mutation", "Tapestri_CCF", "Bulk_CCF")
  ccf_comparison[[x]] %>% 
    ggplot(aes(x = Tapestri_CCF, y = Bulk_CCF)) +
    geom_point() +
    ggtitle(x)
})

names(pp) = names(ccf_comparison)
wrap_plots(pp)

excel_file = readxl::read_excel("~/Desktop/cdslab_lts/CLL/scDNA/h5_file/Tapestry CLL Targeted Panel.xlsx", sheet = 10)

excel_file = excel_file %>% 
  dplyr::mutate(mut_id = paste0(Chr, ":", Start, ":", Ref, "/", Var))

ccf_comparison_v2 = lapply(ccf_comparison %>% names, function(x) {
  tt = ccf_comparison[[x]]
  
  ex = excel_file %>% 
    filter(Case == unlist(strsplit(x, "_")) %>% first())
  
  left_join(tt, ex, by = join_by("mutation" == "mut_id"))
})

names(ccf_comparison_v2) = names(ccf_comparison)

ccf_comparison_v2$CT339_T2_neg %>% 
  ggplot(aes(x = Tapestri_CCF_CT339_T2_neg, y = CCF_T2_neg.x, colour = viber_cluster)) +
  # geom_point() +
  ggtitle("CT339_T2_neg")


ccf_comparison$CT339_T2_neg %>% 
  ggplot(aes(Tapestri_CCF_CT339_T2_neg)) +
  geom_histogram(binwidth = 0.01) +
  xlim(c(-0.01,1.01))
