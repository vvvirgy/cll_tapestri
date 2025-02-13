library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(ggplotify)
library(rcartocolor)
# library(gt)
# library(dendextend)

source("cll_tapestri/00.auxiliary_functions.R")
source("cll_tapestri/clusters_colors.R")

tapestri_data = load_tapestri_data(path = "cll_tapestri/tapestri_results")

# filter the results to get only good mutations and cells, then recalculate the CCFs
tapestri_filtered = lapply(tapestri_data, function(x) {filter_tapestri(x, filter_cells = T, filter_ids = T)})

# compute CCFs after the filtering
tapestri_filtered = lapply(tapestri_filtered, compute_tapestri_CCF)

# load the bulk information
bulk_data = readxl::read_excel("/orfeo/LTS/CDSLab/LT_storage/CLL/scDNA/h5_file/Tapestry CLL Targeted Panel.xlsx", sheet = 10)
# excel_file = readxl::read_excel("~/Desktop/LTS_cdslab/Tapestry CLL Targeted Panel.xlsx", sheet = 10)
# excel_file = readxl::read_excel("~/Desktop/cdslab_lts/vgazziero/Tapestry CLL Targeted Panel.xlsx", sheet = 10)

bulk_data = bulk_data %>% 
  dplyr::mutate(nd...1 = NULL, nd...41 = NULL) %>% 
  dplyr::mutate(mut_id = paste0(Chr, ":", Start, ":", Ref, "/", Var))

cz40_bulk = readRDS("data/CZ40_split_cluster.rds")
cz40_bulk = cz40_bulk %>% 
  tidyr::separate(id, into = c('chr', 'from', 'to', 'ref', 'alt'), sep = ":") %>% 
  dplyr::mutate(id = paste0(chr, ":", from, ":", ref, "/", alt)) %>% 
  dplyr::select(id, viber_cluster, phylo_color) %>% 
  dplyr::mutate(Case = "CZ40") %>% 
  dplyr::rename(new_viber_cluster = viber_cluster)

bulk_data = bulk_data %>% 
  dplyr::group_by(Case) %>% 
  dplyr::group_split()
names(bulk_data) = lapply(bulk_data, function(x) {x$Case %>% unique}) %>% unlist

tt = bulk_data$CZ40 %>% 
  dplyr::inner_join(., cz40_bulk, by = join_by("mut_id" == "id", "Case" == "Case")) %>% 
  dplyr::mutate(viber_cluster = new_viber_cluster) %>% 
  dplyr::mutate(new_viber_cluster = NULL, phylo_color = NULL)

bulk_data$CZ40 = bulk_data$CZ40 %>% 
  filter(!mut_id %in% tt$mut_id) %>% 
  bind_rows(., tt)

bulk_data$RM238 = bulk_data$RM238 %>% 
  dplyr::filter(viber_cluster != "S2")

bulk_data = bulk_data %>% 
  do.call(rbind, .)

# save "S3_T2+"

bulk_data = bulk_data %>% 
  dplyr::mutate(viber_cluster = case_when(viber_cluster == "S2_T2-" ~ "S2", 
                                          viber_cluster == "S1_T2-" ~ "S1", 
                                          viber_cluster == "S2_T2+" ~ "S2", 
                                          viber_cluster == "S2_T1+" ~ "S2", 
                                          viber_cluster == "S1_T1-" ~ "S1", 
                                          viber_cluster == "S1_T1+" ~ "S2", 
                                          viber_cluster == "S3_T2-" ~ "S3", 
                                          viber_cluster == "S3_T1+" ~ "S3", 
                                          viber_cluster == "S5_T3+" ~ "S5", 
                                          .default = viber_cluster))

all_patients_time_points = all_patients = tapestri_filtered %>% 
  names %>% 
  gsub("_neg|_pos", "", .) %>%   
  unique

pt_sampl = dplyr::tibble(
  patient = all_patients_time_points %>% gsub("_T[2|3]$", "", .), 
  samples_t = all_patients_time_points
)

pt_data = lapply(pt_sampl$samples_t, function(x) {
  p = pt_sampl %>% 
    filter(samples_t == x) %>% 
    pull(patient)
  print(p)
  get_patient_data(tapestri_results = tapestri_filtered, patient = p, sample_t = x, bulk_data = bulk_data)
})

names(pt_data) = pt_sampl$samples_t

pt_genotypes = lapply(pt_sampl$samples_t, function(x) {
  get_patient_genotypes(tapestri_results = tapestri_filtered, patient = x)
})
names(pt_genotypes) = pt_sampl$samples_t


lapply(names(pt_data), function(x) {
  
  generate_heatmap(patient_name = x, 
                   patient_data = pt_data[[x]], 
                   patient_genotypes = pt_genotypes[[x]], 
                   clusters_legend = colors_clusters_by_patient[[x]], 
                   ht_cols = cols, 
                   samples_cols = samples_cols, 
                   path_ht = "genotyping_results/reports_heatmaps")

})

# stats 

# dir.create("data/stats")
# ccf_stats = lapply(tapestri_results_filtered_NGT %>% names, function(x) {
#   tt = tapestri_results_filtered_NGT[[x]]$CCF_comparison %>% 
#     dplyr::mutate(diff = Bulk_CCF - Tapestri_CCF) %>% 
#     dplyr::mutate(diff_abs = abs(diff)) %>% 
#     dplyr::mutate(quality_ccf = 
#                     case_when(diff_abs <= 0.1 ~ "less than 0.1 difference", 
#                               diff_abs > 0.1 ~ "more than 0.1 difference", 
#                               is.na(diff_abs) ~ "mutation not found or with bad quality in tapestri data"))
#     
#   
#   pt = tt %>%
#     pivot_longer(cols = c("Tapestri_CCF", "Bulk_CCF"), names_to = "source", values_to = "CCF") %>%
#     ggplot(aes(y = reorder(mutation, +diff), x = CCF)) +
#     # geom_segment() +
#     geom_line(aes(group = mutation, color = quality_ccf), alpha = 0.7, size = 1.5) +
#     geom_point(aes(color = source), size = 2) +
#     theme_bw() + 
#     scale_color_manual(values = c("more than 0.1 difference" = "#DF2E38", 
#                                   "less than 0.1 difference" = "#399918", 
#                                   "mutation not found or with bad quality in tapestri data" = "gainsboro",
#                                   "Bulk_CCF" = "#7E5CAD", 
#                                   "Tapestri_CCF" = "#578FCA")) + 
#     facet_wrap(~quality_ccf, scales = "free_y") + 
#     theme(legend.position = "bottom") +
#     ylab("Delta CCF (Bulk - Tapestri)") + 
#     ggtitle(x)
# 
#   ggsave(filename = paste0("data/stats/", x, "_diff_ccf.pdf"), width = 12, height = 6)  
# })



# get_patient_mutations(x, patient, ccf_table){
#   
#   samples = names(x) %>% grep(x, . , value = T)
#   
#   mutations = ccf_comparison_v2[samples] %>% 
#     bind_rows() %>% 
#     filter(filtered == 0) %>% 
#     select(viber_cluster, mutation) %>% 
#     distinct()
# }