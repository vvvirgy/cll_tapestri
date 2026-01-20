library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(ggplotify)
library(rcartocolor)
# library(gt)
# library(dendextend)

setwd('/orfeo/scratch/cdslab/vgazziero/tapestri_cll')
source("cll_tapestri/R/utils.R")
source("cll_tapestri/R/clusters_colors.R")

# load data and prepare everything for the clustering
tapestri_data_path = "cll_tapestri/tapestri_results" # in which there are the folder by sample with all the extracted infos
bulk_data_path = "/orfeo/LTS/CDSLab/LT_storage/CLL/scDNA/h5_file/Tapestry CLL Targeted Panel.xlsx"
cz40_bulk_path= 'data/CZ40_split_cluster.rds'
source('/orfeo/cephfs/scratch/cdslab/vgazziero/tapestri_cll/cll_tapestri/R/data_preprocessing.R')


# now you can plot the heatmaps

heatmaps = lapply(names(pt_data), function(x) {
  
  generate_heatmap(patient_name = x, 
                   patient_data = pt_data[[x]], 
                   patient_genotypes = pt_genotypes[[x]], 
                   clusters_legend = colors_clusters_by_patient[[x]], 
                   ht_cols = cols, 
                   samples_cols = samples_cols, 
                   # path_ht = "genotyping_results/reports_heatmaps")
                   path_ht = NULL)

})

names(heatmaps) = names(pt_data)

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