library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(ggplotify)
library(rcartocolor)
library(gt)

source("cll_tapestri/00.auxiliary_functions.R")

samples = list.files("cll_tapestri/tapestri_results")

tapestri_results = lapply(samples, function(x) {
  list.files("cll_tapestri/tapestri_results", full.names = T, pattern = x, recursive = T) 
})
names(tapestri_results) = samples

tapestri_results = lapply(tapestri_results, function(x) {
  lapply(x, function(file) {
    read.table(file, header = T, sep = ",")
  })
})

tapestri_results = lapply(tapestri_results, function(x) {
  names(x) = c("cells_info", "CCF_comparison", "AF", "DP", "NGT", "variants_info")
  return(x)
}) 

# filter the results to get only good mutations and cells, then recalculate the CCFs

tapestri_results_filtered_NGT = lapply(tapestri_results, function(x) {
  good_barcodes = x$cells_info %>% filter(filtered == 0) %>% pull(barcode)
  good_variants = x$variants_info %>% filter(filtered == 0) %>% pull(id)
  
  colnames(x$NGT) = replace_dots(colnames(x$NGT))
  x$NGT = x$NGT %>%
    column_to_rownames(var = "X")
  x$NGT = x$NGT[good_barcodes, good_variants]
  
  colnames(x$AF) = replace_dots(colnames(x$AF))
  x$AF = x$AF %>%
    column_to_rownames(var = "X")
  x$AF = x$AF[good_barcodes, good_variants]
  
  colnames(x$DP) = replace_dots(colnames(x$DP))
  x$DP = x$DP %>%
    column_to_rownames(var = "X")
  x$DP = x$DP[good_barcodes, good_variants]
  
  # calculate Tapestri CCF
  filtered_ccfs = apply(x$NGT, 2, function(x) {
    all_locus = length(x)
    
    tot_muts = as.data.frame(table(x)) %>% 
      filter(x %in% c(1,2)) %>% 
      pull(Freq) %>% 
      sum
    tot_muts/all_locus
  }) %>% 
    as.data.frame() %>% 
    rownames_to_column("mut_id") %>% 
    rename(tapestri_ccf = ".")
  
  x$CCF_comparison = x$CCF_comparison %>% 
    dplyr::mutate(X = NULL, Tapestri_id = NULL) %>% 
    dplyr::full_join(., filtered_ccfs, by = join_by("mutation" == "mut_id")) %>% 
    dplyr::select(mutation, starts_with("CCF"), tapestri_ccf) %>% 
    dplyr::rename("Bulk_CCF" = starts_with("CCF")) %>% 
    dplyr::rename("Tapestri_CCF" = tapestri_ccf) %>% 
    dplyr::distinct()
  
  return(x)
})

ccf_comparison = lapply(tapestri_results_filtered_NGT, function(x) {
  
  tt = x$CCF_comparison
  
  infos = x$variants_info %>% 
    # dplyr::filter(filtered == 0) %>% 
    dplyr::mutate(X = NULL)
  
  dplyr::full_join(tt, infos, by = join_by("mutation" == "id")) %>% 
    dplyr::mutate(diff_ccf = Bulk_CCF-Tapestri_CCF)
})

ccf_comparison$CT525_T2_neg = ccf_comparison$CT525_T2_neg %>% 
  mutate(Bulk_AF = NULL)

# is ccf similar between tapestri data and bulk ones?

# plots = lapply(ccf_comparison %>% names, function(x) {
#   [[x]] %>% 
#     mutate(quality = ifelse(filtered == 1, "low", "high")) %>% 
#     mutate(quality = factor(quality, levels = c("low", "high"))) %>% 
#     ggplot(aes(Tapestri_CCF,Bulk_CCF, colour = quality)) + 
#     geom_point() + 
#     theme_bw() + 
#     scale_color_brewer(palette = "Set1") + 
#     ggtitle("Bulk CCF vs scDNA CCF")
# })

# names(plots) = ccf_comparison %>% names

# wrap_plots(plots,guides = "collect") & theme(legend.position = "bottom")

# 
# ccf_comparison$CT339_T2_neg %>%  
#   ggplot(aes(y = reorder(mutation, +CCF), x = CCF)) + 
#   # geom_segment() +
#   geom_line(aes(group = mutation, color = filtered)) + 
#   geom_point(aes(color = type)) + 
#   theme_bw() 

dp_stats = lapply(tapestri_results_filtered_NGT, function(x) {
  dp_mean = apply(x$DP[,-1], 2, mean)
  dp_sd = apply(x$DP[,-1], 2, sd)
  
  tibble(id = names(dp_mean), 
         DP_mean = dp_mean, 
         DP_sd = dp_sd)
})

dp_stats_plot = lapply(dp_stats %>% names, function(x) {
  dp_stats[[x]] %>% 
    ggplot(aes(DP_mean)) +
    geom_histogram(binwidth = 3) + 
    theme_bw() + 
    ggtitle("mean DP")
}) 
names(dp_stats_plot) = names(dp_stats)

full_table = lapply(names(ccf_comparison), function(x) {
  ccf_comparison[[x]] %>% 
    # mutate(tmp_id = gsub("\\:", ".", mutation)) %>% 
    # mutate(tmp_id = gsub("/", ".", tmp_id)) %>% 
    full_join(., dp_stats[[x]], by = join_by("mutation" == "id")) %>% 
    mutate(quality = ifelse(filtered == 1, "low", "high")) %>% 
    mutate(quality = factor(quality, levels = c("low", "high"))) %>% 
    mutate(tmp_id = NULL)
})
names(full_table) = names(ccf_comparison)

# stats
# quality_vs_dp = lapply(full_table, function(x) {
#   x %>% 
#     ggplot(aes(quality, DP_mean, fill = quality)) +
#     geom_violin() +
#     theme_bw() +
#     scale_fill_brewer(palette = "Set1")
# })

excel_file = readxl::read_excel("/orfeo/LTS/CDSLab/LT_storage/CLL/scDNA/h5_file/Tapestry CLL Targeted Panel.xlsx", sheet = 10)
# excel_file = readxl::read_excel("~/Desktop/LTS_cdslab/Tapestry CLL Targeted Panel.xlsx", sheet = 10)
# excel_file = readxl::read_excel("~/Desktop/cdslab_lts/vgazziero/Tapestry CLL Targeted Panel.xlsx", sheet = 10)

excel_file = excel_file %>% 
  dplyr::mutate(nd...1 = NULL, nd...41 = NULL) %>% 
  dplyr::mutate(mut_id = paste0(Chr, ":", Start, ":", Ref, "/", Var))

ccf_comparison_v2 = lapply(full_table %>% names, function(x) {
  # colnames(ccf_comparison[[x]]) = c("mutation", "Tapestri_CCF", "Bulk_CCF")
  tt = full_table[[x]]
  
  ex = excel_file %>% 
    filter(Case == unlist(strsplit(x, "_")) %>% first()) %>% 
    dplyr::distinct(mut_id, .keep_all = TRUE)
  
  left_join(tt, ex, by = join_by("mutation" == "mut_id"))
})

names(ccf_comparison_v2) = names(full_table)

ccf_by_cluster = lapply(ccf_comparison_v2, function(x) {
  x %>% 
    ggplot(aes(x = viber_cluster, y = diff_ccf)) + 
    geom_boxplot(outliers = FALSE) + 
    geom_jitter() +
    scale_color_brewer(palette = "Set1") +
    theme_bw() + 
    ggtitle("diff CCF (Bulk - Tapestri)")
})

  # pivot_longer(cols = c(Tapestri_CCF, Bulk_CCF), names_to = "experiment", values_to = "CCF") %>% 
  
  # summarise(dist_mean = mean(diff_ccf), dist_sd = sd(diff_ccf)) 

# genotype heatmap per sample 

patients = tapestri_results %>% names %>% 
  gsub("_neg|_pos", "", .) %>% 
  # gsub("_T[12]$", "", .) %>% 
  unique

viber_cluster_colors = lapply(ccf_comparison_v2, function(x) {x$viber_cluster %>% unique}) %>% unlist %>% unique
# viber_cluster_colors = excel_file %>% filter(Case %in% patients) %>% pull(viber_cluster) %>% unique
viber_cluster_colors = setNames(carto_pal(n = length(viber_cluster_colors), name = "Prism"), viber_cluster_colors)

x = patients %>% first

gnt_ht = lapply(patients, function(x) {
  
  print(x)

  # patient_muts = excel_file %>% 
  #   filter(Case == x)
  
  samples_to_check = names(tapestri_results_filtered_NGT) %>% grep(x, . , value = T)
  mutations = ccf_comparison_v2[samples_to_check] %>% 
    bind_rows() %>% 
    filter(filtered == 0) %>% 
    select(viber_cluster, mutation) %>% 
    distinct()
  
  mut_order = setNames(nm = mutations %>% arrange(desc(viber_cluster)) %>% pull(mutation), 
                        object = mutations %>% arrange(desc(viber_cluster)) %>% pull(viber_cluster))
  
  data = lapply(samples_to_check, function(s) {
    df = tapestri_results_filtered_NGT[[s]]$NGT %>% 
      mutate(across(everything(), ~ case_when(
        . == 0 ~ "",
        . == 1 ~ "het-mut",
        . == 2 ~ "hom-mut",
        . == 3 ~ "unknown"
      ))) 
    rownames(df) = paste(s, seq(1:nrow(df)), sep = "_")
     
    df %>% t %>% as.data.frame() %>% 
      rownames_to_column(var = "mut")
  }) %>% 
    Reduce(full_join, .)

  data = column_to_rownames(data, var = "mut")
  data = as.matrix(data)
  
  # cols = setNames(object = c("gainsboro", "#C62300", "grey95"), c(0,1,3))
  cols = setNames(object = c("#C62300", "#0A3981", "#F5F5F5"),c("het-mut", "hom-mut", "unknown"))
  # ann_cols = setNames(object = c('#48A6A7', '#3E7B27', '#FF9D23'), unname(mut_order) %>% unique)
  # ann_cols = setNames(object = awtools::mpalette[1:length(mut_order %>% unique)], unname(mut_order) %>% unique)
  ann_cols = viber_cluster_colors[unname(mut_order) %>% unique]
  
  # ccf_values_by_cluster = setNames(nm = patient_muts %>% arrange(desc(viber_cluster)) %>% pull(Bulk_CCF) %>% unique, 
  #                                  object = patient_muts %>% arrange(desc(viber_cluster)) %>% pull(viber_cluster) %>% unique)
  
  ha = rowAnnotation(viber_cluster = unname(mut_order),
                         col = list(viber_cluster = ann_cols)) #, 
                         # annotation_legend_param = list(at = unique(unname(mut_order)), 
                         #                                labels = paste0(unique(unname(mut_order)), " (CCF = ", 
                         #                                               round(ccf_values_by_cluster %>% names %>% as.numeric, 2), ")"))) 
  
  
  all_samples = gsub("_[0-9]{1,5}$", "", colnames(data))
  sample_cols = setNames(c("#7EB5A6", "#2B580C"), unique(all_samples))
  
  sample_ann = HeatmapAnnotation(sample = all_samples, col = list(sample = sample_cols))
  
  data = data[names(mut_order), ]
  pdf("test.pdf")
  ComplexHeatmap::Heatmap(data[1:10,1:100], 
                          # alter_fun = gt_fun, 
                          col = cols, 
                          cluster_columns = T,
                          show_column_dend = T,
                          cluster_rows = T, show_row_dend = T, 
                          # left_annotation = ha, 
                          show_row_names = T, 
                          show_column_names = F)
                          # row_split = factor(unname(mut_order), levels =unname(mut_order) %>% unique), 
                          # column_split = factor(all_samples, levels = all_samples %>% unique), 
                          # top_annotation = sample_ann)
  # dev.off()
  graphics.off()
  
  ComplexHeatmap::Heatmap(data, 
                            # alter_fun = gt_fun, 
                            col = cols, 
                            left_annotation = ha, show_row_names = T, 
                            # row_split = factor(unname(mut_order), levels =unname(mut_order) %>% unique), 
                            # column_split = factor(all_samples, levels = all_samples %>% unique), 
                            top_annotation = sample_ann, cluster_columns = T, cluster_rows = F, 
                            # show_pct = F, 
                            name = "genotypes", 
                            row_title = "mutations",
                            column_title = "cells", 
                            row_names_gp = gpar(fontsize = 9)) 
  
  
  gt_fun = lapply(cols, function(x) { 
    alter_graphic("rect", width = 0.9, height = 0.9, fill = unname(x))
  })

  gt_fun = c("background" = list(alter_graphic("rect", width = 0.9, height = 0.9, fill = "#CCCCCC")), gt_fun)

  # pdf("test.pdf")
  # ComplexHeatmap::Heatmap(data, 
  #                         col = cols, 
  #                         name = "genotypes", 
  #                         # top_annotation = ha, 
  #                         cluster_columns = F, 
  #                         cluster_rows = F, 
  #                         column_title = "mutations", 
  #                         row_title = "cells", 
  #                         column_names_gp = gpar(fontsize = 9))
  # graphics.off()
  
  gt = ComplexHeatmap::oncoPrint(data, 
                                 alter_fun = gt_fun, 
                                 col = cols, 
                                 left_annotation = ha, 
                                 row_split = factor(unname(mut_order), levels =unname(mut_order) %>% unique), 
                                 column_split = factor(all_samples, levels = all_samples %>% unique), 
                                 top_annotation = sample_ann, 
                                 show_pct = F, 
                                 name = "genotypes", 
                                 row_title = "mutations",
                                 column_title = "cells", 
                                 row_names_gp = gpar(fontsize = 9)) 
  
  return(gt)
})

names(gnt_ht) = patients

dir.create("genotyping_results")

lapply(gnt_ht %>% names, function(x) {
  pdf(paste0("genotyping_results/clusters_genotyping_",x,".pdf"), width = 15, height = 12)
  draw(gnt_ht[[x]], heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  graphics.off()
})

gnt_ht_v2 = lapply(gnt_ht, function(x) {
  draw(x, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
})

# CCFS 

summary_CCFs = lapply(ccf_comparison_v2, function(x) {
  
  x %>% 
    dplyr::filter(!is.na(Tapestri_CCF)) %>% 
    dplyr::group_by(viber_cluster, Bulk_CCF) %>% 
    dplyr::summarise(mean_ccf = mean(Tapestri_CCF)) %>% 
    dplyr::ungroup()
})

summary_CCFs = lapply(patients, function(x) {
  
  samples_to_check = names(summary_CCFs) %>% grep(x, . , value = T)
  
  lapply(summary_CCFs[samples_to_check] %>% names, function(x) {
    summary_CCFs[[x]] %>% 
      dplyr::mutate(sample = x)
  }) %>% bind_rows() %>% 
    tidyr::pivot_wider(names_from = sample, values_from = c(Bulk_CCF, mean_ccf))
  
})

names(summary_CCFs) = patients

ccfs_tables = lapply(summary_CCFs %>% names, function(x) {
  summary_CCFs[[x]] %>% 
    gt(rowname_col = "viber_cluster") %>% 
    gt::fmt_number(decimals = 3) %>% 
    gt::tab_stubhead(label = "Cluster") %>% 
    gt::tab_spanner(label = "positive", 
                    columns = ends_with("pos")) %>% 
    gt::tab_spanner(label = "negative", 
                    columns = ends_with("neg")) %>% 
    gt::cols_label_with(starts_with("Bulk_CCF"), fn = function(x) {
      paste(str_extract(string = x, pattern = "Bulk"), "CCF")
    }) %>% 
    gt::cols_label_with(starts_with("mean_ccf"), fn = function(x) {
      paste0("Tapestri CCF (", str_extract(x, pattern = "mean"), ")")
    }) %>% 
    gt::tab_header(title = gsub("_", " ", x)) %>% 
    gt::as_gtable(. , text_grob = gridtext::richtext_grob, plot = T)
})
names(ccfs_tables) = names(summary_CCFs)

layout <- c(
  area(t = 1, l = 2, b = 2, r = 4), 
  area(t = 3, l = 1, b = 7, r = 5)
)

reports = lapply(names(ccfs_tables), function(x) {
  wrap_plots(list(ccfs_tables[[x]], ggplotify::as.ggplot(gnt_ht[[x]])), design = layout) + plot_annotation(title = x)
})

names(reports) = names(ccfs_tables)

lapply(names(reports), function(x) {
  pdf(paste0("genotyping_results/", x, "_report.pdf"), width = 15, height = 10)
  print(reports[[x]])
  dev.off()
})

# lapply(names(ccf_by_cluster), function(x) {
#   
#   # gg = draw(,
#   #           heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
#   
#   pp = list(ggplotify::as.ggplot(gnt_ht[[x]]), ccf_by_cluster[[x]])
#   
#   ll = "
#   AAAB
#   AAA#
#   "
#   
#   p = wrap_plots(pp, design = ll) +
#     plot_annotation(title = x) & theme(legend.position = "bottom")
#   
#   ggsave(paste0("data/reports/patient_", x, "_ccf_comparison.png"), p, width = 15, height = 10)
#   
# })

 
# ggplotify::as.ggplot(gt) + ccf_by_cluster[[1]] +
#   plot_annotation(title = "CT339_T2_neg")

# ggsave("data/test_patient_CT339_T2_neg.png", width = 20, height = 15)

reports = lapply(names(dp_stats_plot), function(x) {
  pp = list(dp_stats_plot[[x]]
            quality_vs_dp[[x]],
            plots[[x]], 
            ccf_by_cluster[[x]]
  )
  
  fp = wrap_plots(pp, guides = "collect") + 
    plot_annotation(title = x) & 
    theme(legend.position = "bottom")
  
  ggsave(filename = paste0("data/", x, "report.png"), plot = fp)
})



pp = lapply(ccf_comparison_v2 %>% names, function(x) {
  # colnames(ccf_comparison[[x]]) = c("mutation", "Tapestri_CCF", "Bulk_CCF")
  ccf_comparison_v2[[x]] %>% 
    ggplot(aes(x = Tapestri_CCF, y = Bulk_CCF, colour = viber_cluster)) +
    geom_point() +
    ggtitle(x) + 
    theme_bw() +
    theme(legend.position = "bottom")
})

names(pp) = names(ccf_comparison_v2)
wrap_plots(pp) & theme(legend.position = 'bottom')

wrap_plots(list(pp$CT339_T2_neg, pp$CT339_T2_pos), guides = "collect") & theme(legend.position = 'bottom')

ggsave("data/ct399_pos_neg_ccfs_tap_excel.png", width = 8, height = 5)

# ccf_comparison_v2$CT339_T2_neg %>% 
#   ggplot(aes(x = Tapestri_CCF_CT339_T2_neg, y = CCF_T2_neg.x, colour = viber_cluster)) +
#   # geom_point() +
#   ggtitle("CT339_T2_neg")
# 
# 
# ccf_comparison$CT339_T2_neg %>% 
#   ggplot(aes(Tapestri_CCF_CT339_T2_neg)) +
#   geom_histogram(binwidth = 0.01) +
#   xlim(c(-0.01,1.01))
