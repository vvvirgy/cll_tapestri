# stats

pt_data = lapply(pt_data, function(x) {
  x %>% 
    dplyr::mutate(diff = Bulk_CCF - Tapestri_CCF) %>% 
    dplyr::mutate(diff_abs = abs(diff)) %>% 
    dplyr::mutate(quality_ccf = 
                    case_when(diff_abs <= 0.1 ~ "good", 
                              diff_abs > 0.1 ~ "bad", 
                              is.na(diff_abs) ~ "unknown")) %>% 
    as_tibble() %>% 
    dplyr::mutate(CHROM = as.character(CHROM))
})

bind_rows(pt_data) %>%
  ggplot(aes(x = viber_cluster, fill = quality_ccf)) +
  geom_histogram(stat = "count", position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("bad" = "#C62300",
                               "good" = "#001A6E",
                               "unknown" = "gainsboro"),
                    labels = c("bad" = "more than 10% difference",
                               "good" = "less than 10% difference",
                               "unknown" = "mutation with bad quality or not present in tapestri")) +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(nrow=3)) +
  facet_wrap(~Case, scales = "free_x") 

ggsave("data/stats/ccfs_all_patients.png", width = 10, height = 10)

bind_rows(pt_data) %>%
  ggplot(aes(x = viber_cluster, fill = quality_ccf)) +
  geom_histogram(stat = "count", position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("bad" = "#C62300",
                               "good" = "#001A6E",
                               "unknown" = "gainsboro"),
                    labels = c("bad" = "more than 10% difference",
                               "good" = "less than 10% difference",
                               "unknown" = "mutation with bad quality or not present in tapestri")) +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(nrow=3)) +
  facet_wrap(~sample, scales = "free_x") 

ggsave("data/stats/ccfs_all_samples.png", width = 10, height = 10)


tot_by_cluster = bind_rows(pt_data) %>% 
  mutate(viber_cluster = case_when(viber_cluster == "S2_T2-" ~ "S2", 
                                   viber_cluster == "S3_T2+" ~ "S3", 
                                   viber_cluster == "S2_T2+" ~ "S2", 
                                   viber_cluster == "S1_T2-" ~ "S1",
                                   viber_cluster == "S3_T2-" ~ "S3", 
                                   .default = viber_cluster)) %>% 
  group_by(viber_cluster) %>% 
  count() %>% 
  rename(tot = n)
tot_by_cluster_quality = bind_rows(pt_data) %>% 
  mutate(viber_cluster = case_when(viber_cluster == "S2_T2-" ~ "S2", 
                                   viber_cluster == "S3_T2+" ~ "S3", 
                                   viber_cluster == "S2_T2+" ~ "S2", 
                                   viber_cluster == "S1_T2-" ~ "S1", 
                                   viber_cluster == "S3_T2-" ~ "S3", 
                                   .default = viber_cluster)) %>% 
  group_by(viber_cluster, quality_ccf, sample) %>% 
  count() 

full_join(tot_by_cluster, tot_by_cluster_quality, by = "viber_cluster") %>% 
  mutate(pt = (n/tot)*100) %>% 
  ggplot(aes(x = viber_cluster, y = pt, fill = quality_ccf)) + 
  geom_boxplot(outliers = F) +
  theme_bw() +
  # geom_jitter() + 
  scale_fill_brewer(palette = "Set1",
                    labels = c("bad" = "more than 10% difference",
                               "good" = "less than 10% difference",
                               "unknown" = "mutation with bad quality or not present in tapestri")) +
  theme(legend.position = "bottom") + 
  guides(fill=guide_legend(nrow=3)) + 
  ylab("% mutations")




x %>% 
  pivot_longer(cols = c(mutant_alleles, wt_alleles), names_to = "type", values_to = "number_of_alleles") %>%
  ggplot() +
  geom_bar(aes(y = driver_label, x = number_of_alleles, fill = type), stat = "identity", position = "stack") +
  facet_grid(gene~sample, scales = "free") +
  theme_light() +
  scale_fill_manual(values = type_colors) +
  geom_point(data = x %>% 
               filter(VAF < 0.05), 
             aes(x = 0, y = driver_label, color = VAF_filter), 
             shape = 8, size = 3) +
  scale_color_manual(values = vaf_col) +
  theme(axis.text.y = element_text(size = 10, vjust = 0.5), 
        legend.position = "bottom", legend)



lapply(pt_data %>% names, function(x) {
  tot = nrow(pt_data[[x]])
  pt_data[[x]] %>% 
    group_by(viber_cluster, quality_ccf) %>% 
    count()
})

