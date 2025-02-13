# to uniform all the strings in the mutation ids

replace_dots <- function(input_string) {
  # Replace all dots with colons
  modified_string <- str_replace_all(input_string, "\\.", ":")
  
  # Replace the last colon with a slash
  modified_string <- str_replace(modified_string, ":(?=[^:]*$)", "/")
  
  return(modified_string)
}

# load the results 
load_tapestri_data = function(path, colnames = c("cells_info", "CCF_comparison", "AF", "DP", "NGT", "variants_info")) {
  
  # get all the samples
  samples = list.files(path)
  
  # get the files
  tapestri_results = lapply(samples, function(x) {
    list.files(path, full.names = T, pattern = x, recursive = T) 
  })
  names(tapestri_results) = samples
  
  # read the data and name the elements of the list
  tapestri_results = lapply(tapestri_results, function(x) {
    data = lapply(x, function(file) {
      read.table(file, header = T, sep = ",")
    })
    names(data) = colnames
    return(data)
  })
  
  return(tapestri_results)
}

# filter the results 
filter_tapestri = function(x, filter_cells = TRUE, filter_ids = TRUE) {
  
  # select cells to keep 
  if(filter_cells) {
    good_barcodes = x$cells_info %>% 
      dplyr::filter(filtered == 0) %>% 
      dplyr::pull(barcode)
  } else {
    good_barcodes = x$cells_info %>% 
      dplyr::pull(barcode)
  }
  
  # select mutations to keep
  if(filter_ids) {
    good_variants = x$variants_info %>% 
      dplyr::filter(filtered == 0) %>% 
      dplyr::pull(id)
  } else {
    good_variants = x$variants_info %>% 
      dplyr::pull(id)
  }
  
  # select genotype matrix 
  # fix the names of the ids
  colnames(x$NGT) = replace_dots(colnames(x$NGT))
  x$NGT = x$NGT %>%
    tibble::column_to_rownames(var = "X")
  
  # select the ids and the cells 
  x$NGT = x$NGT[good_barcodes, good_variants]
  
  # select AF matrix 
  # fix the names of the ids
  colnames(x$AF) = replace_dots(colnames(x$AF))
  x$AF = x$AF %>%
    tibble::column_to_rownames(var = "X")
  
  # select the ids and the cells 
  x$AF = x$AF[good_barcodes, good_variants]
  
  # select DP matrix 
  # fix the names of the ids
  colnames(x$DP) = replace_dots(colnames(x$DP))
  x$DP = x$DP %>%
    tibble::column_to_rownames(var = "X")
  
  # select the ids and the cells 
  x$DP = x$DP[good_barcodes, good_variants]
  
  return(x)
  
}

# compute CCFs
compute_tapestri_CCF = function(x, thr_ccf = 0.1) {
  
  # calculate Tapestri CCF for each mutation
  filtered_ccfs = apply(x$NGT, 2, function(x) {
    
    # get the total number of cells 
    all_locus = length(x)
    
    # count for each locus which is the mutation state
    tot_muts = as.data.frame(table(x)) %>%
      dplyr::filter(x %in% c(1,2)) %>%
      dplyr::pull(Freq) %>%
      sum
    tot_muts/all_locus
  }) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("mut_id") %>%
    dplyr::rename(tapestri_ccf = ".") %>% 
    dplyr::as_tibble()
  
  # join new ccfs and the old ones and compare the bulk and the sc results wrt a defined threshold
  x$CCF_comparison = x$CCF_comparison %>%
    dplyr::mutate(X = NULL, Tapestri_id = NULL) %>%
    dplyr::full_join(., filtered_ccfs, by = join_by("mutation" == "mut_id")) %>%
    dplyr::select(mutation, starts_with("CCF"), tapestri_ccf) %>%
    dplyr::rename("Bulk_CCF" = starts_with("CCF")) %>%
    dplyr::rename("Tapestri_CCF" = tapestri_ccf) %>%
    dplyr::distinct() %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(diff = Bulk_CCF - Tapestri_CCF) %>% 
    dplyr::mutate(diff_abs = abs(diff)) %>% 
    dplyr::mutate(quality_ccf = 
                    case_when(diff_abs <= thr_ccf ~ "good", 
                              diff_abs > thr_ccf ~ "bad", 
                              is.na(diff_abs) ~ "not present in tapestri"))
  
  x$variants_info = x$variants_info %>% 
    dplyr::mutate(X = NULL)
  
  # add the infos on the variant from the sc experiment
  x$CCF_comparison = dplyr::full_join(x$CCF_comparison, x$variants_info, by = join_by("mutation" == "id")) 
  
  # add tapestri quality scores re depth
  x = add_quality_scores(x)
  
  return(x)
}

add_quality_scores = function(x){
  
  # compute mean depth of the locus and its sd
  dp_mean = apply(x$DP[,-1], 2, mean)
  dp_sd = apply(x$DP[,-1], 2, sd)
  
  dp_table = dplyr::tibble(id = names(dp_mean), 
                    DP_mean = dp_mean, 
                    DP_sd = dp_sd)
  
  # add it to the ccf information
  x$CCF_comparison = x$CCF_comparison %>% 
    dplyr::full_join(., dp_table, by = join_by("mutation" == "id")) %>% 
    dplyr::as_tibble()
  
  return(x)
}


get_patient_data = function(tapestri_results, patient, sample_t, bulk_data){
  
  # get the samples names for the specific patient at each time point
  tt = lapply(get_samples(x = tapestri_results, sample_t = sample_t), function(s){
    
    # add the the ccf table of tapestri the bulk information
    ccf = add_bulk_info(ccf = tapestri_results[[s]]$CCF_comparison, 
                        bulk_data = bulk_data, 
                        patient = patient)
    ccf = ccf %>% mutate(sample = s)
    
    # remove weird columns
    if("Bulk_AF" %in% colnames(ccf)) ccf = ccf %>% select(-Bulk_AF)
    ccf
  }) %>% do.call(rbind, .)
  
  return(tt)
}

add_bulk_info = function(ccf, bulk_data, patient){
  
  # select the data of one specific case and bind it to the tapestri data
  ex = bulk_data %>% 
    dplyr::filter(Case == patient) %>% 
    dplyr::distinct(mut_id, .keep_all = TRUE)
  
  x = dplyr::left_join(ccf, ex, by = join_by("mutation" == "mut_id"))
  return(x)
}

# get the samples (pos and neg) based on the beginning of the sample name
get_samples = function(x, sample_t){
  names(x) %>% grep(sample_t, ., value = TRUE)
}

get_patient_genotypes = function(tapestri_results, patient){
  
  # get the samples for the patient time point
  samples = get_samples(x = tapestri_results, sample_t = patient)
  
  # extract genotyping data and do some manipulation on them
  data = lapply(samples, function(s) {
    
    df = tapestri_results[[s]]$NGT %>% 
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ dplyr::case_when(
        . == 0 ~ 0, # wild type mutations stay 0
        . == 1 ~ 1, # mutated locus are counted in the same way
        . == 2 ~ 1,
        . == 3 ~ 0 # not consider the missing loci in the distance
      )))
    rownames(df) = paste(s, seq(1:nrow(df)), sep = "_")
    
    # transpose them
    df %>% 
      t %>%
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "mut") 
  }) %>% 
    Reduce(dplyr::full_join, .) 
  
  data = tibble::column_to_rownames(data, var = "mut")
  
  # remove eventual NAs due to the join and set them to 0 
  data = data %>% 
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ dplyr::case_when(
      is.na(.) ~ 0, 
      .default = .
    )))
  
  # convert data to matrix
  data = as.matrix(data)
  
  data
}

generate_heatmap = function(patient_name, 
                            patient_data, 
                            patient_genotypes, 
                            clusters_legend, 
                            ht_cols, 
                            samples_cols, 
                            path_ht) {
  
  get_ordered_mutations = function(x, clusters_legend){
    
    mutations = x %>% 
      filter(filtered == 0) %>% 
      select(viber_cluster, mutation) %>% 
      distinct() %>% 
      mutate(viber_cluster = factor(viber_cluster, levels = clusters_legend$cluster))
    
    mut_order = setNames(
      nm = mutations %>% arrange(viber_cluster) %>% pull(mutation), 
      object = mutations %>% arrange(viber_cluster) %>% pull(viber_cluster)
    )
    
    mut_order
    
  }
  
  mut_order = get_ordered_mutations(patient_data, clusters_legend)
  patient_genotypes = patient_genotypes[names(mut_order), ]
  ann_cols = setNames(clusters_legend$colors, clusters_legend$cluster)
  
  ha = rowAnnotation(
    viber_cluster = unname(mut_order),
    col = list(viber_cluster = ann_cols)
  )
  
  all_samples = gsub("_[0-9]{1,5}$", "", colnames(patient_genotypes))
  samples = patient_data$sample %>% unique() %>% sort()
  sample_cols = setNames(samples_cols, samples)
  
  sample_ann = HeatmapAnnotation(sample = all_samples, col = list(sample = sample_cols))
  
  ht = ComplexHeatmap::Heatmap(patient_genotypes, 
                               name = "Genotype", 
                               col = cols, 
                               row_names_gp = gpar(size = 3), 
                               cluster_columns = T,
                               column_dend_height = unit(4, "cm"),
                               show_column_dend = T,
                               cluster_rows = T, 
                               show_row_dend = F,
                               left_annotation = ha,
                               show_row_names = F, 
                               show_column_names = F, 
                               heatmap_legend_param = 
                                 list(at = c("0", "1"), labels = c("Wild Type/Missing", "Mutant (hom-het)")),
                               row_split = unname(mut_order), 
                               cluster_row_slices = F, 
                               # column_split = factor(all_samples, levels = all_samples %>% unique), 
                               top_annotation = sample_ann)
  
  pdf(paste0(path_ht, "/", patient_name, "_heatmap_genotypes.pdf"), width = 8, height = 10)
  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  graphics.off()
}
