library(tidyverse)

# library(gt)
# library(dendextend)

# loading data obtained from tapestri
tapestri_data = load_tapestri_data(path = tapestri_data_path)

# filter the results to get only good mutations and cells, then recalculate the CCFs
tapestri_filtered = lapply(tapestri_data, function(x) {filter_tapestri(x, filter_cells = T, filter_ids = T)})

# compute CCFs after the filtering -- kept only good quality cells and mutations
tapestri_filtered = lapply(tapestri_filtered, compute_tapestri_CCF)

# load the bulk information -- clusters assignment and ccf
bulk_data = readxl::read_excel(bulk_data_path, sheet = 10)

# removing useless columns and fixing mut ids
bulk_data = bulk_data %>% 
  dplyr::mutate(nd...1 = NULL, nd...41 = NULL) %>% 
  dplyr::mutate(mut_id = paste0(Chr, ":", Start, ":", Ref, "/", Var))

# loading corrected clusters for a patient
cz40_bulk = readRDS(cz40_bulk_path)
cz40_bulk = cz40_bulk %>% 
  tidyr::separate(id, into = c('chr', 'from', 'to', 'ref', 'alt'), sep = ":") %>% 
  dplyr::mutate(id = paste0(chr, ":", from, ":", ref, "/", alt)) %>% 
  dplyr::select(id, viber_cluster, phylo_color) %>% 
  dplyr::mutate(Case = "CZ40") %>% 
  dplyr::rename(new_viber_cluster = viber_cluster)

bulk_data = bulk_data %>% 
  split(.$Case)

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
