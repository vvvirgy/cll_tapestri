# colors of the clusters 
# doing some clean up to fix some errors with the clusters (manually done)
colors_clusters_by_patient = list(
  
  "CT339_T2" = tibble(
    cluster = factor(c('Truncal', 'S2', 'S3'), 
                     levels = c('Truncal', 'S2', 'S3')), 
    colors = c('steelblue', 'goldenrod', 'purple3')
  ), 
  
  "CT344_T2" = tibble(
    cluster = factor(c('Truncal', 'S2', 'S3', 'S3_T2+'), 
                     levels = c('Truncal', 'S2', 'S3', 'S3_T2+')), 
    colors = c('steelblue', 'goldenrod', 'purple3', 'forestgreen')
  ),
  
  'CT48_T2' = tibble(
    cluster = factor(c('Truncal', 'S1', 'S3', 'S2'), 
                     levels = c('Truncal', 'S1', 'S3', 'S2')), 
    colors = c('steelblue', 'forestgreen', 'purple3', 'goldenrod')
  ),
  
  "CT525_T2" = tibble(
    cluster = factor(c('Truncal', 'S1', 'S2'), 
                     levels = c('Truncal', 'S1', 'S2')), 
    colors = c('steelblue', 'goldenrod', 'purple3')
  ),
  
  "CZ40_T2" = tibble(
    cluster = factor(c('Truncal', 'S3', 'S2', 'H2'), 
                     levels = c('Truncal', 'S3', 'S2', 'H2')), 
    colors = c('steelblue', 'purple3', 'goldenrod', 'forestgreen')
  ),
  
  "RM238_T2" = tibble(
    cluster = factor(c('Truncal', 'S1', 'S3'), 
                     levels = c('Truncal', 'S1', 'S3')), 
    colors = c('steelblue', 'purple3', 'forestgreen')
  ),
  
  "TS187_T2" = tibble(
    cluster = factor(c('Truncal', 'S2', 'S1', 'S5'), 
                     levels = c('Truncal', 'S2', 'S1', 'S5')), 
    colors = c('steelblue','forestgreen', 'goldenrod','purple3')
  ), 
  
  "TS187_T3" = tibble(
    cluster = factor(c('Truncal', 'S2', 'S1', 'S5'), 
                     levels = c('Truncal', 'S2', 'S1', 'S5')), 
    colors = c('steelblue','forestgreen', 'goldenrod','purple3')
  )
)

# heatmap cells colors
cols = setNames(object = c("#F8FAFC", "#40A578"), c(0,1))

samples_cols = c("#82CCE0", "#CE6D6E")
