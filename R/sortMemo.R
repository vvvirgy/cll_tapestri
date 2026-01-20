memoSort <- function(M) {
  geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol);
  sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
  return(M[geneOrder, sampleOrder]);
}

data_bis = tapestri_results[[x]]$NGT[,-1] %>% 
  # mutate(across(everything(), ~ case_when(
  #   . == 3 ~ 0, 
  #   .default = .
  # ))) %>%
  as.matrix() %>% t

rownames(data_bis) = replace_dots(rownames(data_bis))


oncoprint_column_order = function() {
  scoreCol = function(x) {
    score = 0
    for(i in 1:length(x)) {
      if(x[i]) {
        score = score + 2^(length(x)-i*1/x[i])
      }
    }
    return(score)
  }
  scores = apply(count_matrix[row_order, ,drop = FALSE], 2, scoreCol)
  order(scores, decreasing=TRUE)
}

