
run_permanova_custom <- function(my_metadata, my_formula, my_method = "euclidean", permutations = 999, label = NULL){
  stat_sig_table <- NULL
  result <- adonis(my_formula,data = my_metadata, permu=permutations,method= my_method)
  # result <- adonis(my_formula,data = my_metadata, permu=999,method="bray")
  for (r in rownames(result$aov.tab)){
    variable <- r
    Degress_of_freedom <- result$aov.tab[r,]$Df[1]
    SumOfSqs <- round(result$aov.tab[r,]$SumsOfSqs[1], 3)
    meanSqs <- round(result$aov.tab[r,]$MeanSqs[1], 3)
    F.model <- round(result$aov.tab[r,]$F.Model[1], 3)
    R2 <- round(result$aov.tab[r,]$R2[1], 3)
    p_value <- round(result$aov.tab[r,]$`Pr(>F)`[1], 5)
    stat_sig_table <- rbind(stat_sig_table, data.frame(variable,
                                                       Degress_of_freedom,
                                                       SumOfSqs,
                                                       meanSqs,
                                                       F.model,
                                                       R2,
                                                       p_value))
  }
  # my_formula_string <- paste0(as.character(my_formula)[2], as.character(my_formula)[1], as.character(my_formula)[3])
  my_formula_string <- paste0(as.character(my_formula)[1], as.character(my_formula)[3])
  print(paste0("FORMULA: ", my_formula_string))
  print(result)
  names(stat_sig_table) <- c("Term","Df", "SumOfSqs","MeanSqs","F.Model","R2","Pr(>F)")
  stat_sig_table <- stat_sig_table[order(stat_sig_table$"Pr(>F)"),]
  stat_sig_table$Method <- my_method
  stat_sig_table$Formula <- my_formula_string
  if (!is.null(label)){
    stat_sig_table$Label <- label
  }
  stat_sig_table
}



run_permdisp_custom <- function(my_metadata, my_data, my_group, my_method = "euclidean", permutations = 999, label = NULL){
  stat_sig_table <- NULL
  dist_matrix <- vegdist(t(my_data), method = my_method)
  betadisper_object <- with(my_metadata, betadisper(dist_matrix, group = get(my_group)))
  permutest_results <- permutest(betadisper_object, permutations = permutations, parallel = 2)

  for (r in rownames(permutest_results$tab)){
    variable <- r
    Degrees_of_freedom <- permutest_results$tab[r,]$Df[1]
    SumOfSqs <- round(permutest_results$tab[r,]$`Sum Sq`[1],3)
    meanSqs <- round(permutest_results$tab[r,]$`Mean Sq`[1], 3)
    F.model <- round(permutest_results$tab[r,]$F[1], 3)
    N_permutations <- permutest_results$tab[r,]$N.Perm[1]
    p_value <- round(permutest_results$tab[r,]$`Pr(>F)`[1], 5)
    stat_sig_table <- rbind(stat_sig_table, data.frame(variable,
                                                       Degrees_of_freedom,
                                                       SumOfSqs,
                                                       meanSqs,
                                                       F.model,
                                                       N_permutations,
                                                       p_value))
  }
  print(permutest_results)
  names(stat_sig_table) <- c("Term","Df", "SumOfSqs","MeanSqs","F.Model","Permutations","Pr(>F)")
  stat_sig_table <- stat_sig_table[order(stat_sig_table$"Pr(>F)"),]
  stat_sig_table$Method <- my_method
  stat_sig_table$Group <- my_group
  if (!is.null(label)){
    stat_sig_table$Label <- label
  }
  stat_sig_table
}


run_anosim_custom <- function(my_metadata, my_data, my_group, my_method = "euclidean", permutations = 999, label = NULL){
  stat_sig_table <- NULL
  anosim_object <- with(my_metadata, anosim(x = t(my_data),
                                            grouping = get(my_group),
                                            permutations = permutations,
                                            distance = my_method,
                                            parallel = 2))

  stat_sig_table <- rbind(stat_sig_table, data.frame(my_group,
                                                     anosim_object$statistic,
                                                     anosim_object$signif,
                                                     anosim_object$permutations))
  names(stat_sig_table) <- c("Variable","R_statistic", "Significance","Permutations")
  stat_sig_table$Method <- my_method
  if (!is.null(label)){
    stat_sig_table$Label <- label
  }
  stat_sig_table
}

