
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

run_permanova_custom2 <- function(metadata.df,
                                  my_formula,
                                  method = "euclidean",
                                  by = "terms",
                                  permutations = 999,
                                  strata_variable = NULL,
                                  label = NULL,
                                  verbose = T
){
  my_formula <- as.formula(my_formula)
  stat_sig_table <- NULL
  strata <- NULL
  if (!is.null(strata_variable)){
    strata <- metadata.df[,strata_variable]
  }
  result <- adonis2(my_formula,data = metadata.df, permu=permutations,method= method, by = by, strata = strata)
  # result <- adonis(my_formula,data = my_metadata, permu=999,method="bray")
  result <- as.data.frame(result)
  for (r in rownames(result)){
    variable <- r
    Degress_of_freedom <- result[r,]$Df[1]
    SumOfSqs <- round(result[r,]$SumOfSqs[1], 3)
    # meanSqs <- round(result[r,]$MeanSqs[1], 3)
    F.model <- round(result[r,][,"F"][1], 3)
    R2 <- round(result[r,]$R2[1], 3)
    p_value <- round(result[r,]$`Pr(>F)`[1], 5)
    stat_sig_table <- rbind(stat_sig_table, data.frame(variable,
                                                       Degress_of_freedom,
                                                       SumOfSqs,
                                                       # meanSqs,
                                                       F.model,
                                                       R2,
                                                       p_value))
  }
  # my_formula_string <- paste0(as.character(my_formula)[2], as.character(my_formula)[1], as.character(my_formula)[3])
  formula_string <- paste0(as.character(my_formula)[1], as.character(my_formula)[3])
  if (verbose) {
    print(paste0("FORMULA: ", formula_string))
    print(result)
  }
  names(stat_sig_table) <- c("Term","Df", "SumOfSqs","F.Model","R2","Pr(>F)")
  stat_sig_table <- stat_sig_table[order(stat_sig_table$"Pr(>F)"),]
  stat_sig_table$Method <- method
  stat_sig_table$Formula <- formula_string
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

run_permdisp_custom2 <- function(metadata.df, input_data, group, method = "euclidean", permutations = 999, label = NULL, verbose = T){
  stat_sig_table <- NULL
  if (!inherits(input_data, "dist")){
    message("Input data is of not class of 'dist', generating distance object")
    input_data <- vegdist(input_data, method = method) # assume sites (samples) are rows
  } else{
    method = attr(input_data,which = "method")
  }

  # Remove all factor levels
  metadata.df[] <- lapply(metadata.df, function(x) {
    if (is.factor(x)) {
      as.character(x)
    } else {
      x
    }
  })

  betadisper_object <- with(metadata.df, betadisper(input_data, group = get(group)))
  permutest_results <- permutest(betadisper_object, permutations = permutations, parallel = 2, pairwise = T)

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
  if (verbose){
    print(permutest_results)
  }
  names(stat_sig_table) <- c("Term","Df", "SumOfSqs","MeanSqs","F.Model","Permutations","Pr(>F)")
  stat_sig_table <- stat_sig_table[order(stat_sig_table$"Pr(>F)"),]
  stat_sig_table$Method <- method
  stat_sig_table$Group <- group

  pairwise_results.df <- data.frame()
  for (i in names(permutest_results$pairwise$permuted)){
    groups.v <- unlist(strsplit(i, split = "-"))
    value <- permutest_results$pairwise$permuted[[i]]
    pairwise_results.df <- rbind(pairwise_results.df, data.frame("Group_1" = groups.v[1], "Group_2" = groups.v[2], "P_value" = value))
  }
  pairwise_results.df$Method <- method
  pairwise_results.df$Group <- group

  if (!is.null(label)){
    stat_sig_table$Label <- label
    pairwise_results.df$Label <- label
  }
  list("main_result" = stat_sig_table, "pairwise_results" = pairwise_results.df)
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

