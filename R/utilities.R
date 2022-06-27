
#' Convert a dataframe to a matrix. First column becomes row names
#' @param mydataframe input dataframe
df2m <- function(mydataframe){
  mymatrix <- mydataframe
  rownames(mymatrix) <- mydataframe[,1]
  mymatrix[,1] <- NULL
  mymatrix <- as.matrix(mymatrix)
  mymatrix
}

#' Convert a matrix to a dataframe. Rows become first column.
#' @param column_name Name of the new column (defaults to "Row_variable")
m2df <- function(mymatrix, column_name = "Row_variable"){
  mydf <- as.data.frame(mymatrix)
  cur_names <- names(mydf)
  mydf[, column_name] <- rownames(mydf)
  rownames(mydf) <- NULL
  mydf <- mydf[,c(column_name,cur_names)]
  return(mydf)
}

#' Convert a dataframe to a list. Only first two columns are used.
#' @param mydataframe input dataframe
df2l <- function(mydataframe){
  setNames(mydataframe[,2], mydataframe[,1])
}

#' Convenience function to remove item from a vector
#' @param input vector
#' @param to_remove items to remove
remove_from_vector <- function(vector, to_remove){
  vector[!vector %in% to_remove]
}

#' Remove inf and na values from a matrix
#' @param mymatrix input matrix
remove_inf_na_matrix <- function(mymatrix){
  do.call(data.frame, lapply(m2df(mymatrix), function(x) {
    replace(x, is.infinite(x) | is.na(x), 0)
  })
  ) %>% df2m
}

#' Make row names of a dataframe the value of a defined column (default 1st) and then remove the column
#' @param mydataframe input dataframe
#' @param rowname_col column to get row names from (default = 1)
df_column_to_rownames <- function(mydataframe, rowname_col = 1){
  my_clean.df <- mydf
  rownames(my_clean.df) <- my_clean.df[,rowname_col]
  my_clean.df[,1] <- NULL
  return(my_clean.df)
}

#' Transfer factor levels from a dataframe to another dataframe or list.
#' @param mydataframe input dataframe
#' @param recieving dataframe or list receiving factor levels
transfer_factor_levels <- function(original_dataframe, recieving){
  for(cn in colnames(original_dataframe)) {
    if (cn %in% names(recieving) && is.factor(original_dataframe[,cn])){
      if (!inherits(recieving, "data.frame")){
        recieving[[cn]] <- recieving[[cn]][levels(original_dataframe[,cn])]
      } else if (inherits(recieving, "data.frame")){
        recieving[,cn] <- factor(recieving[,cn], levels = levels(original_dataframe[,cn]))
      }
    }
  }
  recieving
}

#' The geometric mean, with some error-protection bits.
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}

#' Center log ratio transform
clr <- function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}


read_counts_and_unique_features <- function(count_matrix, sample_list){
  # Number of features and read counts
  temp <- count_matrix
  temp[temp > 0] <- 1
  read_sums <- colSums(count_matrix)
  feature_sums <- colSums(temp)
  sample_read_counts <- data.frame(sample_list, unlist(lapply(sample_list,
                                                              function(x) ifelse(x %in% names(read_sums), read_sums[x],0))))
  sample_feature_counts <- data.frame(sample_list, unlist(lapply(sample_list,
                                                                 function(x) ifelse(x %in% names(feature_sums), feature_sums[x],0))))
  out <- list("sample_read_counts" = sample_read_counts,
              "sample_feature_counts" = sample_feature_counts)
  out
}





#' Create p-value label for input value
#' @param x integer/flout to create label for
#' @return p-value label string, e.g. "***"
get_pvalue_labels <- function(x){
  ifelse(x <= 0.001, "***", ifelse(x <= 0.01, "**", ifelse(x <= 0.05, "*", "NS")))
}

#' Get prevalence of rows in matrix across columns
#' @param matrix.m input matrix
get_prevalence <- function(matrix.m){
  apply(matrix.m, 1, function(x) {length(which(x > 0)) / length(x)})
}

#' Get the lowest taxonomy level (or 'first') that has been resolved from a taxonomy string
#' @param x taxonomy string
#' @examples
#' first_resolved_taxonomy("d__Bacteria;p__Acidobacteriota;Unassigned")
#' first_resolved_taxonomy("d__Bacteria;p__Acidobacteriota;uncultured")
first_resolved_taxonomy <- function(x) {
  ranks.v <- unlist(strsplit(x, split = ';'))
  for (i in length(ranks.v):1) {
    split.v <- unlist(strsplit(ranks.v[i], split = '__'))
    if (!is.na(split.v[2]) & split.v[2] != "") {
      if (!grepl(split.v[2], pattern = "Unassigned|uncultured")) {
        return(ranks.v[i])
      }
    }
  }
  return(ranks.v[1])
}

#' Customise a taxonomy string. Returns specified taxonomy levels. Assumes standard 7 ranks.
#' @param x taxonomy string
#' @param taxa_levels vector of taxa levels to return.
#' @examples customise_taxonomy_string("d__Bacteria;p__Acidobacteriota", taxa_levels = c("phylum"))
customise_taxonomy_string <- function(x,
                                      taxa_levels = c("domain",
                                                      "phylum",
                                                      "class",
                                                      "order",
                                                      "family",
                                                      "genus",
                                                      "species")){
  if (grepl(";", x)){
    ranks.v <- unlist(strsplit(x, split = ';'))
    domain <- ""
    phylum <- ""
    class <- ""
    order <- ""
    family <- ""
    genus <- ""
    species <- ""
    if (length(ranks.v) == 7){
      domain <- ranks.v[1]
      phylum <- ranks.v[2]
      class <- ranks.v[3]
      order <- ranks.v[4]
      family <- ranks.v[5]
      genus <- ranks.v[6]
      species <- ranks.v[7]
    } else if (length(ranks.v) == 6){
      domain <- ranks.v[1]
      phylum <- ranks.v[2]
      class <- ranks.v[3]
      order <- ranks.v[4]
      family <- ranks.v[5]
      genus <- ranks.v[6]
    } else if (length(ranks.v) == 5){
      domain <- ranks.v[1]
      phylum <- ranks.v[2]
      class <- ranks.v[3]
      order <- ranks.v[4]
      family <- ranks.v[5]
    } else if (length(ranks.v) == 4){
      domain <- ranks.v[1]
      phylum <- ranks.v[2]
      class <- ranks.v[3]
      order <- ranks.v[4]
    } else if (length(ranks.v) == 3){
      domain <- ranks.v[1]
      phylum <- ranks.v[2]
      class <- ranks.v[3]
    } else if (length(ranks.v) == 2){
      domain <- ranks.v[1]
      phylum <- ranks.v[2]
    } else if (length(ranks.v) == 1){
      domain <- ranks.v[1]
    }
    taxa_level_strings.l <-
      c("domain" = domain,
        "phylum" = phylum,
        "class" = class,
        "order" = order,
        "family" = family,
        "genus" = genus,
        "species" = species)
    # print(paste0(as.character(lapply(taxa_levels, function(x) taxa_level_strings.l[x])),collapse = ";"))
    paste0(as.character(lapply(taxa_levels, function(x) get(x))),collapse = ";")
  } else{
    as.character(x)
  }
}

customise_taxonomy_string_vectorized <- function(x, taxa_levels){
  unlist(lapply(x,customise_taxonomy_string, taxa_levels))
}



generate_tax_level_data <- function(feature_count_taxonomy_data.df, sample_ids, tax_string_levels, remove_zero_row_entries = F){
  # This function will take a dataframe with the feature counts (first column ID and remaining columns the sample counts and taxonomy labels),
  # a list of the samples and the taxa labels you wish to sum over
  # The output is a list with the counts and relative abundances at each specified taxonomy level

  df2matrix <- function(mydataframe){
    mymatrix <- mydataframe
    rownames(mymatrix) <- mydataframe[,1]
    mymatrix[,1] <- NULL
    mymatrix <- as.matrix(mymatrix)
    mymatrix
  }

  m2df <- function(mymatrix, column_name = "Row_variable"){
    mydf <- as.data.frame(mymatrix)
    cur_names <- names(mydf)
    mydf[, column_name] <- rownames(mydf)
    rownames(mydf) <- NULL
    mydf <- mydf[,c(column_name,cur_names)]
    return(mydf)
  }
  output <- list()
  for (tax_string_level in tax_string_levels){
    counts.df <- melt(feature_count_taxonomy_data.df[c(tax_string_level, sample_ids)], measure.vars = c(sample_ids))
    counts.df <- dcast(counts.df, get(tax_string_level) ~ variable, sum)
    names(counts.df)[1] <- tax_string_level

    # Now create the relative abundance matrix at the current taxa level
    abundances.m <- counts.df

    rownames(abundances.m) <- abundances.m[[tax_string_level]]
    abundances.m[tax_string_level] <- NULL
    abundances.m <- as.matrix(abundances.m)
    abundances.m <- t(t(abundances.m) / colSums(abundances.m))
    abundances.m[is.nan(abundances.m)] <- 0

    counts.m <- df2matrix(counts.df)
    if (remove_zero_row_entries == T){
      abundances.m <- abundances.m[apply(abundances.m, 1, max) > 0,]
      counts.m <- counts.m[apply(counts.m, 1, max) > 0,]
    }
    output[[tax_string_level]] <- list("counts" = counts.m, "abundances" = abundances.m)
  }
  output
}

abundance_summary <- function(abundance.df, top_n = NULL){
  # Generate summary from abundance dataframe.
  # Assumes first column is the unique taxa/feature ID and other columns are abundances for each sample
  # If top_n specified, will return the top n taxa by mean abundance
  abundance_melt.df <- melt(abundance.df, variable.name = "Sample", value.name = "Relative_abundance")
  abundance_summary.df <-
    abundance_melt.df %>%
    dplyr::mutate(Samples_total = n_distinct(Sample)) %>% # number of unique samples/index
    dplyr::group_by_at(1,.drop = F) %>%
    dplyr::mutate(Samples_present = sum(Relative_abundance > 0, na.rm = TRUE)) %>%
    dplyr::summarise(
      Samples_present = max(Samples_present),
      Samples_total = max(Samples_total),
      # Mean_relative_abundance2 = round(mean(sum(Relative_abundance)/max(Samples_total)), 5),
      Min_relative_abundance = round(min(Relative_abundance),5),
      Max_relative_abundance = round(max(Relative_abundance),5),
      Mean_relative_abundance = round(mean(Relative_abundance), 5),
      Median_relative_abundance = round(median(Relative_abundance), 5),
    ) %>%
    arrange(dplyr::desc(Mean_relative_abundance))
  if (!is.null(top_n)){
    abundance_summary.df <-
      abundance_summary.df %>%
      arrange(dplyr::desc(Mean_relative_abundance)) %>%
      head(top_n)
  }
  as.data.frame(abundance_summary.df)
}


abundance_summary_with_metadata <- function(abundance_data.df,
                                            metadata.df,
                                            variable_to_summarise,
                                            metadata_sample_column = "Sample",
                                            n_taxa = 10) {
  # Summarise taxa abundances within a group/variable
  # Assumes first column in abundance_data.df is the unique taxa/feature ID and other columns are abundances for each sample

  collapse_abundance_dataframe <- function(abundance.df, filtering_taxa.v = NULL){
    # Collapse abundance dataframe.
    # Assumes first column is the unique taxa/feature ID and other columns are abundances for each sample
    # If a list of filtering_taxa.v is provided, taxa not in list are assigned a label 'Other' and their abundance
    # values summed together
    abundance.df <- melt(abundance.df, variable.name = "Sample", value.name = "Relative_abundance")
    if (!is.null(filtering_taxa.v)){
      abundance.df[!abundance.df[,1] %in% filtering_taxa.v,][,1] <- "Other"
    }

    abundance.df <-
      abundance.df %>%
      group_by_at(c(2,1)) %>%
      summarise_all(list(sum)) %>%
      # summarise_all(funs(sum)) %>%
      as.data.frame()
    abundance.df
  }

  abundances_summary.df <-
    collapse_abundance_dataframe(abundance_data.df) %>%
    left_join(., metadata.df, by = c("Sample" = metadata_sample_column)) %>% # combine with metadata
    # dplyr::group_by_at(vars(variable_to_summarise,dplyr::starts_with("taxonomy_"))) %>%
    dplyr::group_by_at(vars(variable_to_summarise,2)) %>%
    dplyr::summarise(
      Variable = variable_to_summarise,
      Mean_relative_abundance = mean(Relative_abundance),
      Median_relative_abundance = median(Relative_abundance)) %>%
    dplyr::group_by_at(variable_to_summarise) %>%
    dplyr::arrange(.by_group = T, dplyr::desc(Mean_relative_abundance)) %>%
    top_n(n_taxa,wt = Mean_relative_abundance) %>%
    relocate(Variable) %>%
    dplyr::rename(., Group = all_of(variable_to_summarise)) %>%
    as.data.frame()

  abundances_summary.df
}


collapse_abundance_dataframe <- function(abundance.df, filtering_taxa.v = NULL){
  # Collapse abundance dataframe.
  # Assumes first column is the unique taxa/feature ID and other columns are abundances for each sample
  # If a list of filtering_taxa.v is provided, taxa not in list are assigned a label 'Other' and their abundance
  # values summed together
  abundance.df <- melt(abundance.df, variable.name = "Sample", value.name = "Relative_abundance")
  if (!is.null(filtering_taxa.v)){
    if (!all(abundance.df[,1] %in% filtering_taxa.v)){
      abundance.df[!abundance.df[,1] %in% filtering_taxa.v,][,1] <- "Other"
    }
  }

  abundance.df <-
    abundance.df %>%
    group_by_at(c(2,1)) %>%
    summarise_all(list(sum)) %>%
    # summarise_all(funs(sum)) %>%
    as.data.frame()
  abundance.df
}



calculate_PC_taxa_contributions <- function(pca_object){
  # Calculate the percentage contribution from each taxa for PC1-3. Requires unscaled values that are squared
  pc1_contribution <- melt(round(100*scores(pca_object, display = "species", scaling = 0)[,1]^2, 3),value.name = "PC1_contribution_percentage")
  pc2_contribution <- melt(round(100*scores(pca_object, display = "species", scaling = 0)[,2]^2, 3),value.name = "PC2_contribution_percentage")
  pc3_contribution <- melt(round(100*scores(pca_object, display = "species", scaling = 0)[,2]^2, 3),value.name = "PC3_contribution_percentage")

  data.frame(pc1_contribution, pc2_contribution, pc3_contribution)
}
