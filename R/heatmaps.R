plot_heatmap_v2 <- function(input_matrix.m,
                                   column_annotation_data.df = NULL,
                                   row_annotation_data.df = NULL,

                                   column_discrete_variables_to_plot.v = NULL,
                                   column_continuous_variables_to_plot.v = NULL,
                                   row_discrete_variables_to_plot.v = NULL,
                                   row_continuous_variables_to_plot.v = NULL,

                                   column_continuous_variable_breaks.l = NULL,
                                   row_continuous_variable_breaks.l = NULL,

                                   output_filename = NULL,
                                   height = 10, width = 10,
                                   row_title = "Gene",
                                   column_title = "Sample",
                                   heatmap_palette_choice = 1,

                                   annotation_name_size = 4,
                                   annotation_legend_label_size = 7,
                                   annotation_legend_title_size = 8,
                                   annotation_legend_length = 2,
                                   annotation_legend_height = 3,

                                   column_split_variable = NULL,
                                   column_label_variable = NULL,
                                   column_labels.df = NULL,

                                   row_split_variable = NULL,
                                   row_label_variable = NULL,
                                   row_labels.df = NULL,

                                   column_title_font_size = 7,
                                   column_labels_font_size = 6,
                                   row_title_font_size = 7,
                                   row_labels_font_size = 6,

                                   column_annotation_padding = 0.2,
                                   row_annotation_padding = 0.2,

                                   cluster_rows = T,
                                   cluster_columns = T,
                                   show_row_dend = T,
                                   show_column_dend = T,

                                   heatmap_legend_padding = 2,
                                   plot_padding = c(1,1,1,1),# bottom, left, top and right

                                   main_legend_title = NULL,

                                   order_heatmap_columns_by_labels = F, # Order the heatmap rows by the column labels
                                   order_heatmap_rows_by_labels = F, # Order the heatmap rows by the row labels
                                   show_cell_values = F,
                                   # If show_cell_values =T, cells less than this will have a black font colour
                                   # and above white
                                   cell_fun_value_col_threshold = 15,
                                   my_cell_fun = NULL # function to apply to the cell values
){
  # if (!is.null(taxa_data.df)){
  #   # Filter taxa_data.df to entries in input_matrix.m.
  #   # Requires first column to be rowname
  #   taxa_data.df %>% dplyr::filter(cur_data()[[1]] %in% rownames(input_matrix.m))
  #
  #   # Ensure order of taxa_data.df matches input matrix.
  #   taxa_data.df <- taxa_data.df %>% dplyr::arrange(factor(cur_data()[[1]], levels = rownames(input_matrix.m)))
  #
  # }
  # Create row and column labels
  get_row_labels <- function(x, labeller_function){
    data.frame(rownames(x), unlist(lapply(rownames(x), labeller_function)))
  }
  # column_labels.df <- metadata.df[,c("Index", "Sample_name")]

  # ---------------------------------------------------------------------------------------------------------
  # ---------------------------------------------------------------------------------------------------------
  # Create colour palettes

  # Choice of heatmap palette
  if (heatmap_palette_choice == 1){
    colours.v <- c("white", rev(viridis(8)))
  } else if (heatmap_palette_choice == 2){
    colours.v <- c("white", rev(magma(8)))
  } else if (heatmap_palette_choice == 3){
    colours.v <-  rev(cividis(8))
    colours.v <- c(colours.v[1:4], "white", colours.v[5:8])
  } else{
    error_message <- paste0("Abundance palette choice ", heatmap_palette_choice,
                            "is invalid. Choose from '1','2' or '3'.")
    stop(error_message)
  }

  abundance_palette <- colorRampPalette(colours.v)
  # TODO dynamic breaks and/or user provided
  abundance_breaks.v <- unlist(lapply(c(0.0001, 0.1, 0.5,5, seq(10,40,10)),log,2))
  # abundance_breaks.v <- unlist(c(0, 0.1, 0.5,5, seq(10,20,10)))
  abundance_col_fun <- circlize::colorRamp2(breaks = abundance_breaks.v, colors = abundance_palette(length(abundance_breaks.v)))

  # ---------------------------------------------------------------------------------------------------------
  # Create continuous variable palettes. Only for those variables that have defined breaks.
  # TODO Add ability to provide colours, currently random
  # TODO Reduce duplicate code
  column_continuous_variables_palettes.l <- list()
  for (variable in column_continuous_variables_to_plot.v){
    if (variable %in% names(column_continuous_variable_breaks.l)){
      column_continuous_variables_palettes.l[[variable]] <- make_continuous_palette(column_annotation_data.df, variable,
                                                                                    annotation_palette = colorRampPalette(c("white", circlize::rand_color(1))),
                                                                                    custom_breaks.v = column_continuous_variable_breaks.l[[variable]])
    }
  }
  row_continuous_variables_palettes.l <- list()
  for (variable in row_continuous_variables_to_plot.v){
    if (variable %in% names(row_continuous_variable_breaks.l)){
      row_continuous_variables_palettes.l[[variable]] <- make_continuous_palette(row_annotation_data.df, variable,
                                                                                 annotation_palette = colorRampPalette(c("white", circlize::rand_color(1))),
                                                                                 custom_breaks.v = row_continuous_variable_breaks.l[[variable]])
    }
  }
  # ---------------------------------------------------------------------------------------------------------
  # Create colour palettes for discrete variables
  # TODO Reduce duplicate code
  column_discrete_variables_palettes.l <- list()
  for (variable in column_discrete_variables_to_plot.v){
    if (paste0(variable, "_colour") %in% colnames(column_annotation_data.df)){
      column_discrete_variables_palettes.l[[variable]] <- colour_list_from_dataframe(column_annotation_data.df, variable)
    }
  }
  row_discrete_variables_palettes.l <- list()
  for (variable in row_discrete_variables_to_plot.v){
    if (paste0(variable, "_colour") %in% colnames(row_annotation_data.df)){
      row_discrete_variables_palettes.l[[variable]] <- colour_list_from_dataframe(row_annotation_data.df, variable)
    }
  }

  # ---------------------------------------------------------------------------------------------------------
  # ---------------------------------------------------------------------------------------------------------
  # Create annotations for variables
  column_annotation <- NULL
  row_annotation <- NULL

  # Create annotation for column variables
  all_column_variables.v <- c(column_discrete_variables_to_plot.v, column_continuous_variables_to_plot.v)
  if (length(all_column_variables.v) > 0){
    all_column_variables.v <- all_column_variables.v[all_column_variables.v %in% colnames(column_annotation_data.df)]

    # Filter column annotation data to those entries present in the heatmap matrix
    # Assumes rownames defined
    column_annotation_data.df <- column_annotation_data.df[rownames(column_annotation_data.df) %in% colnames(input_matrix.m),,drop =F]

    # Order/filter the heatmap matrix to order/entries of metadata
    input_matrix.m <- input_matrix.m[,rownames(column_annotation_data.df),drop = F]

    # Order the heatmap matrix by the variables
    input_matrix.m <- input_matrix.m[,do.call(order, column_annotation_data.df[,all_column_variables.v,drop=F]),drop =F]

    # Order the column_annotation_data.df by the variables
    column_annotation_data.df <- column_annotation_data.df[do.call(order, column_annotation_data.df[,all_column_variables.v,drop=F]),,drop=F]

    # Check that rownames match colnames
    if (!all(rownames(column_annotation_data.df) == colnames(input_matrix.m))){
      stop("Row names in column annotation data do not match column names in matrix")
    }

    # Create final annotation dataframe from the column annotation data
    annotation_data.df <- column_annotation_data.df[,all_column_variables.v,drop = F]

    if (length(column_continuous_variables_to_plot.v) > 0){
      # Remove variables/columns that are all zero/NA
      annotation_data.df <-
        annotation_data.df %>%
        dplyr::select(-names(which(apply(annotation_data.df[,column_continuous_variables_to_plot.v,drop =F], 2, sum, na.rm = T) == 0)))
      column_continuous_variables_to_plot.v <- column_continuous_variables_to_plot.v[column_continuous_variables_to_plot.v %in% names(annotation_data.df)]
    }

    # Build the annotation_legend_param list for each variable
    annotation_legend_param.l <- list()
    for (variable in names(annotation_data.df)){
      if (variable %in% names(column_continuous_variable_breaks.l)){

        annotation_legend_param.l[[variable]] <-
          list(direction = "horizontal",border =T,
               labels_gp = gpar(fontsize = annotation_legend_label_size),
               # labels_rot = 45,
               title_gp = gpar(fontsize = annotation_legend_title_size, fontface = "bold"),
               legend_width = unit(annotation_legend_length,"cm"),
               legend_height = unit(annotation_legend_height,"cm"),
               at = column_continuous_variable_breaks.l[[variable]])
      } else{
        annotation_legend_param.l[[variable]] <-
          list(direction = "horizontal", border =T,
               labels_gp = gpar(fontsize = annotation_legend_label_size),
               # labels_rot = 45,
               title_gp = gpar(fontsize = annotation_legend_title_size, fontface = "bold"),
               legend_width = unit(annotation_legend_length,"cm"),
               legend_height = unit(annotation_legend_height,"cm")
          )
      }
    }

    column_annotation <-
      HeatmapAnnotation(
        df = annotation_data.df,
        col = column_discrete_variables_palettes.l,
        # Uncomment to provide palettes for both discrete and continuous variables
        # col = c(column_discrete_variables_palettes.l,
        #         column_continuous_variables_palettes.l
        # ),
        gap = unit(.1,"cm"),
        gp = gpar(col = "grey30", lwd = .1),
        border = T,
        simple_anno_size = unit(.3, "cm"),
        annotation_name_gp = gpar(fontsize = annotation_name_size),
        annotation_legend_param = annotation_legend_param.l,
        na_col = "grey",
        which = "column"
        # show_legend = F
      )
  }
  # ---------------------------------------------------------------------------------------------------------
  # Create annotation for row variables
  all_row_variables.v <- c(row_discrete_variables_to_plot.v, row_continuous_variables_to_plot.v)
  if (length(all_row_variables.v) > 0){
    all_row_variables.v <- all_row_variables.v[all_row_variables.v %in% colnames(row_annotation_data.df)]

    # Filter row annotation data to those entries present in the heatmap matrix
    # Assumes rownames defined
    row_annotation_data.df <- row_annotation_data.df[rownames(row_annotation_data.df) %in% rownames(input_matrix.m),,drop =F]

    # Order/filter the heatmap matrix to order/entries of metadata
    input_matrix.m <- input_matrix.m[rownames(row_annotation_data.df),,drop = F]

    # Order the heatmap matrix by the variables
    input_matrix.m <- input_matrix.m[do.call(order, row_annotation_data.df[,all_row_variables.v,drop=F]),,drop =F]

    # Order the row_annotation_data.df by the variables
    row_annotation_data.df <- row_annotation_data.df[do.call(order, row_annotation_data.df[,all_row_variables.v,drop=F]),,drop=F]

    # Check that rownames match colnames
    if (!all(rownames(row_annotation_data.df) == rownames(input_matrix.m))){
      stop("Row names in row annotation data do not match row names in matrix")
    }

    # Create final annotation dataframe from the row annotation data
    annotation_data.df <- row_annotation_data.df[,all_row_variables.v,drop = F]

    if (length(row_continuous_variables_to_plot.v) > 0){
      # Remove variables/rows that are all zero/NA
      annotation_data.df <-
        annotation_data.df %>%
        dplyr::select(-names(which(apply(annotation_data.df[,row_continuous_variables_to_plot.v,drop =F], 2, sum, na.rm = T) == 0)))
      row_continuous_variables_to_plot.v <- row_continuous_variables_to_plot.v[row_continuous_variables_to_plot.v %in% names(annotation_data.df)]
    }

    # Build the annotation_legend_param list for each variable
    annotation_legend_param.l <- list()
    for (variable in names(annotation_data.df)){
      if (variable %in% names(row_continuous_variable_breaks.l)){

        annotation_legend_param.l[[variable]] <-
          list(direction = "horizontal",border =T,
               labels_gp = gpar(fontsize = annotation_legend_label_size),
               # labels_rot = 45,
               title_gp = gpar(fontsize = annotation_legend_title_size, fontface = "bold"),
               legend_width = unit(annotation_legend_length,"cm"),
               legend_height = unit(annotation_legend_height,"cm"),
               at = row_continuous_variable_breaks.l[[variable]])
      } else{
        annotation_legend_param.l[[variable]] <-
          list(direction = "horizontal", border =T,
               labels_gp = gpar(fontsize = annotation_legend_label_size),
               # labels_rot = 45,
               title_gp = gpar(fontsize = annotation_legend_title_size, fontface = "bold"),
               legend_width = unit(annotation_legend_length,"cm"),
               legend_height = unit(annotation_legend_height,"cm")
          )
      }
    }

    row_annotation <-
      HeatmapAnnotation(
        df = annotation_data.df,
        col = row_discrete_variables_palettes.l,
        # Uncomment to provide palettes for both discrete and continuous variables
        # col = c(row_discrete_variables_palettes.l,
        #         row_continuous_variables_palettes.l
        # ),
        gap = unit(.1,"cm"),
        gp = gpar(col = "grey30", lwd = .1),
        border = T,
        simple_anno_size = unit(.3, "cm"),
        annotation_name_gp = gpar(fontsize = annotation_name_size),
        annotation_legend_param = annotation_legend_param.l,
        na_col = "grey",
        which = "row"
        # show_legend = F
      )
  }

  # if (!is.null(taxa_data.df)){
  #
  #   row_annotation <- HeatmapAnnotation(
  #     # Phylum = anno_block(gp = gpar(fill = 1:1000),
  #     # labels = c("Oct-19", "Jun-20", "Nov-20", "Jul-21"),
  #     # labels_gp = gpar(col = "black", fontsize = 10)),
  #     df = taxa_data.df %>% dplyr::select(-1,-2) %>% as.data.frame(),
  #     gap = unit(.1,"cm"),
  #     gp = gpar(col = "grey30", lwd = .1),
  #     show_annotation_name = F,
  #     annotation_legend_param = list(border = "black"),
  #     border = T,
  #     name = "Phylum",
  #     which = "row")
  # }


  # ----------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------
  # Create legends
  if (is.null(main_legend_title)){
    main_legend_title <- expression(bold(log[2]~"(Relative abundance%)"))
  }

  abundance_legend <- Legend(
    # TODO Dynamic breaks/labels
    labels = rev(c(c(0, 0.1, 0.5,5, seq(10,30,10)), ">= 40")),
    at = abundance_breaks.v,
    labels_gp = gpar(fontsize = 10),
    legend_gp = gpar(fill = rev(abundance_col_fun(abundance_breaks.v))), # For discrete
    title_position = "leftcenter-rot",
    # title_position = "topcenter",
    # title_position = "topcenter",
    title_gp = gpar(fontsize = 10,fontface = "bold"),
    title = main_legend_title,
    # title = expression(bold(log[2]~"(Estimated gene abundance in community)")),
    direction = "vertical",
    border = "black",
    # nrow = 1,
    #labels_rot = 45,
    # col_fun = col_fun
  )

  # ----------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------
  # Create main heatmap

  # Set global heatmap parameters
  ht_opt$HEATMAP_LEGEND_PADDING <- unit(heatmap_legend_padding,"cm")
  ht_opt$COLUMN_ANNO_PADDING <- unit(column_annotation_padding,"cm")
  ht_opt$ROW_ANNO_PADDING <- unit(row_annotation_padding,"cm")
  ht_opt$legend_gap <- unit(c(1,1),"cm")

  if (!is.null(column_split_variable) & !is.null(column_annotation_data.df)){
    column_split_values.v <- column_annotation_data.df[colnames(input_matrix.m), column_split_variable]
  } else{
    column_split_values.v <- NULL
  }

  if (!is.null(row_split_variable) & !is.null(row_annotation_data.df)){
    row_split_values.v <- row_annotation_data.df[,row_split_variable]
  } else{
    row_split_values.v <- NULL
  }

  if (!is.null(column_label_variable)){
    column_label_values.v <- column_annotation_data.df[colnames(input_matrix.m), column_label_variable]
  } else{
    column_label_values.v <- colnames(input_matrix.m)
  }
  if (!is.null(column_labels.df)){
    column_label_values.v <- as.character(lapply(col_labels.v, function(x) as.character(column_labels.df[column_labels.df[,1] == x,][,2])))
  }

  if (!is.null(row_label_variable)){
    row_label_values.v <- row_annotation_data.df[,row_label_variable]
  } else{
    row_label_values.v <- rownames(input_matrix.m)
  }

  if (!is.null(row_labels.df)){
    row_label_values.v <- as.character(lapply(row_label_values.v, function(x) as.character(row_labels.df[row_labels.df[,1] == x,][,2])))
  }

  # TODO FIXME
  # if (order_heatmap_columns_by_labels == T){
  #   # Order the heatmap rows by the column labels names
  #   input_matrix.m <- input_matrix.m[,order(column_label_values.v),drop =F]
  #   column_label_values.v <- column_label_values.v[order(column_label_values.v)]
  #   column_split_values.v <- column_split_values.v[order(column_label_values.v)]
  # }
  #
  # if (order_heatmap_rows_by_labels == T){
  #   # Order the heatmap rows by the row labels names
  #   input_matrix.m <- input_matrix.m[order(row_label_values.v),,drop =F]
  #   row_label_values.v <- row_label_values.v[order(row_label_values.v)]
  #   row_split_values.v <- row_split_values.v[order(row_label_values.v)]
  # }

  # if show values and no function provided
  if (show_cell_values == T & is.null(my_cell_fun)){
    my_cell_fun <- function(j, i, x, y, width, height, fill) {
      # if(input_matrix.m[i, j] < cell_fun_value_col_threshold & input_matrix.m[i, j] != 0){
      if(input_matrix.m[i, j] < cell_fun_value_col_threshold){
        grid.text(sprintf("%.2f", input_matrix.m[i, j]), x, y, gp = gpar(fontsize = 6, col = "black"))}
      else if(input_matrix.m[i, j] >= cell_fun_value_col_threshold ) {
        grid.text(sprintf("%.2f", input_matrix.m[i, j]), x, y, gp = gpar(fontsize = 6, col = "white"))
      }
    }
  }


  hm = Heatmap(input_matrix.m,
               name = "abundance",

               show_heatmap_legend = F,
               top_annotation = column_annotation,
               left_annotation = row_annotation,
               cluster_rows = cluster_rows,
               cluster_columns = cluster_columns,

               # Colours
               col = abundance_col_fun,
               na_col = "grey",

               border_gp = gpar(lwd = 1, col = "black"),
               rect_gp = gpar(lwd = .2, col = "grey10"),

               # ------------------------------------------
               # Columns
               column_title = column_title,
               column_title_side = "bottom",
               column_title_gp = gpar(fontsize = column_title_font_size, fontface = "bold"),
               # Split columns. IF FACTOR, levels are used to determine order
               column_split = column_split_values.v,
               column_labels = column_label_values.v,
               column_names_side = "bottom",
               column_names_gp = gpar(fontsize = column_labels_font_size),
               column_names_rot = 45,
               column_gap = unit(.2,"cm"),
               show_column_dend = show_column_dend,
               # ------------------------------------------
               # Rows
               row_title = row_title,
               row_title_side = "left",
               row_title_gp = gpar(fontsize = row_title_font_size, fontface = "bold"),
               # Split columns. IF FACTOR, levels are used to determine order
               row_split = row_split_values.v,
               row_labels = row_label_values.v,
               row_names_side = "right",
               row_names_gp = gpar(fontsize = row_labels_font_size),
               row_names_rot = 0,
               row_gap = unit(.2,"cm"),
               show_row_dend = show_row_dend,
               cell_fun = my_cell_fun
  )

  hm_list = hm
  if (!is.null(output_filename)){
    pdf(output_filename,height=height,width=width)
    padding <- unit(plot_padding, "cm")
    draw(hm_list,
         # main_heatmap = "rpkm",
         # legend_grouping = "original",
         # annotation_legend_list = c(abundance_legend),
         # show_annotation_legend = F,
         heatmap_legend_list = c(abundance_legend),
         merge_legend = F,
         ht_gap = unit(.5, "cm"),
         padding = padding,
         annotation_legend_side = "bottom",
         heatmap_legend_side = "right"
    )
    dev.off()
  }
  return(list("heatmap" = hm, "legend" = abundance_legend))
}


# Function to create heatmap
make_heatmap <- function(heatmap.m,
                         metadata.df = NULL, # rownames must match the columns of the heatmap.m
                         filename= NULL, # Location to save heatmap figure

                         # Dataframe with two columns. First must match row entry, second the new label
                         row_labels.df = NULL,
                         col_labels.df = NULL, # Same as row_labels, though with columns
                         height = 10, # Not currently used
                         width = 10, # Not currently used
                         heatmap_height = 10, # Not currently used
                         heatmap_width = 10, # Not currently used
                         plot_height =10,
                         plot_width =10,

                         # e.g. unit(c(2, 2, 2, 20), "mm"))
                         # bottom, left, top, right paddings
                         my_padding = NULL, # Padding around plot

                         column_title_size = 10,
                         row_title_size = 10,
                         annotation_bar_name_size = 10,
                         annotation_variables = NULL, # Annotations
                         continuous_variables = NULL, # List of discrete variables
                         discrete_annotation_palette = NULL,
                         continuous_annotation_palette = NULL,
                         simple_anno_size = unit(.5, "cm"), # size of annotations
                         col_annotation_label_size = 6,
                         col_annotation_title_size = 6,
                         col_annotation_legend_grid_height = .2,
                         col_annotation_legend_grid_width = .2,
                         plot_title = NULL,

                         cluster_columns = T,
                         cluster_rows = T,
                         my_breaks = NULL,
                         heatmap_palette = NULL,
                         default_palette_choice = NULL, # red, blue, purple, bluered, dark_bluered
                         legend_title = NULL,
                         legend_labels = NULL,
                         scale_legend_label_size = 6,
                         scale_legend_title_size = 6,
                         discrete_legend = FALSE, # Whether or not to display continuous legend as discrete

                         column_title = "Sample",
                         row_title = "Taxa",

                         show_column_dend = F,
                         show_row_dend = F,
                         do_not_order = F, # Do not order the heatmap rows by the row labels names
                         show_cell_values = F,
                         # If show_cell_values =T, cells less than this will have a black font colour
                         # and above white
                         cell_fun_value_col_threshold = 15,
                         my_cell_fun = NULL, # function to apply to the cell values
                         show_legend = T,
                         show_top_annotation = T,
                         row_name_size = 6,
                         col_name_size = 6,

                         grid_thickness = 1,
                         grid_colour = "white",

                         draw_border = F, # Can be a colour

                         draw_plot = T, # Draw the plot
                         ...){
  # print(list(...))
  argList<-list(...) # argument list for checking unspecified optional parameters
  # print(argList$cell_fun)
  # return(1)

  colour_palette_206 <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff",
                          "#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d",
                          "#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c",
                          "#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b",
                          "#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff",
                          "#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4",
                          "#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00",
                          "#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8",
                          "#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100",
                          "#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8",
                          "#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f",
                          "#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff",
                          "#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3",
                          "#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff",
                          "#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614",
                          "#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec",
                          "#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd",
                          "#93f6e7","#01b4a4")

  # Assign internal objects
  internal_heatmap_matrix.m <- heatmap.m
  internal_metadata.df <- metadata.df

  ha <- NULL

  if (!is.null(internal_metadata.df)){
    # Order/filter the heatmap matrix to order/entries of metadata
    internal_heatmap_matrix.m <- internal_heatmap_matrix.m[,rownames(internal_metadata.df),drop = F]
    if (!is.null(annotation_variables)){
      # Order the heatmap matrix by the annotation_variables
      internal_heatmap_matrix.m <- internal_heatmap_matrix.m[,do.call(order, internal_metadata.df[,annotation_variables,drop=F]),drop =F]
      # Order the metadata by the annotation_variables
      internal_metadata.df <- internal_metadata.df[do.call(order, internal_metadata.df[,annotation_variables,drop=F]),,drop=F]
      # Create metadata just containing the annotation_variables
      metadata_just_variables <- internal_metadata.df[,annotation_variables, drop = F]
    }
    # Check that rownames match colnames
    if (!all(rownames(internal_metadata.df) == colnames(internal_heatmap_matrix.m))){
      stop("Row names in metadata do not match column names in matrix")
    }

    # Create annotations
    if (is.null(discrete_annotation_palette)){ # TODO make discrete annotation and continuous annotation palette
      discrete_annotation_palette <- colour_palette_206
    }
    if (is.null(continuous_annotation_palette)){
      continuous_annotation_palette <- colorRampPalette(c("#17468a","#ffdd47","#99113a"))
    } else{
      continuous_annotation_palette <- colorRampPalette(continuous_annotation_palette)
    }

    colour_lists <- list()
    ignore_continuous_variable.v <- c()
    for (myvar in annotation_variables){
      # internal_colour_palette_10_distinct <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
      # internal_colour_palette_10_distinct <- my_colour_palette_20_distinct
      if (!is.null(continuous_variables) & myvar %in% continuous_variables){
        variable_values.v <- unique(internal_metadata.df[,myvar])
        internal_breaks <- seq(min(variable_values.v, na.rm = T), max(variable_values.v,na.rm = T), length.out = 10)
        # col_fun <- circlize::colorRamp2(breaks = internal_breaks,
        # colors = continuous_annotation_palette(length(internal_breaks)))
        if (length(unique(internal_breaks)) == 1) {
          warning(paste0("Cannot assign colours to variable: ", myvar, " as there is only one unique value"))
          ignore_continuous_variable.v <- c(myvar, ignore_continuous_variable.v)
          next
        } else{
          col_fun <- colorRamp2(breaks = internal_breaks,
                                colors = continuous_annotation_palette(length(internal_breaks)))
          colour_lists[[myvar]] <- col_fun
        }
      } else{
        # Assumes there is a colour column for each variable in the metadata
        # If there is no colour column, create one and assign from palette
        var_colour_name <- paste0(myvar, "_colour")
        if (!var_colour_name %in% names(internal_metadata.df)){
          myvar_values <- factor(as.character(sort(unique(internal_metadata.df[,myvar]))))
          myvar_colours <- setNames(discrete_annotation_palette[1:length(myvar_values)], myvar_values)
          all_variable_colours <- as.character(lapply(as.character(internal_metadata.df[,myvar]), function(x) myvar_colours[x]))
          internal_metadata.df[,paste0(myvar,"_colour")] <- all_variable_colours
        }

        metadata_subset <- unique(internal_metadata.df[,c(myvar, var_colour_name)])
        # Order by the variable column
        metadata_subset <- metadata_subset[order(metadata_subset[,myvar]),]
        # Factorise the variable column
        metadata_subset[,myvar] <- factor(metadata_subset[,myvar])
        metadata_subset <- metadata_subset[!is.na(metadata_subset[,myvar]),]
        named_colour_list <- setNames(as.character(metadata_subset[, var_colour_name]), as.character(metadata_subset[,myvar]))
        colour_lists[[myvar]] <- named_colour_list
      }
    }
    for (myvar in ignore_continuous_variable.v){
      colour_lists[[myvar]] <- NULL
      metadata_just_variables <- metadata_just_variables %>% select(-all_of(myvar))
    }

    # Appearance of the column annotations
    #HeatmapAnnotation
    if (show_top_annotation == T){
      if(is.null(annotation_variables)){
        print("No annotation variables specified, cannot add annotations")
        show_top_annotation <- F
      } else{
        # TODO make flexible for different types of annotations, e.g. bar
        # Or allow user to provide annotations that are added in, 'provided_ha'
        ha <- columnAnnotation(df = metadata_just_variables,
                               # which = "column",
                               na_col = "grey",
                               col = colour_lists,
                               gp = gpar(col = "black",lwd =.2),
                               gap = unit(.1,"cm"),
                               show_annotation_name = T,
                               # annotation_legend_param, # ?color_mapping_legend for options
                               show_legend = show_legend,
                               simple_anno_size = simple_anno_size,
                               annotation_legend_param = list(labels_gp = gpar(fontsize = col_annotation_label_size),
                                                              title_gp = gpar(fontsize = col_annotation_title_size,fontface = "bold"),
                                                              grid_height = unit(col_annotation_legend_grid_height,"cm"),
                                                              grid_width = unit(col_annotation_legend_grid_width,"cm")),
                               annotation_name_gp = gpar(fontsize = annotation_bar_name_size))
      }
    }
  }

  # TODO - add option for row annotation
  if (is.null(heatmap_palette)){
    if (is.null(default_palette_choice)) {default_palette_choice <- "blue"}
    if (!default_palette_choice %in% c("blue", "purple","red","dark_bluered","bluered")) { default_palette_choice <- "blue"}
    if (default_palette_choice == "blue"){
      heatmap_palette <- colorRampPalette(c("white", "#ffffcc","#cce1b8", "#91cabc", "#61b4c1","#335fa5","#28387a", "#071447"))
    }
    else if (default_palette_choice == "purple"){
      heatmap_palette <- colorRampPalette(c("white", "#f9cdac","#f3aca2", "#ee8b97", "#e96a8d","#db5087","#b8428c", "#973490", "#742796","#5e1f88", "#4d1a70", "#3d1459","#2d0f41"))
    } else if (default_palette_choice == "red"){
      heatmap_palette <- colorRampPalette(c("white", "#fded86","#fde86e", "#f9d063", "#f5b857","#f0a04b","#eb8a40", "#e77235","#e35b2c", "#c74e29","#9d4429","#753c2c","#4c3430"))
    } else if (default_palette_choice == "dark_bluered"){
      heatmap_palette <- colorRampPalette(c("#08306B","#FFD92F","#67001F"))
    } else if (default_palette_choice == "bluered"){
      heatmap_palette <- colorRampPalette(c("#17468a","#ffdd47","#99113a"))
    }
  } else{
    heatmap_palette <- colorRampPalette(heatmap_palette)
  }

  if (!is.null(my_breaks)){
    internal_breaks <- my_breaks
    col_fun <- circlize::colorRamp2(breaks = internal_breaks, colors = heatmap_palette(length(internal_breaks)))
  } else{
    internal_breaks <- seq(min(internal_heatmap_matrix.m), max(internal_heatmap_matrix.m), length.out = 6)
    col_fun <- circlize::colorRamp2(breaks = internal_breaks, colors = heatmap_palette(length(internal_breaks)))
  }

  row_labels.v = rownames(internal_heatmap_matrix.m)
  if (!is.null(row_labels.df)){
    row_labels.v <- as.character(lapply(row_labels.v, function(x) as.character(row_labels.df[row_labels.df[,1] == x,][,2])))
  }
  col_labels.v = colnames(internal_heatmap_matrix.m)
  if (!is.null(col_labels.df)){
    col_labels.v <- as.character(lapply(col_labels.v, function(x) as.character(col_labels.df[col_labels.df[,1] == x,][,2])))
  }

  if (do_not_order != T){
    # Order the heatmap rows by the row labels names
    internal_heatmap_matrix.m <- internal_heatmap_matrix.m[order(row_labels.v),,drop =F]
    row_labels.v <- row_labels.v[order(row_labels.v)]
  }

  # if show values and no function provided
  if (show_cell_values == T & is.null(my_cell_fun)){
    my_cell_fun <- function(j, i, x, y, width, height, fill) {
      # if(internal_heatmap_matrix.m[i, j] < cell_fun_value_col_threshold & internal_heatmap_matrix.m[i, j] != 0){
      if(internal_heatmap_matrix.m[i, j] < cell_fun_value_col_threshold){
        grid.text(sprintf("%.2f", internal_heatmap_matrix.m[i, j]), x, y, gp = gpar(fontsize = 6, col = "black"))}
      else if(internal_heatmap_matrix.m[i, j] >= cell_fun_value_col_threshold ) {
        grid.text(sprintf("%.2f", internal_heatmap_matrix.m[i, j]), x, y, gp = gpar(fontsize = 6, col = "white"))
      }
    }
  }

  # Legend appearance
  if (is.null(legend_labels)){
    my_labels <- internal_breaks
  } else{
    my_labels <- legend_labels
  }
  if (discrete_legend == TRUE){
    hm_legend <- Legend(
      labels = rev(my_labels),
      at = internal_breaks,
      labels_gp = gpar(fontsize = scale_legend_label_size),
      legend_gp = gpar(fill = rev(col_fun(internal_breaks))), # For discrete
      title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = scale_legend_title_size,fontface = "bold"),
      title = legend_title,
      direction = "vertical",
      border = "black"
    )
  } else{
    hm_legend <- Legend(
      col_fun = col_fun, # For continuous
      labels = my_labels,
      at = internal_breaks,
      labels_gp = gpar(fontsize = scale_legend_label_size),
      title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = scale_legend_title_size,fontface = "bold"),
      title = legend_title,
      direction = "vertical",
      border = "black"
    )
  }

  hm <- Heatmap(matrix = internal_heatmap_matrix.m,

                top_annotation = ha,

                # Colours
                col = col_fun,
                na_col = "grey",

                # Sizing
                show_heatmap_legend = F,
                heatmap_legend_param = list(hm_legend),
                row_names_max_width = unit(35,"cm"),
                row_labels = row_labels.v,
                column_labels = col_labels.v,
                # row_names_side = "left",
                # height = unit(height,"cm"),
                # width = unit(width,"cm"),
                # heatmap_height = unit(heatmap_height,"cm"),
                # heatmap_width = unit(heatmap_width,"cm"),
                # heatmap_width = unit(15,"cm"),

                # Titles
                column_title = column_title,
                column_title_side = "bottom",
                column_title_gp = gpar(fontsize = column_title_size),
                row_title = row_title,
                row_title_side = "left",
                row_title_gp = gpar(fontsize = row_title_size),

                # Clustering
                cluster_columns = cluster_columns,
                cluster_rows = cluster_rows,
                clustering_method_columns = "average",
                clustering_method_rows = "average",
                show_column_dend = show_column_dend,
                show_row_dend = show_row_dend,
                # column_dend_height = unit(2, "cm"),
                # row_dend_width = unit(3, "cm"),

                # Borders
                border = draw_border,
                rect_gp = gpar(col = grid_colour, lwd = grid_thickness),

                # Text appearance
                row_names_gp = gpar(fontsize = row_name_size),
                column_names_gp = gpar(fontsize = col_name_size),
                cell_fun = my_cell_fun,
                ...
  )

  # if (show_top_annotation == T & exists("ha")){
  #   hm <- ha %v% hm
  # }
  padding <- unit(c(2, 2, 2, 2), "mm")
  if (!is.null(my_padding)){
    padding <- my_padding
  }

  if (!is.null(plot_title)){
    column_title <- plot_title
  } else{
    column_title <- character(0)
  }
  if (!is.null(filename)){
    pdf(filename,height=plot_height,width=plot_width)
    draw(hm, annotation_legend_list = c(hm_legend),merge_legends =T,padding = padding, column_title = column_title,
         column_title_gp = gpar(hjust = 0))
    dev.off()
  }
  if (draw_plot == T){
    draw(hm, annotation_legend_list = c(hm_legend),merge_legends =T,padding = padding,column_title = column_title)
  }
  return(list("heatmap" = hm, "legend" = hm_legend))

}
