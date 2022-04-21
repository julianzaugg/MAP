
#' Generate correlation matrix. Probably better/faster to instead use Hmisc::rcorr...
calculate_correlation_matrix <- function(mydata.m, method = "pearson", adjust = "BH"){

  # Remove row entries that don't vary across all samples
  internal_data.m <- mydata.m
  zv <- apply(internal_data.m, 1, function(x) length(unique(x)) == 1)
  internal_data.m <- internal_data.m[!zv, ]

  # Take a two vectors and perform a significance/correlation test
  calculate_stats <- function(x, y, dist.name) {
    k <- cor.test(x, y, method=dist.name)
    c(k$estimate, k$stat, k$p.value)
  }

  # cor_result <- corr.test(t(internal_data.m), adjust = adjust, method = method)
  r <- apply(t(internal_data.m), 2,function(col)
    t(apply(t(internal_data.m), 2, calculate_stats, col, dist.name = method))[,1,drop =F])
  p <- apply(t(internal_data.m), 2,function(col)
    t(apply(t(internal_data.m), 2, calculate_stats, col, dist.name = method))[,3,drop =F])
  padj <- apply(t(internal_data.m), 2,function(col)
    p.adjust(t(apply(t(internal_data.m), 2, calculate_stats, col, dist.name = method))[,3,drop =F],
             method = adjust))

  rownames(r) <- colnames(r)
  rownames(p) <- colnames(p)
  rownames(padj) <- colnames(padj)
  r <- signif(r, digits = 5)
  p <- signif(p, digits = 5)
  padj <- signif(padj, digits = 5)

  # r	The matrix of correlations
  # n	Number of cases per correlation
  # t	value of t-test for each correlation
  # p	two tailed probability of t for each correlation. For symmetric matrices, p values adjusted for multiple tests are reported above the diagonal.
  # se	standard error of the correlation
  # ci	the alpha/2 lower and upper values, as well as the (Holm or Bonferroni) adjusted confidence intervals.
  # list(cor_matrix = cor_result$r, cor_pval_matrix = cor_result$p)
  list(cor_matrix = r, cor_pval_matrix = p, cor_padj_matrix = padj)
}


generate_correlation_network_from_count_data <- function(counts.m, p_value_threshold = 0.05, cor_threshold = 0.5, method = "pearson", adjust = "BH"){
  correlation_results <- calculate_correlation_matrix(counts.m, adjust = adjust, method = method)
  cor.m <- correlation_results$cor_matrix
  cor_pval.m <- correlation_results$cor_pval_matrix
  generate_correlation_network(cor.m, cor_pval.m,p_value_threshold = 0.05, cor_threshold = 0.5)
}

#' Flatten a correlation and corresponding p-value
#' @param cormat matrix of the correlation coefficients
#' @param pmat matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  # cormat : matrix of the correlation coefficients
  # pmat : matrix of the correlation p-values
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}




#' Generate a correlation network
generate_correlation_network <- function(cor_matrix.m, p_matrix.m = NULL, p_value_threshold = 0.05, cor_threshold = 0.5,
                                         filename = NULL, relabeller_function = NULL,
                                         node_size = 4, node_colour = "grey20", node_fill = "grey20", node_shape = 21,
                                         node_label_segment_colour = "grey",
                                         label_colour = "black",label_size = 1,
                                         plot_height = 10, plot_width = 10, plot_title = "",
                                         edge_width_min = .3, edge_width_max = 1,edge_alpha =1,
                                         network_layout = "fr", exclude_to_from_df = NULL,
                                         myseed=NULL, edgetype = "link", show_p_label=F,
                                         file_type = "pdf"){
  if (!is.null(myseed)){
    set.seed(myseed)
  }

  cor.m <- as.matrix(cor_matrix.m)
  if (max(cor.m[row(cor.m)!=col(cor.m)],na.rm = T) < cor_threshold){
    print("No entries left after filtering at correlation threshold")
    return()
  }
  if (!is.null(relabeller_function)){
    colnames(cor.m) <- lapply(colnames(cor.m), relabeller_function)
    rownames(cor.m) <- lapply(rownames(cor.m), relabeller_function)
  }

  if (!is.null(p_matrix.m)){
    if (min(p_matrix.m[row(p_matrix.m)!=col(p_matrix.m)],na.rm = T) > p_value_threshold){
      print("No entries left after filtering at p-value threshold")
      return()
    }
    cor_pval.m <- as.matrix(p_matrix.m)
    if (!is.null(relabeller_function)){
      # colnames(cor_pval.m) <- relabeller_function(colnames(cor_pval.m))
      # rownames(cor_pval.m) <- relabeller_function(rownames(cor_pval.m))
      colnames(cor_pval.m) <- lapply(colnames(cor_pval.m), relabeller_function)
      rownames(cor_pval.m) <- lapply(rownames(cor_pval.m), relabeller_function)
    }
  }

  # Melt correlation matrix
  graph.df <- melt(cor.m,value.name = "Correlation",varnames = c("Variable_1", "Variable_2"))
  if (!is.null(p_matrix.m)){
    graph.df$P_value <- melt(cor_pval.m)$value
    graph.df <- graph.df[graph.df$P_value <= p_value_threshold,]
    graph.df$P_value_label <- NA
    if (show_p_label==T){
      graph.df$P_value_label <- lapply(graph.df$P_value, function(x) ifelse(x <= 0.001, "***", ifelse(x <= 0.01, "**", ifelse(x <= 0.05, "*", ""))))
    }
    # graph.df$P_value_label <- round(graph.df$Correlation,2) # For sanity testing

  } else{
    graph.df$P_value_label <- NA
    # graph.df$P_value_label <- round(graph.df$Correlation,2) # For sanity testing
    names(graph.df) <- c("Variable_1", "Variable_2", "Correlation", "P_value_label")
  }

  # Remove bidirectional links
  graph.df$ordered_label <- apply(graph.df, 1, function(x) {paste0(sort(c(as.character(x[1]), as.character(x[2]))),collapse = ":")})
  graph.df <- graph.df[!duplicated(graph.df$ordered_label),]

  # Remove edges with a correlation of 0, regardless what threshold
  graph.df <- graph.df[graph.df$Correlation != 0,]

  # Remove edges that are below the correlation threshold
  graph.df <- graph.df[abs(graph.df$Correlation) >= cor_threshold,]

  # Remove edges specified in supplied dataframe.
  if (!is.null(exclude_to_from_df)){
    graph.df <- graph.df[!paste0(graph.df$Variable_1, "-", graph.df$Variable_2) %in% paste0(edges_to_remove.df[,1], "-", edges_to_remove.df[,2]),]
    graph.df <- graph.df[!paste0(graph.df$Variable_1, "-", graph.df$Variable_2) %in% paste0(edges_to_remove.df[,2], "-", edges_to_remove.df[,1]),]
  }

  # Generate graph object and remove looped edges and isolated nodes
  graph.df <-
    tidygraph::as_tbl_graph(graph.df) %>%
    # Remove loops
    tidygraph::activate(edges) %>%
    filter(!edge_is_loop()) %>%
    # filter(!edge_is_mutual()) %>%
    # filter(!edge_is_multiple()) %>%
    # Remove isolated nodes
    tidygraph::activate(nodes) %>%
    filter(!node_is_isolated())

  # Build plot
  set_graph_style(plot_margin = margin(1,1,1,1))
  correlation_graph_plot <- ggraph(graph.df, layout = network_layout)
  if (edgetype == "link"){
    correlation_graph_plot <- correlation_graph_plot +
      geom_edge_link(aes(colour = Correlation, width=abs(Correlation), label = P_value_label),  show.legend = T, alpha = edge_alpha,
                     angle_calc = "along", label_colour = "black",
                     label_dodge=unit(1,"mm"),label_push=unit(-1,"mm"))
  } else if (edgetype == "fan"){
    correlation_graph_plot <- correlation_graph_plot +
      geom_edge_fan(aes(colour = Correlation, width=abs(Correlation), label = P_value_label), show.legend = T, alpha = edge_alpha,
                    angle_calc = "along", label_colour = "black",
                    label_dodge=unit(1,"mm"),label_push=unit(-1,"mm"))
  } else if (edgetype == "elbow"){
    correlation_graph_plot <- correlation_graph_plot +
      geom_edge_elbow(aes(colour = Correlation, width=abs(Correlation), label = P_value_label), show.legend = T, alpha = edge_alpha, strength = 1,
                      angle_calc = "along", label_colour = "black",
                      label_dodge=unit(1,"mm"),label_push=unit(-1,"mm"))
  } else if (edgetype == "bend"){
    correlation_graph_plot <- correlation_graph_plot +
      geom_edge_bend(aes(colour = Correlation, width=abs(Correlation), label = P_value_label), show.legend = T, alpha = edge_alpha, strength = 1,
                     angle_calc = "along", label_colour = "black",
                     label_dodge=unit(1,"mm"),label_push=unit(-1,"mm"))
  } else if (edgetype == "hive"){
    correlation_graph_plot <- correlation_graph_plot +
      geom_edge_hive(aes(colour = Correlation, width=abs(Correlation), label = P_value_label), show.legend = T, alpha = edge_alpha, strength = 1,
                     angle_calc = "along", label_colour = "black",
                     label_dodge=unit(1,"mm"),label_push=unit(-1,"mm"))
  }
  correlation_graph_plot <- correlation_graph_plot +
    scale_edge_width_continuous(name="Correlation", range = c(edge_width_min,edge_width_max),
                                breaks = seq(-1,1,.2)) +
    scale_edge_colour_gradientn(colours = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D",
                                                                 "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                                 "#4393C3", "#2166AC", "#053061")))(11),
                                limits = c(-1,1), # limit colours to full corrrelation range
                                breaks = seq(-1,1,.2), # Break colours from -1 to 1 in steps of 0.2
                                guide = guide_edge_colourbar(barwidth = 0.5, barheight = 10)) +
    geom_node_point(colour = node_colour, fill = node_fill, pch = node_shape, size =node_size)  +
    geom_node_text(aes(label = name), colour = label_colour,
                   size = label_size, nudge_y = -.01, #fontface = "bold",
                   point.padding = unit(0.3, "lines"),
                   segment.size = 0.3,
                   segment.colour = node_label_segment_colour,
                   # family = "serif",
                   repel = T) +
    ggtitle(label = plot_title) +
    theme_graph(background = "white")
  # Override edge colour to show coloured width, entries for width need to match the
  # breaks defined in the colour gradient (at least in length)
  # guides(edge_color = guide_legend(override.aes = list(edge_width = abs(seq(-1,1,.2)*1.2))),
  break_length <- length(seq(-1,1,.2))
  edge_widths <- abs(c(rev(rev(seq(-edge_width_max, -edge_width_min, length.out = break_length/2))[-1]),
                       rev(seq(edge_width_max, edge_width_min, length.out = break_length/2))))
  # edge_widths <- abs(seq(-1,1,.2))
  # edge_widths <- abs(edge_widths *seq(-1,1,.2))
  correlation_graph_plot <- correlation_graph_plot +
    guides(edge_color = guide_legend(override.aes = list(edge_width = edge_widths)),
           edge_width = F)

  # Save plot to file
  if (!is.null(filename)){
    if (file_type == "pdf"){
      cairo_pdf(filename = filename,height = plot_height, width = plot_width)
      plot(correlation_graph_plot)
      dev.off()
    } else if (file_type == "svg"){
      # Cairo::CairoSVG(file = filename,width = plot_width,height = plot_height)
      # svg(filename = filename,height = plot_height, width = plot_width)
      svglite(file = filename,height = plot_height, width = plot_width)
      plot(correlation_graph_plot)
      dev.off()
    }
  }
  list(network_data = graph.df, network_plot = correlation_graph_plot)
}

#' Generate a dot correlation plot generated by corrplot
#' @param correlation_matrix
plot_corrplot <- function(correlation_matrix,
                          p_value_matrix = NULL,
                          plot_title = "",
                          plot_title_size = .6,
                          plot_height = 10, plot_width = 10,
                          p_value_threshold = 0.05,
                          label_size = 1,
                          relabeller_function = NULL,
                          insig = "blank", insig_pch = 4, insig_pch_cex = 1, insig_pch_col = "black",
                          make_insig_na = F,
                          to_exclude = NULL,
                          method = "circle",
                          outline = T,
                          label_colour = "black",
                          colour_label_size = 1,
                          grid_colour = "black",
                          pairs_to_na = NULL, # Should be two column dataframe (row column / column row)
                          order = "hclust",
                          col = NULL,
                          filename = NULL,
                          file_type = "pdf",
                          ...
){

  cor.m <- correlation_matrix
  cor_pval.m <- p_value_matrix

  # Entries to remove completely
  if (!is.null(to_exclude)){
    cor.m <- cor.m[!rownames(cor.m) %in% to_exclude, !colnames(cor.m) %in% to_exclude]
    if (!is.null(p_value_matrix)){
      cor_pval.m <- cor_pval.m[!rownames(cor_pval.m) %in% to_exclude, !colnames(cor_pval.m) %in% to_exclude]
    }

  }
  # Entries to make NA. Should be a two column dataframe
  if (!is.null(pairs_to_na)){
    if (order == "hclust"){
      # stop("Cannot use hclust ordering with NA values in matrix")
      # So before inserting NA values, first order the matrix
      # print(summary(hclust(dist(cor.m),method = "average")))
      ord <- hclust(dist(cor.m),method = "average")$order
      cor.m <- cor.m[ord,ord]
      if (!is.null(p_value_matrix)){
        cor_pval.m <- cor_pval.m[ord,ord]
      }
      order <- "original"
    }
    for (row in 1:nrow(pairs_to_na)){
      a <- as.character(pairs_to_na[row,1])
      b <- as.character(pairs_to_na[row,2])
      if (a %in% rownames(cor.m) & b %in% rownames(cor.m)){
        cor.m[a, b] <- NA
        cor.m[b, a] <- NA
        if (!is.null(p_value_matrix)){
          cor_pval.m[a, b] <- NA
          cor_pval.m[b, a] <- NA
        }
      }
    }
  }

  if (make_insig_na == T){
    cor.m[which(cor_pval.m >= p_value_threshold)] <- NA
    cor_pval.m[cor_pval.m >= p_value_threshold] <- NA
  }

  if (!is.null(relabeller_function)){
    # print(unlist(lapply(colnames(cor.m)[1:10], relabeller_function)))
    # print(unlist(lapply(rownames(cor.m)[1:10], relabeller_function)))
    colnames(cor.m) <- unlist(lapply(colnames(cor.m), relabeller_function))
    rownames(cor.m) <- unlist(lapply(rownames(cor.m), relabeller_function))
    if (!is.null(p_value_matrix)){
      colnames(cor_pval.m) <- unlist(lapply(colnames(cor_pval.m), relabeller_function))
      rownames(cor_pval.m) <- unlist(lapply(rownames(cor_pval.m), relabeller_function))
      # colnames(cor_pval.m) <- relabeller_function(colnames(cor_pval.m))
      # rownames(cor_pval.m) <- relabeller_function(rownames(cor_pval.m))
    }
  }
  if (is.null(col)){
    col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D",
                                  "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                  "#4393C3", "#2166AC", "#053061")))(200)
  }

  if (!is.null(filename)){
    # pdf(file = filename,height = plot_height, width = plot_width)
    if (file_type == "pdf"){
      pdf(file = filename, height = plot_height, width = plot_width)
    } else if (file_type == "svg"){
      # Cairo::CairoSVG(file = filename,width = plot_width,height = plot_height)
      # svg(filename = filename,height = plot_height, width = plot_width)
      svglite(file = filename,height = plot_height, width = plot_width)
    }
  }

  corrplot(corr = cor.m,
           method = method,
           outline = outline,
           tl.col = label_colour,
           tl.cex = label_size,
           addgrid.col = grid_colour,
           # tl.srt = 45,
           # title = plot_title,
           col = col,
           # col = brewer.pal(n = 8, name = "RdYlBu"),
           type = "lower",
           diag = F,
           na.label = "square",
           na.label.col = "grey",
           order = order,
           hclust.method = "average",
           p.mat = cor_pval.m,
           sig.level = p_value_threshold,
           # insig = "blank",
           insig = insig,
           pch = insig_pch,
           pch.cex = insig_pch_cex,
           pch.col = insig_pch_col,
           cl.pos = 'r',
           cl.cex = colour_label_size,
           mar=c(1,0,3,1),
           ...)
  title(main = plot_title,cex.main = plot_title_size)
  if (!is.null(filename)){
    dev.off()
  }
}
