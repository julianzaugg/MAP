

#' Darken a colour
#' @param colour colour name, hexadecimal string or positive integer \code{i} meaning \code{palette()[i]}
#' @param factor value to scale rgb colour vector by division
darken <- function(colour, factor=1.4){
  col <- col2rgb(colour)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

#' Lighten a colour
#' @param colour colour name, hexadecimal string or positive integer \code{i} meaning \code{palette()[i]}
#' @param factor value to scale rgb colour vector by multiplication
lighten <- function(colour, factor=1.4){
  col <- col2rgb(colour)
  col <- col*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


#' Create colour palettes for continuous variables
#' @param data.df dataframe containing continuous variables to
#' @param variables.v variables/columns to make continuous
#' @return Returns a list of functions, for each variable. Functions accept a vector of numeric values
#' and returns interpolated colours.
make_continuous_palette <- function(data.df, variables.v){

  colorRamp2 = function(breaks, colors, transparency = 0, space = "LAB") {
    # Taken from https://github.com/jokergoo/circlize/

    if(length(breaks) != length(colors)) {
      rlang::abort("Length of `breaks` should be equal to `colors`.\n")
    }

    colors = colors[order(breaks)]
    breaks = sort(breaks)

    l = duplicated(breaks)
    breaks = breaks[!l]
    colors = colors[!l]

    if(length(breaks) == 1) {
      rlang::abort("You should have at least two distinct break values.")
    }


    if(! space %in% c("RGB", "HSV", "HLS", "LAB", "XYZ", "sRGB", "LUV")) {
      rlang::abort("`space` should be in 'RGB', 'HSV', 'HLS', 'LAB', 'XYZ', 'sRGB', 'LUV'")
    }

    colors = t(col2rgb(colors)/255)

    attr = list(breaks = breaks, colors = colors, transparency = transparency, space = space)

    if(space == "LUV") {
      i = which(apply(colors, 1, function(x) all(x == 0)))
      colors[i, ] = 1e-5
    }

    transparency = 1-ifelse(transparency > 1, 1, ifelse(transparency < 0, 0, transparency))[1]
    transparency_str = sprintf("%X", round(transparency*255))
    if(nchar(transparency_str) == 1) transparency_str = paste0("0", transparency_str)

    fun = function(x = NULL, return_rgb = FALSE, max_value = 1) {
      if(is.null(x)) {
        rlang::abort("Please specify `x`\n")
      }

      att = attributes(x)
      if(is.data.frame(x)) x = as.matrix(x)

      l_na = is.na(x)
      if(all(l_na)) {
        return(rep(NA, length(l_na)))
      }

      x2 = x[!l_na]

      x2 = ifelse(x2 < breaks[1], breaks[1],
                  ifelse(x2 > breaks[length(breaks)], breaks[length(breaks)],
                         x2
                  ))
      ibin = .bincode(x2, breaks, right = TRUE, include.lowest = TRUE)
      res_col = character(length(x2))
      for(i in unique(ibin)) {
        l = ibin == i
        res_col[l] = .get_color(x2[l], breaks[i], breaks[i+1], colors[i, ], colors[i+1, ], space = space)
      }
      res_col = paste(res_col, transparency_str[1], sep = "")

      if(return_rgb) {
        res_col = t(col2rgb(as.vector(res_col), alpha = TRUE)/255)
        return(res_col)
      } else {
        res_col2 = character(length(x))
        res_col2[l_na] = NA
        res_col2[!l_na] = res_col

        attributes(res_col2) = att
        return(res_col2)
      }
    }

    attributes(fun) = attr
    return(fun)
  }

  # set.seed(10)
  continuous_annotation_palette <- colorRampPalette(c("#17468a","#ffdd47","#99113a"))
  palettes.l <- list()
  for (variable in variables.v){
    variable_breaks.v <- seq(min(data.df[,variable], na.rm = T), max(data.df[,variable],na.rm = T), length.out = 10)
    col_fun <- colorRamp2(breaks = variable_breaks.v,
                          colors = continuous_annotation_palette(length(variable_breaks.v)))
    palettes.l[[variable]] <- col_fun
  }
  palettes.l
}
