###############################################################################
## Plotting Utilities

get_lims <- function(x) {
  ## Get axis limits for a set of points
  ## INPUTS:
  ##  x = numeric vector
  ## RETURNS:
  range_x <- max(x) - min(x)
  buff <- 0.05*range_x
  return(c(min(x) - buff, max(x) + buff))
}

#' generic scatterplot function for data frame
#' @param df = data.frame
#' @param var1 = x axis variable
#' @param var2 = y axis variable
#' @param color_by = variable to color points by
#' @param filter_by = subset data by a particular categorical variable
#' @param filter_class = set of allowed classes for subset
#' @param pt.size
#' @param alpha
#' @param legend_pt_size
#' @param legend_pt_shape
#' @param xlim
#' @param ylim
#' @value  ggplot2 object
#' @export

scatter_plot <- function(df, var1, var2, color_by,
                        filter_by = NULL, filter_class = NULL,
                        pt.size = 0.5, pt.shape =  16, alpha = 0.4,
                        legend_pt_size = 20, legend_pt_shape = 15,
                        xlim = NULL, ylim = NULL, legend_title = FALSE) {
  ## generic scatterplot function for data frame
  ## INPUTS:
  ##  df = data.frame
  ##  var1 = x axis variable
  ##  var2 = y axis variable
  ##  color_by = variable to color points by
  ##    filter_by = subset data by a particular categorical variable
  ##    filter_class = set of allowed classes for subset
  ##    pt.size, alpha, legend_pt_size, legend_pt_shape, xlim, ylim: arguments
  ##    passed to ggplot2 plotting functions. Note that if xlim and ylim are set to null,
  ##    by default the limits are taken as (min(x) - 0.05*range(x), max(x) + 0.05*range(x))
  ##    for the x limits and the same for the y limits but replacing x with y.
  ##    the limits are calculated with the whole set of points prior to filtering
  ##    Note that argument alpha is ignored if we have a continuous scale AND the range
  ##    of values includes negative and positive values (i.e crosses 0)
  ## RETURNS:
  ##  ggplot2 object
  # set x and y limits if NULL

  if (is.null(xlim)) {
    xlim <- get_lims(df[,var1])
  }
  if (is.null(ylim)) {
    ylim <- get_lims(df[,var2])
  }

  tryCatch({df[,color_by]},
           error = function(e) {
             # since warnings not going to stderr file for Rmd, if feature is not found in
             # df we specify an error that says the name of the column not found
             stop(paste('could not find feature', color_by, 'in metadata'))
           })
  # decide if dealing with categorical data, set appropriate settings
  is_categorical <- !is.numeric(df[,color_by])

  # make plot
  if (is_categorical) {
    # categorical values
    discrete_col_mapping <- discrete_color_mapping(df[,color_by])
    # filter data as desired after getting appropriate color mapping
    df <- filter_df(df = df, filter_by = filter_by, filter_class = filter_class)
    p <- ggplot(data = df, mapping = aes_(x = as.name(var1), y = as.name(var2))) +
      geom_point(aes_(color = as.name(color_by)), size = pt.size, alpha = alpha, pch = pt.shape)
    p <- p + scale_color_manual(values = discrete_col_mapping)
    p <- p + guides(colour = guide_legend(override.aes = list(size=legend_pt_size,
                                                              shape=legend_pt_shape),
                                          title = ifelse(legend_title, color_by, '')))
  } else {
    # numeric values
    qtl_vals <- quantile(df[,color_by], seq(0,1, 0.05))
    # set values for color scale
    min_val <- qtl_vals['5%']
    max_val <- qtl_vals['95%']
    if (sign(min_val) == -sign(max_val)) {
      # continuous_colors <- c('blue', 'white', 'red')
      continuous_breaks <- c(qtl_vals['5%'], 0, qtl_vals['95%'])
      zero_center <- TRUE
    } else {
      continuous_colors <- c('red', 'yellow')
      continuous_breaks <- c(qtl_vals['5%'], qtl_vals['95%'])
      zero_center <- FALSE
    }
    # set color values
    df$color_vals <- df[,color_by]
    df[df[,color_by] >= max_val, 'color_vals'] <- max_val
    df[df[,color_by] <= min_val, 'color_vals'] <- min_val
    # if (zero_center) {
    #   # set alpha values to absolute value of feature of interest, using top 5%/bottom 25% as
    #   # max/min values for alpha, then scale to be between 0 and 1
    #   alpha_vals <- abs(df[,color_by])
    #   alpha_quantiles <- quantile(alpha_vals, seq(0,1, 0.05))
    #   pct5_alpha <- alpha_quantiles['5%']
    #   pct95_alpha <- alpha_quantiles['95%']
    #   alpha_vals[alpha_vals > pct95_alpha] <- pct95_alpha
    #   alpha_vals[alpha_vals < pct5_alpha] <- pct5_alpha
    #   alpha_vals <- min(alpha_vals/pct95_alpha, 0.5*pct95_alpha)
    #   df$alpha_vals <- alpha_vals
    # }
    # filter data as desired after deciding appropriate color mapping
    df <- filter_df(df = df, filter_by = filter_by, filter_class = filter_class)
    # do the plot
    p <- ggplot(data = df, mapping = aes_(x = as.name(var1), y = as.name(var2)))
    if (zero_center) {
      p <- p + geom_point(aes_(color = df$color_vals), size = pt.size, alpha = alpha, pch = pt.shape)
      # p <- p + geom_point(aes_(color = df$color_vals), size = pt.size, alpha = alpha, pch = pt.shape)
      # p <- p + scale_color_gradient2(midpoint = 0, limits = c(min_val, max_val), low = 'turquoise', mid = 'white', high = 'orange')
      p <- p + scale_color_gradient2(midpoint = 0, low = 'blue', mid = 'darkgrey', high = 'red')
      # p <- p + theme_dark()
      p <- p + theme_light()
    } else {
      p <- p + geom_point(aes_(color = df$color_vals), size = pt.size, alpha = alpha, pch = pt.shape)
      # p <- p + scale_color_gradient2(low = 'turquoise', mid = 'white', high = 'orange',
      #                                midpoint = mean(df$color_vals),
      #                                limits = c(min_val, max_val))
      p <- p + scale_color_gradient2(low = 'darkblue', mid = 'darkgrey', high = 'red', midpoint = mean(df$color_vals))
      # p <- p + theme_dark()
      p <- p + theme_light()
    }
    if (!legend_title) {
      p <- p + guides(color = guide_colorbar(title = NULL))
    } else {
      p <- p + guides(color = guide_colorbar(title = color_by))
    }
    # p <- p + scale_color_gradientn(colors = continuous_colors, breaks = continuous_breaks)

  }
  p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  return(p)
}

#' Plot Reduced Dims
#' @param x SummExpDR object
#' @param dim1
#' @param dim2
#' @param color_by
#' @param key key for reduced dims to use
#' @param assay required if pulling ouot feature values for assay data
#' @param ... other arguments to scatter_plot, with var1 and var2 overridden by dim1 and dim2
#' @value ggplot2 object
#' @export

setGeneric('plotDR', function(x, dim1, dim2, color_by, key, assay = NULL, ...) standardGeneric('plotDR'))

setMethod('plotDR',
          signature = 'SummExpDR',
          function(x, dim1, dim2, color_by, key, assay = NULL, ...) {
            vars_fetch <- c(dim1, dim2, color_by)
            fetched_data <- fetchData(x, varsFetch = vars_fetch, assayKey = assay, redDimKeys = key, mode = 'sample_wise')
            # renaming is to avoid weird erros when inputs have same name as arguments in functiton
            col <- color_by
            p <- scatter_plot(var1 = dim1, var2 = dim2, color_by = col, ...)
            return(p)
          })

#' Scree Plot of Explained Variance
#'
#' @inheritParams varianceExplained
#' @param dims_use dims to plot. if NULL plot all
#' @value ggplot2 object
#' @export

setGeneric('screePlot', function(x, key, feats_use = NULL, dims_use = NULL) standardGeneric('screePlot'))

setMethod('screePlot',
          signature = 'SummExpDR',
          function(x, key, feats_use = NULL, dims_use = NULL) {
            var_expl <- varianceExplained(x, key, NULL, feats_use)$r2_by_dim[dims_use]
            plot_df <- data.frame(dim = names(var_expl), pct_var = var_expl)
            p <- ggplot(data = plot_df, mapping = aes(dim, pct_var)) + geom_line(color = 'blue') +
              geom_point(size = 4.0)
            return(p)
          })

