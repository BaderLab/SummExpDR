###############################################################################
## Plotting Utilities

filter_df <- function(df, filter_by = NULL, filter_class = NULL) {
  ## filter by value for a particular class
  ## INPUTS:
  ##    df = data.frame
  ##    filter_by = NULL or string denoting class by which data is to be subsetted
  ##    filter_class = set of valid classes in filter_by categorical variable to use in subset
  ##    Note that if filter_by is set to NULL, no filtering is done
  ##    If filter_class is left unspecified and filter_by is specified as a valid column name,
  ##    an error will be thrown
  ## RETURNS:
  ##  either the same data.frame or if valid non-NULL arguments provided to filter_by and
  ##  filter_class, a subset of that data consisting of all members of the valid classes
  ##  specified in filter_class
  if (is.character(filter_by)) {
    if (is.character(filter_class)) {
      filter_vect <- df[,filter_by]
      if (is.numeric(filter_vect)) {
        warning('specified filter feature is of numeric value, coercing to character')
      }
      filter_vect <- as.character(filter_vect)
      if (all(filter_class %in% filter_vect)) {
        valid_cell_inds <- filter_vect %in% filter_class
        df <- df[which(valid_cell_inds),]
      } else {
        stop('Not all specified classes are present in the specified filter_by column')
      }
    } else {
      warning('filter_class left unspecified, continuing to plot without filtering data')
    }
  } else if (!is.null(filter_by)) {
    warning('filter_by must be specified as string (character vector of length 1) if non-NULL\
    in order to filter. Ignoring argument')
  }
  return(df)
}

get_lims <- function(x) {
  ## Get axis limits for a set of points
  ## INPUTS:
  ##  x = numeric vector
  ## RETURNS:
  range_x <- max(x) - min(x)
  buff <- 0.05*range_x
  return(c(min(x) - buff, max(x) + buff))
}

discrete_color_mapping <- function(categorical_vals) {
  unique_vals <- unique(categorical_vals)
  unique_vals <- unique_vals[order(unique_vals)]
  col_map <- scales::hue_pal()(length(unique_vals))
  names(col_map) <- unique_vals
  return(col_map)
}

#' generic scatterplot function for data frame
#' @param df = data.frame
#' @param var1 = x axis variable
#' @param var2 = y axis variable
#' @param color_by = variable to color points by
#' @param label = variable to label samples by
#' @param filter_by = subset data by a particular categorical variable
#' @param filter_class = set of allowed classes for subset
#' @param pt.size
#' @param alpha
#' @param legend_pt_size
#' @param legend_pt_shape
#' @param xlim
#' @param ylim
#' @param lab_size
#' @value  ggplot2 object
#' @export

scatter_plot <- function(df, var1, var2, color_by, label = NULL,
                        filter_by = NULL, filter_class = NULL,
                        pt.size = 0.5, pt.shape =  16, alpha = 0.4,
                        legend_pt_size = 20, legend_pt_shape = 15,
                        xlim = NULL, ylim = NULL, legend_title = FALSE,
                        lab_size = 0.25) {
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
    p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes_(x = as.name(var1), y = as.name(var2))) +
      ggplot2::geom_point(ggplot2::aes_(color = as.name(color_by)), size = pt.size, alpha = alpha, pch = pt.shape)
    p <- p + ggplot2::scale_color_manual(values = discrete_col_mapping)
    p <- p + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=legend_pt_size,
                                                              shape=legend_pt_shape),
                                          title = ifelse(legend_title, color_by, '')))
  } else {
    # numeric values
    color_num <- df[,color_by]
    qtl_vals <- quantile(color_num[!is.na(color_num)], seq(0,1, 0.05))
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
    df$color_vals <- color_num
    color_num_not_na <- color_num[!is.na(color_num)]
    color_num_not_na[color_num_not_na >= max_val] <- max_val
    color_num_not_na[color_num_not_na <= min_val] <- min_val
    df[!is.na(color_num), 'color_vals'] <- color_num_not_na
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
    p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes_(x = as.name(var1), y = as.name(var2)))
    if (zero_center) {
      p <- p + ggplot2::geom_point(ggplot2::aes_(color = df$color_vals), size = pt.size, alpha = alpha, pch = pt.shape)
      # p <- p + geom_point(aes_(color = df$color_vals), size = pt.size, alpha = alpha, pch = pt.shape)
      # p <- p + scale_color_gradient2(midpoint = 0, limits = c(min_val, max_val), low = 'turquoise', mid = 'white', high = 'orange')
      p <- p + ggplot2::scale_color_gradient2(midpoint = 0, low = 'blue', mid = 'darkgrey', high = 'red')
      # p <- p + theme_dark()
      p <- p + ggplot2::theme_light()
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes_(color = df$color_vals), size = pt.size, alpha = alpha, pch = pt.shape)
      # p <- p + scale_color_gradient2(low = 'turquoise', mid = 'white', high = 'orange',
      #                                midpoint = mean(df$color_vals),
      #                                limits = c(min_val, max_val))
      p <- p + ggplot2::scale_color_gradient2(low = 'darkblue', mid = 'darkgrey', high = 'red', midpoint = mean(df$color_vals[!is.na(df$color_vals)]))
      # p <- p + theme_dark()
      p <- p + ggplot2::theme_light()
    }
    if (!legend_title) {
      p <- p + ggplot2::guides(color = guide_colorbar(title = NULL))
    } else {
      p <- p + ggplot2::guides(color = guide_colorbar(title = color_by))
    }
    # p <- p + scale_color_gradientn(colors = continuous_colors, breaks = continuous_breaks)

  }
  if (!is.null(label)) {
    p <- p + ggrepel::geom_text_repel(mapping = ggplot2::aes_(label = as.name(label)), size = lab_size)
  }
  p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  p <- p + ggplot2::ggtitle(color_by)
  return(p)
}

#' Plot Reduced Dims
#' @param x SummExpDR object
#' @param dim1
#' @param dim2
#' @param color_by
#' @param label
#' @param key key for reduced dims to use
#' @param assay required if pulling ouot feature values for assay data
#' @param ... other arguments to scatter_plot, with var1 and var2 overridden by dim1 and dim2
#' @value ggplot2 object
#' @export

setGeneric('plotDR', function(x, dim1, dim2, key, color_by, label = NULL, assay = NULL, ...) standardGeneric('plotDR'))

setMethod('plotDR',
          signature = 'SummExpDR',
          function(x, dim1, dim2, key, color_by, label = NULL, assay = NULL, ...) {
            vars_fetch <- c(dim1, dim2, color_by)
            if (!is.null(label)) {
              vars_fetch <- c(vars_fetch, label)
            }
            fetched_data <- fetchData(x, varsFetch = vars_fetch, assayKey = assay, redDimKeys = key, mode = 'sample_wise')
            # renaming is to avoid weird errors when inputs have same name as arguments in functiton
            col <- color_by
            lab <- label
            p <- scatter_plot(fetched_data, var1 = dim1, var2 = dim2, color_by = col, label = lab, ...)
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
            var_expl <- varianceExplained(x, key, dims_use, feats_use)$r2_by_dim
            plot_df <- data.frame(dim = names(var_expl), var_expl = var_expl, index = 1:length(var_expl))
            p <- ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(index, var_expl)) + ggplot2::geom_line(color = 'blue') +
              ggplot2::geom_point(size = 4.0, color = 'blue') + ggplot2::ggtitle(key) +
              ggplot2::xlab('Dimension') + ggplot2::ylab('% Variance Explained') +
              ggplot2::scale_x_continuous(breaks = plot_df$index)
            return(p)
          })


#' 2 variable correlation plot
#' @inheritParams scatter_plot
#' @value ggplot2 object
#' @export

corplot_2_var <- function(df, var1, var2, color_by, method = 'pearson',
                          filter_by = NULL, filter_class = NULL, pt.size = 1, alpha = 0.4,
                          legend_pt_size = 20, legend_pt_shape = 15,
                          xlim = NULL, ylim = NULL, legend_title = TRUE) {
  ## Show correlation of two numeric variables for either a dataframe or seurat object
  ## INPUTS:
  ##  df = data.frame
  ##  var1, var2 = variables to plot on x and y axis, respectively
  ##  slot = slot to pull any gene expression data from if input is a Seurat object
  ##  pt_size = point size for plot
  df_subs <- filter_df(df, filter_by = filter_by, filter_class = filter_class)
  # calculate correlation
  cor_result <- cor.test(df_subs[,var1], df_subs[,var2], method = method)
  cor_val <- round(cor_result$estimate, 2)
  cor_pval <- cor_result$p.value
  if (cor_pval < .Machine$double.eps) {
    cor_pval <- paste0(' < ', signif(.Machine$double.eps, 2))
  } else {
    cor_pval <- paste0(' = ', signif(cor_pval, 2))
  }
  # make plot
  p <- scatter_plot(df,
                   var1 = var1,
                   var2 = var2,
                   color_by = color_by,
                   filter_by = filter_by,
                   filter_class = filter_class,
                   pt.size = pt.size,
                   alpha = alpha,
                   legend_pt_size = legend_pt_size,
                   legend_pt_shape = legend_pt_shape,
                   xlim = xlim,
                   ylim = ylim,
                   legend_title = legend_title)
  title_use <- paste(var1, 'vs', var2)
  subtitle_use <- paste(method, 'correlation:', cor_val, '; p', cor_pval)
  p <- p + labs(title = title_use, subtitle = subtitle_use)
  return(p)
}

#' correlation plot method,
#' @inheritParams plotDR
#' @param feat1 x variable to plot
#' @param feat2 y variable to plot
#' @param ... other arguments to corplot_2_var
#' @value ggplot2 object
#' @export

setGeneric('plotCorrelation',
           function(x, feat1, feat2, color_by, key = NULL, assay = NULL, ...) standardGeneric('plotCorrelation'))

setMethod('plotCorrelation',
          signature = 'SummExpDR',
          function(x, feat1, feat2, color_by, key = NULL, assay = NULL, ...) {
            feats_fetch <- c(feat1, feat2, color_by)
            fetched_data <- fetchData(x, feats_fetch, key, assay, mode = 'sample_wise')
            p <- corplot_2_var(fetched_data, feat1, feat2, color_by, ...)
            return(p)
          })

#' Plot variance explained across all features and per datatype
#' @param x multiExp object
#' @param key key for reduced dims
#' @param dims_use
#' @param do_print logical, do print stuff
#' @value list of ggplot2 objects, or nothing of do_print set to TRUE
#' @export

setGeneric('plotVarianceExplained',
           function(x, key, dims_use = NULL, do_print = TRUE) standardGeneric('plotVarianceExplained'))

setMethod('plotVarianceExplained',
          signature = 'multiExp',
          function(x, key, dims_use = NULL, do_print = TRUE) {
            # get feature metadata
            row_data <- rowData(x)
            # subset feature metadata by features used in reduced dims
            loadings_mat <- getLoadings(x, key, dims_use)
            row_data <- row_data[row_data$feat_id %in% colnames(loadings_mat), ]
            # get unique experiments
            unique_expts <- sort(unique(row_data$expt))
            # variance explained across all features
            var_expl_all <- varianceExplained(x, key, dims_use)
            var_expl_mat <- matrix(var_expl_all$r2_by_dim, nrow = length(var_expl_all$r2_by_dim))
            total_var <- var_expl_all$r2_all
            # in loop, get variance explained for features in each expt
            for (expt in unique_expts) {
              var_expl_expt <- varianceExplained(x, key, dims_use, feats_use = row_data[row_data$expt == expt, 'feat_id'])
              var_expl_mat <- cbind(var_expl_mat, var_expl_expt$r2_by_dim)
              total_var <- c(total_var, var_expl_expt$r2_all)
            }
            colnames(var_expl_mat) <- c('All', unique_expts)
            rownames(var_expl_mat) <- rownames(loadings_mat[dims_use, , drop = FALSE])
            var_expl_df <- reshape2::melt(var_expl_mat)
            colnames(var_expl_df) <- c('Dim', 'DataType', 'var_expl')
            var_expl_df$var_expl <- signif(var_expl_df$var_expl, 2)
            # tile plot showing variance explained per datatype
            p1 <- ggplot2::ggplot(data = var_expl_df, mapping = ggplot2::aes(Dim, DataType)) +
              geom_tile(ggplot2::aes(fill = var_expl)) +
              scale_fill_gradient(low = 'white', high = 'red') +
              geom_text(ggplot2::aes(label = formatC(var_expl, format = "e", digits = 1)))
            # barplot showing total variance explained per datatype
            p2 <- ggplot2::ggplot(data = data.frame(total_var = total_var, features = c('All', unique_expts))) +
              ggplot2::geom_col(ggplot2::aes(x = features, y = total_var))
            if (do_print) {
              print(p1)
              print(p2)
            } else {
              return(list(p1, p2))
            }
          })

