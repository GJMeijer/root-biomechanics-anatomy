## FUNCTIONS RELATED TO GENERATING PLOTS
## G. J. Meijer, 06/06/2023


#' Add a subplot annotation to an existing ggplot
#' 
#' @description 
#' Add some text annotation to an existing ggplot object
#' 
#' @param plt ggplot object
#' @param label vector with character strings to add
#' @param location two-part vector indicating the relative location in the 
#'   plot for placing labels. Follows same logic as legend placement in ggplot.
#' @param justification two-part vector indicating the horizontal and 
#'   vertical justification of the labels. Follwos the same logic as legend
#'   justification in ggplot
#' @param offset relative distance between labels and sides
#' @param label_spacing if `label` consists of multiple elements, the relative
#'   `y`-spacing between each element
#' @param size plot size of text
#' @param color plot color of text
#' @param parse if `TRUE`, labels are parsed when plotted
#' @return ggplot object with added annotation
annotate_subplot <- function(
  plt, 
  label,
  location = c(0, 1), 
  justification = c(0, 1),
  offset = 0.02,
  label_spacing = 0.10,
  size = 3,
  color = "black",
  parse = FALSE
) {
  # get current ggplot limits
  xlim <- plt$coordinates$limits$x
  ylim <- plt$coordinates$limits$y
  # plus or minus
  pm <- 1 - 2*justification
  # get coordinates for plotting
  x <- xlim[1] + (location[1] + pm[1]*offset)*diff(xlim)
  y <- ylim[1] + (location[2] + pm[2]*(offset + label_spacing*(seq(length(label)) - 1)))*diff(ylim)
  # add annotation
  plt_out <- plt + ggplot2::annotate(
    "text",
    x = x, y = y, label = label,
    hjust = justification[1], vjust = justification[2],
    color = color, size = size,
    parse = parse
  )
  # return
  return(plt_out)
}


#' Plot a ggplot object with best power-law (or linear) fit
#' 
#' @description 
#' Generate a nice-looking ggplot for linear or power-law (log-log transformed
#' linear fitting) fits
#' 
#' @param ft linear fitting object
#' @param loglog if `TRUE`, the input data in `ft` was log-log transformed 
#'   prior to fitting.
#' @param label vector with character strings to annotate the plot with. Used
#'   for showing best fit coefficients
#' @param xlab,ylab x and y-axis labels
#' @param xlim,ylim x and y-axis limits. If not defined, these are
#'   automatically chosen
#' @param color colors to use for (fitting objects, data objects)
#' @param fi
#' @return ggplot objects
#' 
linear_fit_ggplot <- function(
  ft,
  loglog = TRUE,
  label = NULL,
  xlab = "x",
  ylab = "y",
  xlim = c(NA, NA),
  ylim = c(NA, NA),
  color = RColorBrewer::brewer.pal(3, "Set1")
) {
  # points
  if (loglog == TRUE) {
    dp <- data.frame(
      x = exp(ft$model[, 2]),
      y = exp(ft$model[, 1])
    )
  } else {
    dp <- data.frame(
      x = ft$model[, 2],
      y = ft$model[, 1]
    )
  }
  # fit and confidence interval
  dc <- linear_fit_confidence(ft, loglog = loglog)
  # plot
  plt <- ggplot2::ggplot() + 
    ggplot2::theme_bw() +
    ggplot2::geom_ribbon(
      data = dc,
      ggplot2::aes(x = x, ymin = y_lower, ymax = y_upper),
      color = NA,
      fill = color[1],
      alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = dc,
      ggplot2::aes(x = x, y = y),
      color = color[1]
    ) +
    ggplot2::geom_point(
      data = dp,
      ggplot2::aes(x = x, y = y),
      color = color[2]
    ) + 
    ggplot2::xlab(xlab) + 
    ggplot2::ylab(ylab) + 
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE)
  # add annotations (if required)
  if (!is.null(label)) {
    plt <- annotate_subplot(
      plt,
      label = label,
      location = c(1, 1),
      justification = c(1, 1),
      parse = TRUE,
      color = color[1]
    )
  }
  # return
  return(plt)
}
