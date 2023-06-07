## FUNCTIONS RELATED TO BIOMECHANICS MODELS
## G. J. Meijer, 06/06/2023


#' Calculate intersection point of 2 lines, based on 2 points + gradients
#' 
#' @description 
#' Calculate the intersection between two straight lines. The first line has
#' a gradient `dydx1` and passes through the point (`x1`, `y1`). The second line
#' has gradient `dydx2` and passes through the point (`x2`, `y2`).
#' 
#' @param x1,y1,dydx1 coordinates and gradient of line 1
#' @param x2,y2,dydx2 coordinates and gradient of line 2
#' @return dataframe with intersection coordinates `x` and `y`
#' @examples
#' intersection_coordinates(0, 0, 2, 5, 4, 0)
intersection_coordinates <- function(x1, y1, dydx1, x2, y2, dydx2) {
  xi <- (y2 - y1 - x2*dydx2 + x1*dydx1)/(dydx1 - dydx2)
  yi <- y1 + (xi - x1)*dydx1
  data.frame(x = xi, y = yi)
}


#' Fit tensile stress--strain trace with linear segments
#' 
#' @md
#' @description
#' Fit a measured tensile stress--strain curve with three linear segments:
#' 1. Initial stiffness
#' 2. Elastic stiffness
#' 3. Plastic stiffness
#' 
#' ELASTIC STIFFNESS
#' The elastic stiffness is defined as the largest gradient dy/dx along the
#' trace. To reduce the effect of noise in the data, the data is smoothed
#' using a spline smoothing function with a user-selected value of the `spar`
#' parameter. Larger values of `spar` result in smoother traces 
#' (`smooth.spline(x, y, spar = spar)`). This procedure gives the location 
#' along the curve where the gradient is largest (`x1`, `y1`) as well as the 
#' gradient itself (`dydx1`).
#' 
#' PLASTIC STIFFNESS
#' To obtain the plastic stiffness, we first determine the point `x2`,`y2` 
#' where the `y` value is largest. This is the peak reinforcement. Subsequently,
#' we find a point (xi,yi) on the domain 
#' [x1 <= x <= x1 + (1 - f)*(y2 - y1)/dydx1 + f*(x2 - x1)]
#' that minimises the secant stiffness: (y2 - yi)/(x2 - xi). This is considered
#' the plastic stiffness
#' 
#' INITIAL STIFFNESS
#' The initial part of the curve (`min(x)` <= `x` <= `x1`) is fitted with 
#' another smoothing spline. This time, the number of knots is selected by the 
#' user. The gradient at `x = min(x)` is taken as the initial stiffness of the 
#' curve. We thus obtain another point `x0 = min(x)`, `y0` and a gradient 
#' `dxdy0`.
#' 
#' It is assumed that the force is zeroed before the root is installed. 
#' If the force y > 0 at x == 0, the root is under tension prior to 
#' starting the test. If y < 0 at x == 0, the root is under compression directly
#' after placing the root in the test equipment.
#' 
#' The function return a list with the x,y coordinates for the following points:
#' - 1: Intersection initial stiffness and y = 0
#' - 2: Point where initial stiffness determined
#' - 3: Intersection elastic stiffness and y = 0
#' - 4: Intersection initial and elastic stiffness
#' - 5: Point where elastic stiffness determined
#' - 6: Intersection elastic and plastic stiffness
#' - 7: Point where plastic stiffness determined
#' - 8: Point at end of plastic search window
#' - 9: Ultimate point
#' 
#' @param x,y arrays with measured forces and displacements
#' @param spar spar argument used in smooth spline fitting `smooth.spline(x, y)`
#' @param nknots number of knots for spline fitting of initial section of 
#'   curve (min(x) <= x <= x5)
#' @param n number of points to use on x-domain for spline predictions
#' @param f additional area for plastic stiffness search domain
#' @return dataframe with all eight fitting 
#'   points:
#'   - `displacement`: x-positions
#'   - `force`: y-positions
#'   - `gradient`: force/displacement gradient at point
#'   
fit_tensiletest_linear <- function(
  x, y, 
  f = 1/3,
  spar = 0.75, 
  nknots = 4,
  n = 1001
) {
  # find maximum force value - ultimate point (point 8)
  i9 <- which.max(y)
  x9 <- x[i9]
  y9 <- y[i9]
  # fit smoothing spline to data - up to x=x_u
  spl <- stats::smooth.spline(x[1:i9], y[1:i9], spar = spar)
  xs <- seq(min(x), x9, l = n)
  ys <- stats::predict(spl, xs)$y
  dysdxs <- stats::predict(spl, x = xs, deriv = 1)$y
  # find 'elastic' stiffness - point where dy/dx is maximum (point 5)
  i5 <- which.max(dysdxs)
  x5 <- xs[i5]
  y5 <- predict(spl, x = x5)$y
  dysdxs5 <- dysdxs[i5]
  # 'plastic' stiffness
  x8 <- x5 + (1 - f)*(y9 - y5)/dysdxs5 + f*(x9 - x5)
  sel <- (xs >= x5) & (xs <= x8)
  i7 <- which.min((ys[sel] - y9)/(xs[sel] - x9))
  x7 <- xs[sel][i7]
  y7 <- ys[sel][i7]
  dysdxs9 <- (y9 - y7)/(x9 - x7)
  y8 <- y9 - (x9 - x8)*dysdxs9
  # 'initial' stiffness - fit smoothing spline to only first part of the curve
  ind <- (x <= x5)
  if (length(unique(x[ind])) < 4) {
    ind[1:(2*4 - length(unique(x[ind])))] <- TRUE
  }
  spl0 <- stats::smooth.spline(x[ind], y[ind], nknots = nknots)
  xs0 <- seq(min(x), x5, l = n)
  ys0 <- stats::predict(spl0, x = xs0)$y
  range2 <- 0.05  # average initial stiffness over the first 5% of points between the initial and elastic point
  if (min(ys0) >= 0) {
    # root initially under tension - or neutral
    dysdxs2 <- mean(stats::predict(spl0, x = xs0[1:ceiling(range2*n)], deriv = 1)$y)
    x2 <- xs0[1]
    y2 <- stats::predict(spl0, x = x2)$y
  } else {
    # root initially under compression
    i0 <- max(which(ys0 < 0))
    y2 <- 0
    x2 <- stats::approx(ys0[i0:n], xs0[i0:n], xout = y2)$y
    xleft <- max(min(xs0), x2 - 0.5*range2*(max(xs0) - min(xs0)))
    xright <- min(max(xs0), x2 + 0.5*range2*(max(xs0) - min(xs0)))
    xav <- xs0[(xs0 > xleft) & (xs0 <= xright)]
    dysdxs2 <- mean(stats::predict(spl0, x = xav, deriv = 1)$y)
  }
  # intersection initial stiffness and y=0
  xy1 <- intersection_coordinates(0, 0, 0, x2, y2, dysdxs2)
  x1 <- xy1$x
  y1 <- xy1$y
  # intersection elastic stiffness and y=0
  xy3 <- intersection_coordinates(0, 0, 0, x5, y5, dysdxs5)
  x3 <- xy3$x
  y3 <- xy3$y
  # intersection initial and elastic stiffness
  xy4 <- intersection_coordinates(x2, y2, dysdxs2, x5, y5, dysdxs5)
  x4 <- xy4$x
  y4 <- xy4$y
  # intersection elastic and plastic stiffness
  xy6 <- intersection_coordinates(x5, y5, dysdxs5, x9, y9, dysdxs9)
  x6 <- xy6$x
  y6 <- xy6$y
  # return dataframe with key values
  data.frame(
    point = seq(9),
    displacement = c(x1, x2, x3, x4, x5, x6, x7, x8, x9),
    force = c(y1, y2, y3, y4, y5, y6, y7, y8, y9),
    gradient = c(NA, dysdxs2, NA, NA, dysdxs5, NA, NA, NA, dysdxs9)
  )
}


#' Plot tensile test fitting process
#' 
#' @description
#' Generates a ggplot to clarify the concepts behind the root tensile test
#' linear fitting process. Fitting process is determined by function
#' `fit_tensiletest_linear()`.
#' 
#' @inheritParams fit_tensiletest_linear
#' @param xf,yf vectors with x and y positions of key fitting points, 
#'   determined using the function `fit_tensiletest_linear()`
#' @param color 4-element vector with colors for: data, elastic, plastic
#'   and initial curve
#' @param line 4-element vector with linetypes for: data, elastic, plastic and
#'   initial curve
#' @param label_color color for annotations
#' @param label_size font size for annotations
#' @param label_linesize line thickness for annotations
#' @param alpha transparancy value for plastic search zone
#' @return ggplot object
#' 
plot_tensiletest_linear <- function(
  x, y,
  xf, yf,
  color = c("black","#E41A1C", "#377EB8", "#4DAF4A"),
  line = c(1, 2, 4, 5),
  label_color = "grey40",
  label_size = 2.5,
  label_linesize = 0.3,
  alpha = 0.2
) {
  # plot
  ggplot2::ggplot() + 
    # set default theme
    ggplot2::theme_bw() + 
    # plastic stiffness zone
    ggplot2::geom_rect(
      ggplot2::aes(xmin = xf[5], xmax = xf[8], ymin = -Inf, ymax = Inf, fill = "Plastic"),
      color = NA, alpha = alpha, show.legend = FALSE
    ) +
    # data traces - measured
    ggplot2:: geom_path(aes(x = x, y = y, color = "Data", linetype = "Data")) +
    # plastic stiffness
    ggplot2::geom_point(ggplot2::aes(x = xf[7], y = yf[7], color = "Plastic")) + 
    ggplot2::geom_point(ggplot2::aes(x = xf[9], y = yf[9], color = "Plastic")) + 
    ggplot2::geom_segment(ggplot2::aes(
      x = xf[5], 
      y = yf[9] - (xf[9] - xf[5])*(yf[9] - yf[6])/(xf[9] - xf[6]), 
      xend = xf[9], 
      yend = yf[9],
      color = "Plastic", linetype = "Plastic"
    )) + 
    # initial stiffness
    ggplot2::geom_point(ggplot2::aes(x = xf[2], y = yf[2], color = "Initial")) + 
    ggplot2::geom_segment(ggplot2::aes(
      x = xf[1], y = yf[1], 
      xend = xf[6], 
      yend = yf[1] + (xf[6] - xf[1])*(yf[4] - yf[1])/(xf[4] - xf[1]),
      color = "Initial", linetype = "Initial"
    )) + 
    # elastic stiffness
    ggplot2::geom_point(ggplot2::aes(x = xf[5], y = yf[5], color = "Elastic")) + 
    ggplot2::geom_segment(ggplot2::aes(
      x = xf[3], y = yf[3], 
      xend = xf[3] + (yf[9] - yf[3])*(xf[5] - xf[3])/(yf[5] - yf[3]), 
      yend = yf[9],
      color = "Elastic", linetype = "Elastic"
    )) + 
    # annotate text and arrows for yield point
    ggplot2::annotate(
      "text",
      x = xf[5],
      y = yf[6] + 0.5*(yf[9] - yf[6]),
      label = "Yield\npoint",
      vjust = 0, hjust = 0.5,
      size = label_size,
      color = label_color
    ) +
    ggplot2::geom_curve(
      ggplot2::aes(x = xf[5], y = yf[6] + 0.5*(yf[9] - yf[6]), xend = xf[6], yend = yf[6]),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.08, "inch")), 
      size = label_linesize,
      color = label_color, 
      curvature = 0.3
    ) +
    # annotate text and arrows for ultimate point
    ggplot2::annotate(
      "text",
      x = xf[8] + 0.7*(xf[9] - xf[8]),
      y = 2/3*yf[9],
      label = "Ultimate\npoint",
      vjust = 1, hjust = 0.5,
      size = label_size,
      color = label_color
    ) +
    ggplot2::geom_curve(
      aes(x = xf[8] + 0.7*(xf[9] - xf[8]), y = 2/3*yf[9], xend = xf[9], yend = yf[9]),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.08, "inch")), 
      linewidth = label_linesize,
      color = label_color, 
      curvature = -0.2
    ) +
    # annotate u_{r,t}
    ggplot2::annotate(
      "segment",
      x = xf[1], y = 0.5*yf[5], xend = xf[3], yend = 0.5*yf[5],
      color = label_color,
      arrow = arrow(length = unit(0.08, "inch"), ends = "both"),
      linewidth = label_linesize
    ) + 
    ggplot2::annotate(
      "text",
      x = 0.5*(xf[1] + xf[3]), y = 0.52*yf[5],
      label = "u[r*','*t]",
      hjust = 0.5, vjust = 0,
      parse = TRUE,
      size = label_size, color = label_color
    ) + 
    ggplot2::annotate(
      "segment",
      x = xf[1], y = 0, xend = xf[1], yend = 0.6*yf[5],
      color = label_color,
      linetype = 6,
      linewidth = label_linesize
    ) + 
    ggplot2::annotate(
      "segment",
      x = xf[3], y = 0, xend = xf[3], yend = 0.6*yf[5],
      color = label_color,
      linetype = 6,
      linewidth = label_linesize
    ) +
    # axis labels
    ggplot2::xlab(expression("Displacement"~u[r]~"[mm]")) + 
    ggplot2::ylab(expression("Tensile force"~T[r]~"[N]")) + 
    # axes limits
    ggplot2::coord_cartesian(xlim = c(min(0, xf[1]), NA), ylim = c(0, NA)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, .05))) + 
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .05))) +
    # colors
    ggplot2::scale_color_manual(name = NULL, values = color) +
    ggplot2::scale_linetype_manual(name = NULL, values = line) +
    ggplot2::scale_fill_manual(name = NULL, values = color[4])
}


#' Convert root properties to cortex/stele properties
#' 
#' @description
#' Take root yield and ultimate stress and strain, and convert this to 
#' cortex and stele properties using known stele/cortex areas.
#' 
#' The stele--cortex interface is assumed as strong (so it will never debond) 
#' in `model = 1`, while the interface is extremely weak (`r_i = 0`) when
#' `model = 2`.
#' 
#' @param t_ru,eps_ru root peak strength and strain
#' @param t_ry,eps_ry root yield strength and strain
#' @param A_r root cross-sectional area
#' @param A_s stele cross-sectional area
#' @param A_c (optional) cortex cross-sectiona area. By default, 
#'   `A_c = A_r - A_s`
#' @param model either `model = 1` or `model = 2`, see Description
#' @return a dataframe with stele properties `t_su`, `eps_su` and cortex
#'   properties `t_cu`, `eps_cu`
#' @examples
#' # define some parameters
#' t_ry <- 5
#' eps_ry <- 0.05
#' t_ru <- 10
#' eps_ru <- 0.20
#' A_r <- 2
#' A_s <- 0.5
#' 
#' # calculate
#' convert_root2biomech(t_ru, eps_ru, t_ry, eps_ry, A_r, A_s, model = 1)
#' convert_root2biomech(t_ru, eps_ru, t_ry, eps_ry, A_r, A_s, model = 2)
convert_root2biomech <- function(
  t_ru, eps_ru, 
  t_ry, eps_ry, 
  A_r, A_s, A_c = A_r - A_s, 
  model = 1
) {
  if (model == 1) {
    r_eps <- eps_ry*t_ry^2*sqrt(2*t_ru - t_ry) / 
      (t_ru*(t_ry*eps_ru - t_ru*eps_ry)*sqrt(t_ry)*acosh(t_ru/(t_ru - t_ry)) + 
         eps_ry*t_ry*t_ru*sqrt(2*t_ru - t_ry))
  } else if (model == 2) {
    r_eps <- eps_ry/eps_ru
  } else {
    warning("No available model selected")
  }
  r_T <- t_ry/t_ru - r_eps
  data.frame(
    t_su = (A_r/A_s)*t_ru,
    eps_su = eps_ry/r_eps,
    t_cu = r_T*(A_r/A_c)*t_ru,
    eps_cu = eps_ry
  )
}


#' Convert cortex/stele properties to root properties
#' 
#' @description
#' Take stele and cortex properties and convert to root properies (yield point,
#' ultimate point) using known stele/cortex areas.
#' 
#' The stele--cortex interface is assumed as strong (so it will never debond) 
#' in `model = 1`, while the interface is extremely weak (`r_i = 0`) when
#' `model = 2`.
#' 
#' @inheritParams convert_root2biomech
#' @return a dataframe with root ultimate properties `t_ru`, `eps_ru` and
#'   yield properties `t_ry`, `eps_ry`
#' @examples
#' # define some root parameters
#' t_ry <- 5
#' eps_ry <- 0.05
#' t_ru <- 10
#' eps_ru <- 0.20
#' A_r <- 2
#' A_s <- 0.5
#' 
#' # calculate biomechanical properties
#' df1 <- convert_root2biomech(t_ru, eps_ru, t_ry, eps_ry, A_r, A_s, model = 1)
#' df2 <- convert_root2biomech(t_ru, eps_ru, t_ry, eps_ry, A_r, A_s, model = 2)
#' 
#' # check if function correct by reverting back to root properties
#' with(df1, convert_biomech2root(t_su, eps_su, t_cu, eps_cu, A_r, A_s, model = 1))
#' with(df2, convert_biomech2root(t_su, eps_su, t_cu, eps_cu, A_r, A_s, model = 2))
convert_biomech2root <- function(
  t_su, eps_su, 
  t_cu, eps_cu, 
  A_r, A_s, A_c = A_r - A_s, 
  model = 1
) {
  if (model == 1) {
    r_T <- (A_c*t_cu)/(A_s*t_su)
    r_eps <- eps_cu/eps_su
    eps_ru <- eps_su*(
      r_eps/(r_T + r_eps) + r_T*sqrt(2 - r_T - r_eps)/
        (sqrt(r_T + r_eps)*acosh(1/(1 - r_T - r_eps))))
  } else if (model == 2) {
    eps_ru <- eps_su
  } else {
    warning("No available model selected")
  }
  data.frame(
    t_ry = (A_s/A_r)*(eps_cu/eps_su)*t_su + (A_c/A_r)*t_cu,
    eps_ry = eps_cu,
    t_ru = (A_s/A_r)*t_su,
    eps_ru = eps_ru
  )
}

