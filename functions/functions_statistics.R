## FUNCTIONS RELATED TO FITTING AND STATISTICS
## G. J. Meijer, 06/06/2023


#' Function to generate a significance indicator for given p-values
#' 
#' @md
#' @description 
#' Generate a character string indicating the statistical significance.
#' * "***": p <= 0.001
#' * "**": 0.001 > p <= 0.01 
#' * "*": 0.01 > p <= 0.05
#' * "n.s.": p > 0.05, not significant
#' 
#' @param p array or scalar with p-values
#' @return array or scalar with character strings
#' @examples
#' p <- c(0.0002, 0.03, 0.12)
#' significance_indicator(p)
#' 
significance_indicator <- function(p) {
  ifelse(
    p <= 0.001,
    "***",
    ifelse(
      p <= 0.01,
      "**",
      ifelse(
        p <= 0.05,
        "*",
        "n.s."
      )
    )
  )
}


#' Calculate the coefficient of determination
#'
#' @description
#' Calculate the coefficient of determination
#' 
#' @param y measured values
#' @param y_predict model values
#' @param loglog if `TRUE`, log-transform `x` and `y` first
#' @return R^2 values
R2 <- function(y, y_predict, loglog = TRUE) {
  if (loglog == TRUE) {
    1 - sum((log(y) - log(y_predict))^2)/sum((log(y) - mean(log(y_predict)))^2)
  } else {
    1 - sum((y - y_predict)^2)/sum((y - mean(y_predict))^2)  
  }
}


#' Linear fit
#' 
#' @description 
#' Generate a linear fitting object. if `loglog == TRUE`, the `x` and `y`
#' values are both log-transformed prior to fitting, resulting in a power-law
#' fit
#' 
#' @param x,y arrays with x and y-values
#' @param loglog if `TRUE`, log-transform `x` and `y` prior to fitting
#' @return a `lm` fit object#' 
linear_fit <- function(x, y, loglog = TRUE) {
  if (loglog == TRUE) {
    stats::lm(log(y) ~ log(x))
  } else {
    stats::lm(y ~ x)
  }
}


#' Obtain dataframe with fitting coefficients
#' 
#' @description 
#' Generate a dataframe with key fitting results
#' 
#' @md
#' @param ft `lm` fitting object
#' @param loglog if `TRUE`, the values put into the `ft` object are assumed to
#'   have been log-transformed first.
#' @return 
#' dataframe with fields:
#' - "alpha": intercept, or if `loglog == TRUE` the power-law multiplication
#'   coefficient
#' - "beta": gradient, or if `loglog == TRUE` the power-law coefficient 
#' - "p_alpha": p-value of alpha
#' - "p_beta" p-value of beta
#' - "R2": coefficient of determination
linear_fit_coefficients <- function(ft, loglog = TRUE) {
  data.frame(
    # coefficients
    alpha = if (loglog == TRUE) {
        as.numeric(exp(ft$coef[1]))
      } else {
        as.numeric(ft$coef[1])
      },
    beta = as.numeric(ft$coef[2]),
    # p-values
    p_alpha = as.numeric(summary(ft)$coef[1, 4]),
    p_beta = as.numeric(summary(ft)$coef[2, 4]),
    # R^2 value - already logtransformed, so loglog = FALSE
    R2 = R2(ft$model[, 1], stats::predict(ft), loglog = FALSE)
  )
}


#' Generate confidence interval for linear model
#' 
#' @description 
#' Generate the best fit line and confidence interval for a fit object
#' 
#' @md
#' @inheritParams linear_fit_coefficients
#' @param level confidence level for confidence interval
#' @param n number of equally spaced points to use for fit
#' @return dataframe with
#'   - "x": x-values
#'   - "y": y-values of best
#'   - "y_upper": y-values of upper bound confidence interval
#'   - "y_lower": y-values of lower bound confidence interval
linear_fit_confidence <- function(ft, loglog = TRUE, level = 0.95, n = 101) {
  # generate range
  if (loglog == TRUE) {
    df <- data.frame(
      x = seq(exp(min(ft$model[, 2])), exp(max(ft$model[, 2])), l = n)
    )
  } else {
    df <- data.frame(
      x = seq(min(ft$model[, 2]), max(ft$model[, 2]), l = n)
    )
  }
  # generate confidence interval
  conf <- stats::predict(
    ft, 
    newdata = df, 
    interval = "confidence", 
    level = level
  )
  # return dataframe
  if (loglog == TRUE) {
    data.frame(
      x = df$x,
      y = exp(conf[, 1]),
      y_lower = exp(conf[, 2]),
      y_upper = exp(conf[, 3])
    )
  } else {
    data.frame(
      x = df$x,
      y = conf[, 1],
      y_lower = conf[, 2],
      y_upper = conf[, 3]
    )
  }
}


#' Round a number to a number of significant digits and convert to character
#' 
#' @description
#' Round a number to a number of significant digits and convert to character
#' string. Always returns `digits` number of significant digits
#' 
#' @param x numeric value to convert to character string
#' @param digits number of significant digits
#' @return character string
#' @examples
#' round_number(pi, digits = 4)
#' round_number(sqrt(5)/100000, digits = 3)
round_number <- function(x, digits = 3) {
  gsub(
    "\\.$", "", 
    formatC(
      signif(x, digits = digits), 
      digits = digits, 
      format = "fg", 
      flag = "#"
    )
  )
}


#' Generate fit annotation labels for fit
#' 
#' @description 
#' Generate a vector with character strings to show the results of the best 
#' fit and the significance of the fitting parameters.
#' 
#' It returns the alpha and beta parameters + their significance, as well
#' as the coefficient of determination `R^2`. These are generated using the
#' function `linear_fit_coefficients()`.
#' 
#' @param ft_coef dataframe with fitting coefficients, generated by the 
#'   function `linear_fit_coefficients()`
#' @param digits number of significant digits to use for each parameter
#' @param labels to use for parameters
#' @return vector with character strings
#' 
linear_fit_annotations <- function(
  ft_coef, 
  digits = 3, 
  label = c("alpha", "beta", "R^2")
) {
  # labels for fitting parameters
  lab_alpha <- paste0(
    label[1], "=='",
    round_number(ft_coef$alpha, digits), " ",
    significance_indicator(ft_coef$p_alpha), "'"
  )
  lab_beta <- paste0(
    label[2], "=='",
    round_number(ft_coef$beta, digits), " ", 
    significance_indicator(ft_coef$p_beta), "'"
  )
  lab_R2 <- paste0(
    label[3], "==", "'",
    round_number(ft_coef$R2, digits), "'"
  )
  # return vector
  c(lab_alpha, lab_beta, lab_R2)
}

