# load packages
#library("tidyverse")
library("dplyr")
library("readr")
library("stringr")
library("magrittr")
library("ggplot2")
library("gridExtra")

# load functions
source("functions/functions_statistics.R")
source("functions/functions_biomech.R")
source("functions/functions_plotting.R")

# switches
plot_save <- TRUE  # plot various fits, and save to file
data_save <- TRUE   # save R^2 for each model to a file


## LOAD INTERPRETED DATA -----------------------------------------------------

# load data
df <- read_csv("data_interpreted/data_rootproperties.csv") %>%
  rename_with(~str_trim(str_remove(.x, "\\[.*\\]")))  # trim units from headers
# calculate diameters
df <- df %>%
  mutate(
    d_r = sqrt(4/pi*A_r),
    d_s = sqrt(4/pi*A_s),
    d_c = 0.5*(d_r - d_s)
  )


## INITIATE FITTING DATAFRAME ----------------------------------------

# initiate dataframe for all fits
dfit <- data.frame(
  model = NULL, t_ru = NULL, eps_ru = NULL, t_ry = NULL, eps_ry = NULL
)


## MODEL 0a ------------------------------------------------------------

# Model 0a - fit as function of root diameter
ft0a_tru_dr <- linear_fit(df$d_r, df$t_ru, loglog = TRUE)
ft0a_epsru_dr <- linear_fit(df$d_r, df$eps_ru, loglog = TRUE)
ft0a_try_dr <- linear_fit(df$d_r, df$t_ry, loglog = TRUE)
ft0a_epsry_dr <- linear_fit(df$d_r, df$eps_ry, loglog = TRUE)
# fit coefficients
coef0a_tru_dr <- linear_fit_coefficients(ft0a_tru_dr)
coef0a_epsru_dr <- linear_fit_coefficients(ft0a_epsru_dr)
coef0a_try_dr <- linear_fit_coefficients(ft0a_try_dr)
coef0a_epsry_dr <- linear_fit_coefficients(ft0a_epsry_dr)
# add R^2 to results dataframe
dfit <- bind_rows(dfit, data.frame(
  model = "0a", 
  t_ru = coef0a_tru_dr$R2, eps_ru = coef0a_epsru_dr$R2, 
  t_ry = coef0a_try_dr$R2, eps_ry = coef0a_epsry_dr$R2
))


## MODEL 0b ------------------------------------------------------------

# Model 0b - fit as function of stele diameter
ft0b_tru_ds <- with(df, linear_fit(d_s, t_ru, loglog = TRUE))
ft0b_epsru_ds <- with(df, linear_fit(d_s, eps_ru, loglog = TRUE))
ft0b_try_ds <- with(df, linear_fit(d_s, t_ry, loglog = TRUE))
ft0b_epsry_ds <- with(df, linear_fit(d_s, eps_ry, loglog = TRUE))
# fit coefficients
coef0b_tru_ds <- linear_fit_coefficients(ft0b_tru_ds)
coef0b_epsru_ds <- linear_fit_coefficients(ft0b_epsru_ds)
coef0b_try_ds <- linear_fit_coefficients(ft0b_try_ds)
coef0b_epsry_ds <- linear_fit_coefficients(ft0b_epsry_ds)
# add R^2 to results dataframe
dfit <- bind_rows(dfit, data.frame(
  model = "0b", 
  t_ru = coef0b_tru_ds$R2, eps_ru = coef0b_epsru_ds$R2, 
  t_ry = coef0b_try_ds$R2, eps_ry = coef0b_epsry_ds$R2
))


## MODEL 1 - NO SLIPPAGE -----------------------------------------------------

# obtain biomechanical properties
df1 <- bind_cols(
  df, 
  with(df, convert_root2biomech(t_ru, eps_ru, t_ry, eps_ry, A_r = A_r, A_s = A_s, model = 1))
)

# Model 1a - fit as function of root diameter
ft1a_tsu_dr <- with(df1, linear_fit(d_r, t_su, loglog = TRUE))
ft1a_epssu_dr <- with(df1, linear_fit(d_r, eps_su, loglog = TRUE))
ft1a_tcu_dr <- with(df1, linear_fit(d_r, t_cu, loglog = TRUE))
ft1a_epscu_dr <- with(df1, linear_fit(d_r, eps_cu, loglog = TRUE))
# predict
tsu_1a <- exp(predict(ft1a_tsu_dr))
epssu_1a <- exp(predict(ft1a_epssu_dr))
tcu_1a <- exp(predict(ft1a_tcu_dr))
epscu_1a <- exp(predict(ft1a_epscu_dr))
# convert back to root predictions
df1a <- with(df1, convert_biomech2root(tsu_1a, epssu_1a, tcu_1a, epscu_1a, A_r = A_r, A_s = A_s, model = 1))
# add R^2 to results dataframe
dfit <- bind_rows(dfit, data.frame(
  model = "1a",
  t_ru = R2(df1$t_ru, df1a$t_ru, loglog = TRUE),
  eps_ru = R2(df1$eps_ru, df1a$eps_ru, loglog = TRUE),
  t_ry = R2(df1$t_ry, df1a$t_ry, loglog = TRUE),
  eps_ry = R2(df1$eps_ry, df1a$eps_ry, loglog = TRUE)
))

# Model 1b - fit as function of stele/cortex diameter
ft1b_tsu_ds <- with(df1, linear_fit(d_s, t_su, loglog = TRUE))
ft1b_epssu_ds <- with(df1, linear_fit(d_s, eps_su, loglog = TRUE))
ft1b_tcu_dc <- with(df1, linear_fit(d_c, t_cu, loglog = TRUE))
ft1b_epscu_dc <- with(df1, linear_fit(d_c, eps_cu, loglog = TRUE))
# fit coefficients
coef1b_tsu_ds <- linear_fit_coefficients(ft1b_tsu_ds)
coef1b_epssu_ds <- linear_fit_coefficients(ft1b_epssu_ds)
coef1b_tcu_dc <- linear_fit_coefficients(ft1b_tcu_dc)
coef1b_epscu_dc <- linear_fit_coefficients(ft1b_epscu_dc)
# predict
tsu_1b <- exp(predict(ft1b_tsu_ds))
epssu_1b <- exp(predict(ft1b_epssu_ds))
tcu_1b <- exp(predict(ft1b_tcu_dc))
epscu_1b <- exp(predict(ft1b_epscu_dc))
# convert back to root predictions
df1b <- with(df1, convert_biomech2root(tsu_1b, epssu_1b, tcu_1b, epscu_1b, A_r = A_r, A_s = A_s, model = 1))
# add R^2 to results dataframe
dfit <- bind_rows(dfit, data.frame(
  model = "1b",
  t_ru = R2(df1$t_ru, df1b$t_ru, loglog = TRUE),
  eps_ru = R2(df1$eps_ru, df1b$eps_ru, loglog = TRUE),
  t_ry = R2(df1$t_ry, df1b$t_ry, loglog = TRUE),
  eps_ry = R2(df1$eps_ry, df1b$eps_ry, loglog = TRUE)
))


## MODEL 2 - FULL SLIPPAGE -----------------------------------------------------

# obtain biomechanical properties
df2 <- bind_cols(
  df, 
  with(df, convert_root2biomech(t_ru, eps_ru, t_ry, eps_ry, A_r = A_r, A_s = A_s, model = 2))
)

# Model 2a - fit as function of root diameter
ft2a_tsu_dr <- with(df2, linear_fit(d_r, t_su, loglog = TRUE))
ft2a_epssu_dr <- with(df2, linear_fit(d_r, eps_su, loglog = TRUE))
ft2a_tcu_dr <- with(df2, linear_fit(d_r, t_cu, loglog = TRUE))
ft2a_epscu_dr <- with(df2, linear_fit(d_r, eps_cu, loglog = TRUE))
# predict
tsu_2a <- exp(predict(ft2a_tsu_dr))
epssu_2a <- exp(predict(ft2a_epssu_dr))
tcu_2a <- exp(predict(ft2a_tcu_dr))
epscu_2a <- exp(predict(ft2a_epscu_dr))
# convert back to root predictions
df2a <- with(df2, convert_biomech2root(tsu_2a, epssu_2a, tcu_2a, epscu_2a, A_r = A_r, A_s = A_s, model = 2))
# add R^2 to results dataframe
dfit <- bind_rows(dfit, data.frame(
  model = "2a",
  t_ru = R2(df2$t_ru, df2a$t_ru, loglog = TRUE),
  eps_ru = R2(df2$eps_ru, df2a$eps_ru, loglog = TRUE),
  t_ry = R2(df2$t_ry, df2a$t_ry, loglog = TRUE),
  eps_ry = R2(df2$eps_ry, df2a$eps_ry, loglog = TRUE)
))

# Model 2b - fit as function of stele/cortex diameter
ft2b_tsu_ds <- with(df2, linear_fit(d_s, t_su, loglog = TRUE))
ft2b_epssu_ds <- with(df2, linear_fit(d_s, eps_su, loglog = TRUE))
ft2b_tcu_dc <- with(df2, linear_fit(d_c, t_cu, loglog = TRUE))
ft2b_epscu_dc <- with(df2, linear_fit(d_c, eps_cu, loglog = TRUE))
# fit coefficients
coef2b_tsu_ds <- linear_fit_coefficients(ft2b_tsu_ds)
coef2b_epssu_ds <- linear_fit_coefficients(ft2b_epssu_ds)
coef2b_tcu_dc <- linear_fit_coefficients(ft2b_tcu_dc)
coef2b_epscu_dc <- linear_fit_coefficients(ft2b_epscu_dc)
# predict
tsu_2b <- exp(predict(ft2b_tsu_ds))
epssu_2b <- exp(predict(ft2b_epssu_ds))
tcu_2b <- exp(predict(ft2b_tcu_dc))
epscu_2b <- exp(predict(ft2b_epscu_dc))
# convert back to root predictions
df2b <- with(df2, convert_biomech2root(tsu_2b, epssu_2b, tcu_2b, epscu_2b, A_r = A_r, A_s = A_s, model = 2))
# add R^2 to results dataframe
dfit <- bind_rows(dfit, data.frame(
  model = "2b",
  t_ru = R2(df2$t_ru, df2b$t_ru, loglog = TRUE),
  eps_ru = R2(df2$eps_ru, df2b$eps_ru, loglog = TRUE),
  t_ry = R2(df2$t_ry, df2b$t_ry, loglog = TRUE),
  eps_ry = R2(df2$eps_ry, df2b$eps_ry, loglog = TRUE)
))


## SAVE FITTING RESULTS --------------------------------------------

# save to file
if (data_save == TRUE) {
  write_csv(
    dfit,
    "data_interpreted/results_models_R2.csv"
  )
}


## DEFINE LABELS AND RANGES FOR PLOTTING ------------------------------------

# labels
label_dr <- expression("Root diameter"~d[r]~"[mm]")
label_ds <- expression("Stele diameter"~d[s]~"[mm]")
label_dc <- expression("Cortex thickness"~d[c]~"[mm]")
label_try <- expression("Yield strength"~t[r*","*y]~"[MPa]")
label_epsry <- expression("Yield strain"~epsilon[r*","*y]~"[mm/mm]")
label_tru <- expression("Root strength"~t[r*","*u]~"[MPa]")
label_epsru <- expression("Root peak strain"~epsilon[r*","*u]~"[mm/mm]")
label_tsu <- expression("Stele strength"~t[s*","*u]~"[MPa]")
label_epssu <- expression("Stele peak strain"~epsilon[s*","*u]~"[mm/mm]")
label_tcu <- expression("Cortex strength"~t[c*","*u]~"[MPa]")
label_epscu <- expression("Cortex peak strain"~epsilon[c*","*u]~"[mm/mm]")
label_Er <- expression("Root stiffness"~E[r]~"[MPa]")
label_Es <- expression("Stele stiffness"~E[s]~"[MPa]")
label_Ec <- expression("Cortex stiffness"~E[c]~"[MPa]")
# ranges
range_dr <- c(0, 4)
range_ds <- c(0, 2.5)
range_dc <- c(0, 0.8)
range_try <- c(0, 12)
range_epsry <- c(0, 0.15)
range_tru <- c(0, 20)
range_epsru <- c(0, 0.40)
range_tsu <- c(0, 100)
range_epssu <- c(0, 0.40)
range_tcu <- c(0, 12)
range_epscu <- c(0, 0.15)
range_Er <- c(0, 800)
range_Es <- c(0, 800)
range_Ec <- c(0, 800)


## generate plots
if (plot_save == TRUE) {
  
  ## GENERATE FIT PLOTS -------------------------------------------------------
  
  # model 0a
  p0a_tru_dr <- linear_fit_ggplot(
    ft0a_tru_dr, loglog = TRUE, 
    xlab = label_dr, ylab = label_tru,
    xlim = range_dr, ylim = range_tru,
    label = linear_fit_annotations(coef0a_tru_dr)
  )
  p0a_epsru_dr <- linear_fit_ggplot(
    ft0a_epsru_dr, loglog = TRUE, 
    xlab = label_dr, ylab = label_epsru,
    xlim = range_dr, ylim = range_epsru,
    label = linear_fit_annotations(coef0a_epsru_dr)
  )
  p0a_try_dr <- linear_fit_ggplot(
    ft0a_try_dr, loglog = TRUE, 
    xlab = label_dr, ylab = label_try,
    xlim = range_dr, ylim = range_try,
    label = linear_fit_annotations(coef0a_try_dr)
  )
  p0a_epsry_dr <- linear_fit_ggplot(
    ft0a_epsry_dr, loglog = TRUE, 
    xlab = label_dr, ylab = label_epsry,
    xlim = range_dr, ylim = range_epsry,
    label = linear_fit_annotations(coef0a_epsry_dr)
  )
  # model 1b
  p1b_tsu_ds <- linear_fit_ggplot(
    ft1b_tsu_ds, loglog = TRUE, 
    xlab = label_ds, ylab = label_tsu,
    xlim = range_ds, ylim = range_tsu,
    label = linear_fit_annotations(coef1b_tsu_ds)
  )
  p1b_epssu_ds <- linear_fit_ggplot(
    ft1b_epssu_ds, loglog = TRUE, 
    xlab = label_ds, ylab = label_epssu,
    xlim = range_ds, ylim = range_epssu,
    label = linear_fit_annotations(coef1b_epssu_ds)
  )
  p1b_tcu_dc <- linear_fit_ggplot(
    ft1b_tcu_dc, loglog = TRUE, 
    xlab = label_dc, ylab = label_tcu,
    xlim = range_dc, ylim = range_tcu,
    label = linear_fit_annotations(coef1b_tcu_dc)
  )
  p1b_epscu_dc <- linear_fit_ggplot(
    ft1b_epscu_dc, loglog = TRUE, 
    xlab = label_dc, ylab = label_epscu,
    xlim = range_dc, ylim = range_epscu,
    label = linear_fit_annotations(coef1b_epscu_dc)
  )
  # model 2b
  p2b_tsu_ds <- linear_fit_ggplot(
    ft2b_tsu_ds, loglog = TRUE, 
    xlab = label_ds, ylab = label_tsu,
    xlim = range_ds, ylim = range_tsu,
    label = linear_fit_annotations(coef2b_tsu_ds)
  )
  p2b_epssu_ds <- linear_fit_ggplot(
    ft2b_epssu_ds, loglog = TRUE, 
    xlab = label_ds, ylab = label_epssu,
    xlim = range_ds, ylim = range_epssu,
    label = linear_fit_annotations(coef2b_epssu_ds)
  )
  p2b_tcu_dc <- linear_fit_ggplot(
    ft2b_tcu_dc, loglog = TRUE, 
    xlab = label_dc, ylab = label_tcu,
    xlim = range_dc, ylim = range_tcu,
    label = linear_fit_annotations(coef2b_tcu_dc)
  )
  p2b_epscu_dc <- linear_fit_ggplot(
    ft2b_epscu_dc, loglog = TRUE, 
    xlab = label_dc, ylab = label_epscu,
    xlim = range_dc, ylim = range_epscu,
    label = linear_fit_annotations(coef2b_epscu_dc)
  )
  
  # merge subplots together into single plot
  grd <- arrangeGrob(
    annotate_subplot(p0a_tru_dr, "a)"),
    annotate_subplot(p1b_tsu_ds, "b)"),
    annotate_subplot(p2b_tsu_ds, "c)"),
    annotate_subplot(p0a_epsru_dr, "d)"),
    annotate_subplot(p1b_epssu_ds, "e)"),
    annotate_subplot(p2b_epssu_ds, "f)"),
    annotate_subplot(p0a_try_dr, "g)"),
    annotate_subplot(p1b_tcu_dc, "h)"),
    annotate_subplot(p2b_tcu_dc, "i)"),
    annotate_subplot(p0a_epsry_dr, "j)"),
    annotate_subplot(p1b_epscu_dc, "k)"),
    annotate_subplot(p2b_epscu_dc, "l)"),
    ncol = 3
  )
  #save_manually - height 9", width 7"
  ggsave(
    "plots/plot_model_fits.pdf",
    grd,
    width = 7, height = 9, units = "in"
  )
  
  
  ## GENERATE STIFFNESS PLOTS -------------------------------
  
  # generate fits - assume model 0a, 1b and 2b
  ft0_Er_dr <- with(df, linear_fit(d_r, t_ry/eps_ry, loglog = TRUE))
  ft2_Es_ds <- with(df2, linear_fit(d_s, t_su/eps_su, loglog = TRUE))
  ft2_Ec_dc <- with(df2, linear_fit(d_c, t_cu/eps_cu, loglog = TRUE))
  # coefficients
  coef0_Er_dr <- linear_fit_coefficients(ft0_Er_dr, loglog = TRUE)
  coef2_Es_ds <- linear_fit_coefficients(ft2_Es_ds, loglog = TRUE)
  coef2_Ec_dc <- linear_fit_coefficients(ft2_Ec_dc, loglog = TRUE)
  # generate plots  
  p0_Er_dr <- linear_fit_ggplot(
    ft0_Er_dr, loglog = TRUE, 
    xlab = label_dr, ylab = label_Er,
    xlim = range_dr, ylim = range_Er,
    label = linear_fit_annotations(coef0_Er_dr)
  )
  p2_Es_ds <- linear_fit_ggplot(
    ft2_Es_ds, loglog = TRUE, 
    xlab = label_ds, ylab = label_Es,
    xlim = range_ds, ylim = range_Es,
    label = linear_fit_annotations(coef2_Es_ds)
  )
  p0_Ec_dc <- linear_fit_ggplot(
    ft2_Ec_dc, loglog = TRUE, 
    xlab = label_dc, ylab = label_Ec,
    xlim = range_dc, ylim = range_Ec,
    label = linear_fit_annotations(coef2_Ec_dc)
  )
  # arrange together
  grd_E <- arrangeGrob(
    annotate_subplot(p0_Er_dr, "a)"),
    annotate_subplot(p2_Es_ds, "b)"),
    annotate_subplot(p0_Ec_dc, "c)"),
    ncol = 3
  )
  # save plot
  ggsave(
    "plots/plot_stiffness.pdf",
    grd_E,
    width = 8, height = 3, units = "in"
  )
}
