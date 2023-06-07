# load packages
library("dplyr")
library("readr")
library("stringr")
library("magrittr")
library("ggplot2")

# load functions
source("functions/functions_statistics.R")
source("functions/functions_biomech.R")
source("functions/functions_plotting.R")


## LOAD AND ANALYSE LAT DATA --------------------------------------------------

# load LAT data
da <- read_csv("data/data_anatomical.csv") %>%
  rename_with(~str_trim(str_remove(.x, "\\[.*\\]"))) %>% # trim units from headers
  mutate(sample_id = as.character(sample_id))
# average area per sample
da_sum <- da %>%
  filter(test_succesful == TRUE) %>%
  group_by(sample_id) %>%
  summarize(
    root_area = mean(root_area, na.rm = TRUE),
    stele_area = mean(stele_area, na.rm = TRUE),
    cortex_area = mean(cortex_area, na.rm = TRUE)
  ) 
# fit stele/root area ratio as function of diameter, and plot
ft <- with(da_sum, linear_fit(sqrt(4/pi*root_area), stele_area/root_area, loglog = FALSE))
coef <- linear_fit_coefficients(ft, loglog = FALSE)
labels <- linear_fit_annotations(coef, label = c("intercept", "gradient", "R^2"))
linear_fit_ggplot(
  ft, 
  loglog = FALSE,
  label = labels,
  xlab = expression("Root diameter"~d[r]~"[mm]"),
  ylab = expression(A[s]*"/"*A[r]~"[mm"^2*"/mm"^2*"]"),
  xlim = c(0, 5),
  ylim = c(0, 1)
)
# save plot
ggsave(
  "plots/plot_stelerootarearatio.pdf", 
  width = 4, height = 3, units = "in"
)


## META-DATA ------------------------------------------------------------------

# load metadata for each tensile test
dm <- read_csv("data/data_tensiletest_metadata.csv") %>%
  rename_with(~str_trim(str_remove(.x, "\\[.*\\]"))) %>% # trim units from headers
  mutate(sample_id = as.character(sample_id))


## FIT TRACES - OBTAIN YIELD STRESS AND STRENGTH ------------------------------

# load all tensile data - only keep successful tests for which LAT data available
dt <- read_csv("data/data_tensiletests.csv") %>%
  rename_with(~str_trim(str_remove(.x, "\\[.*\\]"))) %>% # trim units from headers
  mutate(sample_id = as.character(sample_id)) %>%
  filter(
    sample_id %in% da$sample_id,
    sample_id %in% dm$sample_id[dm$test_succesful == TRUE]
  )
# fit all traces - find all key points along curve
df <- dt %>%
  group_by(sample_id) %>%
  reframe(fit_tensiletest_linear(displacement, force))
# merge root length and area to obtain stresses
df2 <- df %>%
  left_join(dm %>% select(sample_id, span_length, diameter), by = "sample_id") %>%
  left_join(da_sum %>% select(sample_id, root_area), by = "sample_id")
# get yield and ultimate stress and strain
df3 <- df2 %>% 
  group_by(sample_id, diameter) %>%
  summarize(
    u_rt = displacement[3] - displacement[1],
    tortuosity = u_rt/span_length[1],
    t_ru = force[9]/root_area[1],
    eps_ru = (displacement[9] - u_rt)/(span_length[1] + u_rt),
    t_ry = force[6]/root_area[1],
    eps_ry = (displacement[6] - u_rt)/(span_length[1] + u_rt)
  )


## PLOT ALL FORCE-DISPLACEMENT CURVES -----------------------------------------

# measured traces, normalised by span length and root area (diameter-based)
dtp <- dt %>% 
  left_join(dm %>% select(sample_id, diameter, span_length), by = "sample_id") %>%
  mutate(
    uL = displacement/span_length,
    t = force/(pi/4*diameter^2)
  )
# fitted traces
dfit <- df2 %>% 
  filter(point %in% c(1, 4, 6, 9)) %>%
  mutate(
    uL = displacement/span_length,
    t = force/(pi/4*diameter^2)
  )
# plot
ggplot() + 
  theme_bw() + 
  geom_path(
    data = dtp,
    aes(x = uL, y = t, color = "Measured data", linetype = "Measured data")
  ) + 
  geom_path(
    data = dfit,
    aes(x = uL, y = t, color = "Fit", linetype = "Fit")
  ) +
  geom_point(
    data = dfit,
    aes(x = uL, y = t, color = "Fit"),
    size = 0.75
  ) +
  facet_wrap(~sample_id) + 
  xlab(expression("Normalised displacement"~u[r]/L[s]~"[mm/mm]")) +
  ylab(expression("Tensile stress"~t[r]~"[MPa]")) +
  scale_color_brewer(name = NULL, palette = "Set1") +
  scale_linetype_manual(name = NULL, values = c(5, 1)) + 
  coord_cartesian(
    xlim = c(0, 0.25), 
    ylim = c(0, 20), 
    expand = FALSE
  ) + 
  scale_x_continuous(breaks = c(0, 0.10, 0.20)) +
  theme(
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines")
  )
# plot all mreasured traces + fitted linear traces
ggsave(
  "plots/plot_stressstrain_traces_fits.pdf", 
  width = 7, 
  height = 9, 
  units  = "in"
)

# example stress-strain fit
sample_select <- "578"
dt_select <- dt %>% filter(sample_id == sample_select)
df_select <- df2 %>% filter(sample_id == sample_select)
plot_tensiletest_linear(
  dt_select$displacement, dt_select$force, 
  df_select$displacement, df_select$force,
)
ggsave(
  "plots/plot_example_stressstrain_interpretation.pdf",
  width = 5, height = 3, units = "in"
)


## MERGE ALL DATA TOGETHER INTO A SINGLE DATASET, AND SAVE --------------------

# merge tensile points data with anatomical data
df_out <- df3 %>%
  select(sample_id, diameter, eps_ry, eps_ru, t_ry, t_ru, tortuosity) %>% 
  left_join(da_sum, by = "sample_id")
# save data to file
write_csv(
  df_out %>% rename(
    "diameter [mm]" = diameter,
    "A_r [mm^2]" = root_area,
    "A_s [mm^2]" = stele_area,
    "A_c [mm^2]" = cortex_area,
    "eps_ry [mm/mm]" = eps_ry,
    "t_ry [MPa]" = t_ry,
    "eps_ru [mm/mm]" = eps_ru,
    "t_ru [MPa]" = t_ru,
    "tortuosity [mm/mm]" = tortuosity
  ),
  "data_interpreted/data_rootproperties.csv"
)


