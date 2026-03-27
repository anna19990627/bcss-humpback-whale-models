
################################################################################
# Humpback whale occurrence model (sdmTMB) — BCSS megafauna surveys
# Bazaruto Archipelago (Mozambique), 2020–2024
#
# Workflow:
#   1) Load survey data, project to UTM, and standardize candidate predictors
#   2) Build the spatial mesh and coastline barrier mesh
#   3) Screen model components and candidate predictors with AIC
#   4) Build the final sequential-week time index with explicit gap weeks
#   5) Validate the final model and save diagnostic outputs
#   6) Plot marginal effects on unscaled axes
#   7) Produce present-day prediction maps
#   8) Project occurrence under SST scenarios and export rasters
#   9) Summarize projection uncertainty above probability thresholds
#  10) Check mesh-cutoff sensitivity
#
# Outputs:
#   figs/sdm_env/         diagnostics, effect plots, maps
#   tables/sdm_env/       AIC tables, validation, uncertainty, sensitivity
#   predictions/sdm_env/  GeoTIFFs
#   models/sdm_env/       fitted models and intermediate objects
################################################################################


# ==============================================================================
# 0) Packages and global settings
# ==============================================================================
library(tidyverse)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

library(sf)
library(sfheaders)
library(lubridate)

library(sdmTMB)
library(sdmTMBextra)
library(fmesher)

library(geodata)
library(terra)

library(ggplot2)
library(scales)
library(pROC)

library(mgcv)
library(visreg)
library(patchwork)

set.seed(123)


# ==============================================================================
# 1) Paths and output folders
# ==============================================================================
data_path <- "clean_dfs/obs_plot_short.rds"

fig_dir   <- file.path("figs", "sdm_env")
tab_dir   <- file.path("tables", "sdm_env")
pred_dir  <- file.path("predictions", "sdm_env")
model_dir <- file.path("models", "sdm_env")

dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(pred_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

WS_FILE <- file.path(model_dir, "SDM_HW.RData")
if (file.exists(WS_FILE)) load(WS_FILE)

# Scenario NetCDF paths
nc_path_ssp585_2090 <- "Data/original/climate_predictions/thetao_ssp585_2020_2100_depthsurf_f973_3101_d79e_U1748787933596.nc"
nc_path_ssp245_2090 <- "Data/original/climate_predictions/thetao_ssp245_2020_2100_depthsurf_f973_3101_d79e_U1748787928828.nc"
nc_path_ssp585_2040 <- "Data/original/climate_predictions/thetao_ssp585_2020_2100_depthsurf_6832_f862_ae1a_U1748787924215.nc"
nc_path_ssp245_2040 <- "Data/original/climate_predictions/thetao_ssp245_2020_2100_depthsurf_6832_f862_ae1a_U1748787912787.nc"
nc_path_baseline    <- "Data/original/climate_predictions/thetao_baseline_2000_2019_depthsurf_74ff_39fa_9adc_U1748789830446.nc"


# ==============================================================================
# 2) Data preparation: sf -> df, UTM (km), standardized predictors
# ==============================================================================
obs_test <- readRDS(data_path)

# Harmonize known column-name typos
obs_test <- obs_test %>%
  rename_with(~ "current_speed",   .cols = any_of("currend_speed")) %>%
  rename_with(~ "max_wave_height", .cols = any_of("max_wave_hight"))

# Convert sf to a regular data frame and add UTM coordinates in km
obs_sdm <- obs_test |>
  sfheaders::sf_to_df(fill = TRUE) |>
  select(-any_of(c("sfg_id", "point_id"))) |>
  rename(longitude = x, latitude = y) |>
  add_utm_columns(
    ll_names  = c("longitude", "latitude"),
    ll_crs    = 4326,
    utm_names = c("x", "y"),
    utm_crs   = 32736,
    units     = "km"
  )

# Standardize candidate predictors but keep raw columns
obs_sdm <- obs_sdm %>%
  mutate(
    across(
      .cols = c(
        SST,
        depth,
        wind_speed,
        swell_period,
        max_wave_height,
        current_speed,
        tidal_amplitude,
        dist_m
      ),
      .fns   = ~ as.numeric(scale(.x)),
      .names = "{.col}_scaled"
    ),
    year    = as.factor(year),
    doy     = as.integer(lubridate::yday(Date_Time)),
    ToD     = round(lubridate::hour(Date_Time) + lubridate::minute(Date_Time) / 60, digits = 1),
    yearint = as.integer(year)
  )

saveRDS(obs_sdm, file.path(model_dir, "obs_sdm_prepped.rds"))


# ==============================================================================
# 3) Mesh construction: base mesh and barrier mesh
# ==============================================================================
# Base mesh (coordinates are in km)
mesh1 <- make_mesh(
  obs_sdm,
  c("x", "y"),
  cutoff = 2
)

# Mozambique boundary (level 1), cropped to the broader Bazaruto region
mz <- geodata::gadm("Mozambique", level = 1, path = ".") |> sf::st_as_sf()

bb <- sf::st_bbox(
  c(xmin = 35, xmax = 36.0, ymin = -22.5, ymax = -21.0),
  crs = sf::st_crs(4326)
)
mz_sf <- sf::st_crop(mz, bb)

# Transform to UTM metres, then to UTM km for the barrier mesh
mz_utm_m <- sf::st_transform(mz_sf, crs = 32736)

mz_utm_km <- mz_utm_m |>
  sf::st_transform(
    crs = fm_crs_set_lengthunit(
      crs = fm_CRS(SRS_string = "EPSG:32736"), "km"
    )
  )

# Barrier mesh: proj_scaling = 1000 because x/y are in km
mesh2 <- sdmTMBextra::add_barrier_mesh(
  mesh1,
  mz_utm_km,
  range_fraction = 0.1,
  proj_scaling   = 1000,
  plot           = FALSE
)


# ==============================================================================
# 4) Predictor screening: collinearity check
# ==============================================================================
num_vars <- obs_sdm %>%
  select(
    SST_scaled,
    depth_scaled,
    dist_m_scaled,
    wind_speed_scaled,
    swell_period_scaled,
    max_wave_height_scaled,
    current_speed_scaled,
    tidal_amplitude_scaled,
    ToD
  ) %>%
  drop_na()

cor_mat <- cor(num_vars, use = "pairwise.complete.obs", method = "pearson")
cor_tbl <- as.data.frame(cor_mat) %>%
  tibble::rownames_to_column("variable")
write_csv(cor_tbl, file.path(tab_dir, "predictor_correlations.csv"))


# ==============================================================================
# 5) Exploratory model checks (AIC screening; weeknr index)
# ==============================================================================
# 5.1 AR1 vs RW spatiotemporal structure
m_int_ar1 <- sdmTMB(
  observed ~ 1,
  data           = obs_sdm,
  mesh           = mesh2,
  time           = "weeknr",
  family         = binomial(link = "logit"),
  spatial        = "on",
  spatiotemporal = "ar1",
  reml           = TRUE
)

m_int_rw <- sdmTMB(
  observed ~ 1,
  data           = obs_sdm,
  mesh           = mesh2,
  time           = "weeknr",
  family         = binomial(link = "logit"),
  spatial        = "on",
  spatiotemporal = "rw",
  reml           = TRUE
)

aic_ar1_rw <- AIC(m_int_ar1, m_int_rw) |>
  as.data.frame() %>%
  tibble::rownames_to_column("model") %>%
  arrange(AIC)
write_csv(aic_ar1_rw, file.path(tab_dir, "AIC_intercept_AR1_vs_RW.csv"))

rho_tbl <- tidy(m_int_ar1, effects = "ran_pars") |> filter(term == "rho")
write_csv(rho_tbl, file.path(tab_dir, "rho_intercept_AR1.csv"))

# 5.2 Year random intercept vs none
m_int <- sdmTMB(
  observed ~ 1,
  data           = obs_sdm,
  mesh           = mesh2,
  time           = "weeknr",
  family         = binomial(link = "logit"),
  spatial        = "on",
  spatiotemporal = "ar1",
  reml           = TRUE
)

m_intre <- sdmTMB(
  observed ~ 1 + (1 | year),
  data           = obs_sdm,
  mesh           = mesh2,
  time           = "weeknr",
  family         = binomial(link = "logit"),
  spatial        = "on",
  spatiotemporal = "ar1",
  reml           = TRUE
)

aic_year_re <- AIC(m_int, m_intre) |>
  as.data.frame() %>%
  tibble::rownames_to_column("model") %>%
  arrange(AIC)
write_csv(aic_year_re, file.path(tab_dir, "AIC_intercept_yearRE.csv"))

# 5.3 Linear vs smooth distance-to-shore effect
m_distlin <- sdmTMB(
  observed ~ dist_m_scaled + (1 | year),
  data           = obs_sdm,
  mesh           = mesh2,
  time           = "weeknr",
  family         = binomial(link = "logit"),
  spatial        = "on",
  spatiotemporal = "ar1",
  reml           = TRUE
)

m_distsm <- sdmTMB(
  observed ~ s(dist_m_scaled, bs = "ts") + (1 | year),
  data           = obs_sdm,
  mesh           = mesh2,
  time           = "weeknr",
  family         = binomial(link = "logit"),
  spatial        = "on",
  spatiotemporal = "ar1",
  reml           = TRUE
)

m_distsm2 <- sdmTMB(
  observed ~ s(dist_m_scaled, bs = "cs") + (1 | year),
  data           = obs_sdm,
  mesh           = mesh2,
  time           = "weeknr",
  family         = binomial(link = "logit"),
  spatial        = "on",
  spatiotemporal = "ar1",
  reml           = TRUE
)

aic_distance <- AIC(m_distlin, m_distsm, m_distsm2) |>
  as.data.frame() %>%
  tibble::rownames_to_column("model") %>%
  arrange(AIC)
write_csv(aic_distance, file.path(tab_dir, "AIC_distance_linear_vs_smooth.csv"))

# 5.4 Need for spatial and spatiotemporal random fields
m_intre_nospatem <- sdmTMB(
  observed ~ 1 + (1 | year),
  data           = obs_sdm,
  mesh           = mesh2,
  family         = binomial(link = "logit"),
  spatial        = "on",
  spatiotemporal = "off",
  reml           = TRUE
)

m_intre_nospa <- sdmTMB(
  observed ~ 1 + (1 | year),
  data           = obs_sdm,
  mesh           = mesh2,
  family         = binomial(link = "logit"),
  spatial        = "off",
  spatiotemporal = "off",
  reml           = TRUE
)

aic_spatial <- AIC(m_intre, m_intre_nospatem, m_intre_nospa) |>
  as.data.frame() %>%
  tibble::rownames_to_column("model") %>%
  arrange(AIC)
write_csv(aic_spatial, file.path(tab_dir, "AIC_spatial_components.csv"))

# 5.5 Single-predictor and candidate models
m_sst <- sdmTMB(
  observed ~ s(SST_scaled, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_tod <- sdmTMB(
  observed ~ s(ToD, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_tod_cyclic <- sdmTMB(
  observed ~ s(ToD, bs = "cc", k = 8) + (1 | year),
  data           = obs_sdm,
  mesh           = mesh2,
  time           = "weeknr",
  family         = binomial(link = "logit"),
  spatial        = "on",
  spatiotemporal = "ar1",
  knots          = list(ToD = c(0, 24)),
  reml           = TRUE
)

m_depth <- sdmTMB(
  observed ~ s(depth_scaled, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_windspeed <- sdmTMB(
  observed ~ s(wind_speed_scaled, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_swellperiod <- sdmTMB(
  observed ~ s(swell_period_scaled, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_maxwave <- sdmTMB(
  observed ~ s(max_wave_height_scaled, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_currentspeed <- sdmTMB(
  observed ~ s(current_speed_scaled, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_tidalampl <- sdmTMB(
  observed ~ s(tidal_amplitude_scaled, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_dist_tidalampl <- sdmTMB(
  observed ~ s(dist_m_scaled, bs = "ts") + s(tidal_amplitude_scaled, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_distXtidalampl <- sdmTMB(
  observed ~ t2(dist_m_scaled, tidal_amplitude_scaled, bs = c("ts", "ts")) + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_depthXsst <- sdmTMB(
  observed ~ t2(depth_scaled, SST_scaled, bs = c("ts", "ts")) + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_depth_sst <- sdmTMB(
  observed ~ s(depth_scaled, bs = "ts") + s(SST_scaled, bs = "ts") + (1 | year),
  data           = obs_sdm,
  mesh           = mesh2,
  time           = "weeknr",
  family         = binomial(link = "logit"),
  spatial        = "on",
  spatiotemporal = "ar1",
  reml           = TRUE
)

m_sst_dist <- sdmTMB(
  observed ~ s(SST_scaled, bs = "ts") + s(dist_m_scaled, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_sstXdist <- sdmTMB(
  observed ~ t2(SST_scaled, dist_m_scaled, bs = c("ts", "ts")) + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_depth_dist <- sdmTMB(
  observed ~ s(depth_scaled, bs = "ts") + s(dist_m_scaled, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_depthXdist <- sdmTMB(
  observed ~ t2(depth_scaled, dist_m_scaled, bs = c("ts", "ts")) + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_depth_sst_dist <- sdmTMB(
  observed ~ s(depth_scaled, bs = "ts") + s(SST_scaled, bs = "ts") + s(dist_m_scaled, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_depthXsst_dist <- sdmTMB(
  observed ~ t2(depth_scaled, SST_scaled, bs = c("ts", "ts")) + s(dist_m_scaled, bs = "ts") + (1 | year),
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

# Time-varying effects
m_TVCdepth2 <- sdmTMB(
  observed ~ 1 + (1 | year),
  time_varying = ~depth_scaled,
  time_varying_type = "rw0",
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_TVCdepth3 <- sdmTMB(
  observed ~ 1 + (1 | year),
  time_varying = ~depth_scaled,
  time_varying_type = "ar1",
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_depth_TVCdepth <- sdmTMB(
  observed ~ 1 + s(depth_scaled, bs = "ts") + (1 | year),
  time_varying = ~depth_scaled,
  time_varying_type = "rw0",
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_TVCsst <- sdmTMB(
  observed ~ 1 + (1 | year),
  time_varying = ~SST_scaled,
  time_varying_type = "rw0",
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

m_sst_TVCsst <- sdmTMB(
  observed ~ 1 + s(SST_scaled, bs = "ts") + (1 | year),
  time_varying = ~SST_scaled,
  time_varying_type = "rw0",
  data = obs_sdm,
  mesh = mesh2,
  time = "weeknr",
  family = binomial(link = "logit"),
  spatial = "on",
  spatiotemporal = "ar1",
  reml = TRUE
)

aic_weeknr <- AIC(
  m_intre,
  m_distsm,
  m_sst,
  m_tod_cyclic,
  m_depth,
  m_windspeed,
  m_swellperiod,
  m_maxwave,
  m_currentspeed,
  m_tidalampl,
  m_dist_tidalampl,
  m_distXtidalampl,
  m_depth_sst,
  m_depthXsst,
  m_sst_dist,
  m_sstXdist,
  m_depth_dist,
  m_depthXdist,
  m_depth_sst_dist,
  m_depthXsst_dist,
  m_TVCdepth2,
  m_depth_TVCdepth,
  m_TVCsst,
  m_sst_TVCsst
) %>%
  as.data.frame() %>%
  rownames_to_column("model") %>%
  arrange(AIC)

write_csv(aic_weeknr, file.path(tab_dir, "AIC_weeknr_candidate_set.csv"))
print(aic_weeknr)

sdm_model_labels <- tibble::tribble(
  ~model,               ~Model,
  "m_depth_sst_dist",   "Depth + SST + Distance",
  "m_depth_sst",        "Depth + SST",
  "m_sstXdist",         "SST × Distance",
  "m_depthXsst_dist",   "Depth × SST + Distance",
  "m_depthXsst",        "Depth × SST",
  "m_sst",              "SST only",
  "m_sst_dist",         "SST + Distance",
  "m_sst_TVCsst",       "SST + time-varying SST",
  "m_TVCsst",           "Time-varying SST",
  "m_depth_TVCdepth",   "Depth + time-varying Depth",
  "m_depth",            "Depth only",
  "m_depth_dist",       "Depth + Distance",
  "m_depthXdist",       "Depth × Distance",
  "m_distXtidalampl",   "Distance × Tide amplitude",
  "m_dist_tidalampl",   "Distance + Tide amplitude",
  "m_TVCdepth2",        "Time-varying Depth (rw0)",
  "m_intre",            "Intercept + year RE",
  "m_tod_cyclic",       "Time of day (cyclic)",
  "m_distsm",           "Distance (smooth)",
  "m_currentspeed",     "Current speed",
  "m_swellperiod",      "Swell period",
  "m_maxwave",          "Max wave height",
  "m_windspeed",        "Wind speed",
  "m_tidalampl",        "Tide amplitude"
)

tab_sdm_weeknr <- aic_weeknr %>%
  mutate(
    dAIC = AIC - min(AIC),
    wAIC = exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC))
  ) %>%
  left_join(sdm_model_labels, by = "model") %>%
  mutate(
    Model = dplyr::coalesce(Model, model)
  ) %>%
  arrange(AIC) %>%
  transmute(
    Model,
    AIC  = round(AIC, 2),
    dAIC = round(dAIC, 2),
    wAIC = round(wAIC, 2)
  )

knitr::kable(
  tab_sdm_weeknr,
  booktabs = TRUE,
  caption = paste0(
    "Table Sx. Habitat model selection for humpback whale occurrence using AIC (ranked by increasing AIC). ",
    "All candidate models were fitted with the same spatiotemporal structure (binomial logit; spatial Matérn field on a 2-km barrier mesh; ",
    "AR(1) spatiotemporal random field on weeknr; year random intercept). ",
    "dAIC = AIC difference to the best model; wAIC = Akaike weight."
  )
)

readr::write_csv(tab_sdm_weeknr, file.path(tab_dir, "TableSx_sdm_model_selection_weeknr.csv"))


# ==============================================================================
# 6) Final time index: sequential week with explicit gaps
# ==============================================================================
origin <- as.Date(min(obs_sdm$Date_Time, na.rm = TRUE))

obs_sdm2 <- obs_sdm %>%
  mutate(
    week_global = as.integer(floor(as.numeric(as.Date(Date_Time) - origin) / 7) + 1L)
  ) %>%
  filter(
    !is.na(week_global),
    !is.na(x), !is.na(y),
    !is.na(observed),
    !is.na(depth_scaled),
    !is.na(SST_scaled),
    !is.na(year)
  )

all_weeks <- seq.int(min(obs_sdm2$week_global), max(obs_sdm2$week_global))

stopifnot(is.atomic(all_weeks))
stopifnot(is.atomic(obs_sdm2$week_global))

m_depth_sst_ar1gap <- sdmTMB(
  observed ~ s(depth_scaled, bs = "ts") + s(SST_scaled, bs = "ts") + (1 | year),
  data           = obs_sdm2,
  mesh           = mesh2,
  time           = "week_global",
  extra_time     = all_weeks,
  family         = binomial(link = "logit"),
  spatial        = "on",
  spatiotemporal = "ar1",
  reml           = TRUE
)

saveRDS(m_depth_sst_ar1gap, file.path(model_dir, "m_depth_sst_ar1gap.rds"))

obs_sdm_timeid <- obs_sdm %>%
  mutate(
    season = lubridate::year(Date_Time),
    week_in_season = lubridate::week(Date_Time),
    time_id = interaction(season, week_in_season, drop = TRUE)
  )

m_depth_sst_iid <- sdmTMB(
  observed ~ s(depth_scaled, bs = "ts") + s(SST_scaled, bs = "ts") + (1 | year),
  data           = obs_sdm_timeid,
  mesh           = mesh2,
  time           = "time_id",
  family         = binomial(link = "logit"),
  spatial        = "on",
  spatiotemporal = "iid",
  reml           = TRUE
)

aic_time_compare <- AIC(m_depth_sst, m_depth_sst_iid, m_depth_sst_ar1gap) |>
  as.data.frame() %>%
  tibble::rownames_to_column("model") %>%
  arrange(AIC)
write_csv(aic_time_compare, file.path(tab_dir, "AIC_final_time_structures.csv"))

aic_time <- tibble::tibble(
  Model = c(
    "Depth + SST (weeknr, AR1)",
    "Depth + SST (time_id, IID)",
    "Depth + SST (week_global + extra_time, AR1) (final)"
  ),
  AIC = c(
    as.numeric(AIC(m_depth_sst)),
    as.numeric(AIC(m_depth_sst_iid)),
    as.numeric(AIC(m_depth_sst_ar1gap))
  )
) %>%
  mutate(
    dAIC = AIC - min(AIC),
    wAIC = exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC))
  ) %>%
  arrange(AIC) %>%
  mutate(across(c(AIC, dAIC, wAIC), ~ round(.x, 2)))

knitr::kable(
  aic_time,
  booktabs = TRUE,
  caption = paste0(
    "Table Sy. Time-structure comparison for the final occurrence model using AIC. ",
    "All models use the same fixed effects (s(Depth) + s(SST) + year random intercept) and the same spatial barrier mesh; ",
    "models differ only in the time index and spatiotemporal structure. ",
    "The sequential-week AR(1) model with extra_time was retained as the final structure."
  )
)

readr::write_csv(aic_time, file.path(tab_dir, "TableSy_sdm_time_structure_selection.csv"))


# ==============================================================================
# 7) Model validation and diagnostics
# ==============================================================================
FIT <- m_depth_sst_ar1gap

# Dunn & Smyth randomized quantile residuals
obs_sdm2$resids <- residuals(FIT)

png(file.path(fig_dir, "diag_QQ_residuals.png"), width = 1400, height = 900, res = 180)
qqnorm(obs_sdm2$resids)
qqline(obs_sdm2$resids)
dev.off()

png(file.path(fig_dir, "diag_residuals_map_by_year.png"), width = 1600, height = 900, res = 180)
ggplot(obs_sdm2, aes(x, y, col = resids)) +
  geom_point(size = 1) +
  facet_wrap(~ year) +
  coord_fixed() +
  scale_colour_gradient2() +
  theme_minimal()
dev.off()

# Fixed-effects-only probabilities for validation metrics
p_fixed <- plogis(predict(FIT, re_form = NA, re_form_iid = NA)$est)
roc_c   <- pROC::roc(obs_sdm2$observed, p_fixed, quiet = TRUE)

png(file.path(fig_dir, "diag_ROC_fixed_effects.png"), width = 1400, height = 900, res = 180)
plot(roc_c, legacy.axes = TRUE, main = "ROC curve (fixed effects only)")
abline(a = 0, b = 1, lty = 2, col = "grey50")
dev.off()

auc_val <- as.numeric(pROC::auc(roc_c))
auc_ci  <- as.numeric(pROC::ci.auc(roc_c))

tjur_D <- mean(p_fixed[obs_sdm2$observed == 1], na.rm = TRUE) -
  mean(p_fixed[obs_sdm2$observed == 0], na.rm = TRUE)

eps <- 1e-8
p_clip <- pmin(pmax(p_fixed, eps), 1 - eps)
brier   <- mean((obs_sdm2$observed - p_clip)^2, na.rm = TRUE)
logloss <- -mean(obs_sdm2$observed * log(p_clip) + (1 - obs_sdm2$observed) * log(1 - p_clip), na.rm = TRUE)

val_tbl <- tibble(
  AUC = auc_val,
  AUC_CI_low = auc_ci[1],
  AUC_CI_high = auc_ci[3],
  Tjur_D = tjur_D,
  Brier = brier,
  LogLoss = logloss
)
write_csv(val_tbl, file.path(tab_dir, "sdm_validation_metrics.csv"))

# Calibration intercept and slope
logit_p <- qlogis(p_clip)
cal_int   <- glm(observed ~ 1, family = binomial(), offset = logit_p, data = obs_sdm2)
cal_slope <- glm(observed ~ logit_p, family = binomial(), data = obs_sdm2)
cal_tbl <- tibble(
  cal_intercept = coef(cal_int)[1],
  cal_slope_intercept = coef(cal_slope)[1],
  cal_slope = coef(cal_slope)[2]
)
write_csv(cal_tbl, file.path(tab_dir, "sdm_calibration.csv"))

save.image(file = WS_FILE)

N_total <- nrow(obs_sdm2)
N_pres  <- sum(obs_sdm2$observed == 1, na.rm = TRUE)
N_abs   <- sum(obs_sdm2$observed == 0, na.rm = TRUE)
c(N_total = N_total, N_pres = N_pres, N_abs = N_abs)

# N_total  N_pres   N_abs
#  2876     831    2045


# ==============================================================================
# 8) Effect plots on unscaled axes
# ==============================================================================
FIT <- m_depth_sst_ar1gap

# Scaling constants from the fitted data
mu_depth <- mean(obs_sdm2$depth, na.rm = TRUE)
sd_depth <- sd(obs_sdm2$depth,   na.rm = TRUE)

mu_SST <- mean(obs_sdm2$SST, na.rm = TRUE)
sd_SST <- sd(obs_sdm2$SST,   na.rm = TRUE)

ref_year <- factor(levels(obs_sdm2$year)[1], levels = levels(obs_sdm2$year))
ref_week <- as.integer(round(median(obs_sdm2$week_global, na.rm = TRUE)))

# Some sdmTMB versions require x/y even when random fields are excluded
ref_x <- mean(obs_sdm2$x, na.rm = TRUE)
ref_y <- mean(obs_sdm2$y, na.rm = TRUE)

# Convert simulated predictions to an n x nsim matrix
extract_sims_matrix <- function(pred_obj, nsim) {
  if (is.matrix(pred_obj)) return(pred_obj)

  sims_attr <- attr(pred_obj, "sims")
  if (!is.null(sims_attr) && is.matrix(sims_attr)) return(sims_attr)

  if (is.data.frame(pred_obj)) {
    sim_cols <- grep("^sim", names(pred_obj), value = TRUE)
    if (length(sim_cols) == nsim) return(as.matrix(pred_obj[, sim_cols, drop = FALSE]))
  }

  if (is.list(pred_obj) && !is.null(pred_obj$sims) && is.matrix(pred_obj$sims)) {
    return(pred_obj$sims)
  }

  stop("Could not extract a simulation matrix from predict(..., nsim = ...).")
}

# Fixed-effects-only curve predictions with simulation-based uncertainty
predict_curve_prob <- function(var = c("depth", "SST"),
                               n = 200,
                               trim = c(0.05, 0.95),
                               nsim = 500,
                               seed = 123) {
  var <- match.arg(var)

  x_raw <- obs_sdm2[[var]]
  x_seq <- seq(
    quantile(x_raw, trim[1], na.rm = TRUE),
    quantile(x_raw, trim[2], na.rm = TRUE),
    length.out = n
  )

  depth_raw <- if (var == "depth") x_seq else rep(mu_depth, n)
  sst_raw   <- if (var == "SST")   x_seq else rep(mu_SST,   n)

  nd <- tibble(
    x = ref_x, y = ref_y,
    year = ref_year,
    week_global = ref_week,
    depth = depth_raw,
    SST   = sst_raw,
    depth_scaled = (depth_raw - mu_depth) / sd_depth,
    SST_scaled   = (sst_raw   - mu_SST)   / sd_SST
  )

  set.seed(seed)
  sim_eta <- predict(
    FIT,
    newdata = nd,
    type = "link",
    nsim = nsim,
    re_form = NA,
    re_form_iid = NA
  )

  sim_mat <- extract_sims_matrix(sim_eta, nsim)
  sim_p   <- plogis(sim_mat)

  tibble(
    x = x_seq,
    prob_med = apply(sim_p, 1, median, na.rm = TRUE),
    prob_lo  = apply(sim_p, 1, quantile, probs = 0.025, na.rm = TRUE),
    prob_hi  = apply(sim_p, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}

df_depth <- predict_curve_prob("depth", nsim = 500)
df_sst   <- predict_curve_prob("SST",   nsim = 500)

pmax <- max(df_sst$prob_med, na.rm = TRUE)

peak_band <- df_sst %>%
  filter(prob_med >= 0.95 * pmax) %>%
  summarise(
    sst_low  = min(x, na.rm = TRUE),
    sst_high = max(x, na.rm = TRUE)
  )

peak_band
#   sst_low sst_high
#   23.3     24.3

pmax_d <- max(df_depth$prob_med, na.rm = TRUE)

peak_band_depth <- df_depth %>%
  filter(prob_med >= 0.95 * pmax_d) %>%
  summarise(
    depth_low  = min(x, na.rm = TRUE),
    depth_high = max(x, na.rm = TRUE)
  )

peak_band_depth
#  depth_low depth_high
#      92.6       122.

df_sst2 <- df_sst %>% arrange(x)
pmax <- max(df_sst2$prob_med, na.rm = TRUE)

i_peak <- which.max(df_sst2$prob_med)
i_half <- which(df_sst2$prob_med <= 0.5 * pmax & seq_len(nrow(df_sst2)) >= i_peak)[1]
sst_half <- df_sst2$x[i_half]

c(peak_sst = df_sst2$x[i_peak], sst_half = sst_half)
# peak_sst sst_half
# 23.61156 27.13216

df_depth2 <- df_depth %>% arrange(x)
pmax_d <- max(df_depth2$prob_med, na.rm = TRUE)

i_peak_d <- which.max(df_depth2$prob_med)

i_half_d <- which(
  df_depth2$prob_med <= 0.5 * pmax_d &
    seq_len(nrow(df_depth2)) >= i_peak_d
)[1]

depth_half_drop <- df_depth2$x[i_half_d]

c(peak_depth = df_depth2$x[i_peak_d], depth_half_drop = depth_half_drop)
#     peak_depth depth_half_drop
#    122.25              NA

p_depth <- ggplot(df_depth, aes(x = x, y = prob_med)) +
  geom_ribbon(aes(ymin = prob_lo, ymax = prob_hi), alpha = 0.2) +
  geom_line(linewidth = 1.1) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Depth (m)",
    y = "Predicted probability",
    title = "HW occurrence along a depth gradient",
    subtitle = "Median and 95% intervals (fixed effects only)"
  ) +
  theme_classic(base_size = 13)

p_sst <- ggplot(df_sst, aes(x = x, y = prob_med)) +
  geom_ribbon(aes(ymin = prob_lo, ymax = prob_hi), alpha = 0.2) +
  geom_line(linewidth = 1.1) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Sea surface temperature (°C)",
    y = "Predicted probability",
    title = "HW occurrence along an SST gradient",
    subtitle = "Median and 95% intervals (fixed effects only)"
  ) +
  theme_classic(base_size = 13)

p_depth
p_sst

ggsave(file.path(fig_dir, "FigS_effect_depth_unscaled.png"), p_depth, width = 7.5, height = 4.8, dpi = 300)
ggsave(file.path(fig_dir, "FigS_effect_sst_unscaled.png"),   p_sst,   width = 7.5, height = 4.8, dpi = 300)

make_surface <- function(n_depth = 120, n_sst = 120, trim = c(0.05, 0.95)) {
  depth_seq <- seq(
    quantile(obs_sdm2$depth, trim[1], na.rm = TRUE),
    quantile(obs_sdm2$depth, trim[2], na.rm = TRUE),
    length.out = n_depth
  )

  sst_seq <- seq(
    quantile(obs_sdm2$SST, trim[1], na.rm = TRUE),
    quantile(obs_sdm2$SST, trim[2], na.rm = TRUE),
    length.out = n_sst
  )

  grid <- tidyr::expand_grid(depth = depth_seq, SST = sst_seq) %>%
    mutate(
      x = ref_x, y = ref_y,
      year = ref_year,
      week_global = ref_week,
      depth_scaled = (depth - mu_depth) / sd_depth,
      SST_scaled   = (SST   - mu_SST)   / sd_SST
    )

  pr <- predict(
    FIT,
    newdata = grid,
    type = "link",
    re_form = NA,
    re_form_iid = NA
  )

  grid$prob <- plogis(pr$est)
  grid
}

surf <- make_surface()

p_2d <- ggplot(surf, aes(x = depth, y = SST, fill = prob)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), oob = scales::squish) +
  labs(
    x = "Depth (m)",
    y = "SST (°C)",
    fill = "Prob.",
    title = "HW occurrence surface (Depth × SST)",
    subtitle = "Fixed effects only"
  ) +
  theme_classic(base_size = 13)

p_2d
ggsave(file.path(fig_dir, "FigS_effect_depth_x_sst_surface.png"), p_2d, width = 7.5, height = 5.5, dpi = 300)


# ==============================================================================
# 9) Present-day prediction map
# ==============================================================================
nx <- 200
ny <- 200

pred_grid0 <- expand.grid(
  x = seq(min(obs_sdm2$x, na.rm = TRUE), max(obs_sdm2$x, na.rm = TRUE), length.out = nx),
  y = seq(min(obs_sdm2$y, na.rm = TRUE), max(obs_sdm2$y, na.rm = TRUE), length.out = ny)
)

pred_pts <- st_as_sf(pred_grid0, coords = c("x", "y"), crs = st_crs(mz_utm_km))
on_land  <- lengths(st_intersects(pred_pts, st_union(mz_utm_km))) > 0
pred_grid <- pred_grid0[!on_land, , drop = FALSE] %>%
  mutate(
    cell_id      = row_number(),
    depth_scaled = mean(obs_sdm2$depth_scaled, na.rm = TRUE),
    SST_scaled   = mean(obs_sdm2$SST_scaled,   na.rm = TRUE)
  )

times_present <- obs_sdm2 %>%
  distinct(year, week_global) %>%
  arrange(year, week_global) %>%
  mutate(
    time_id = row_number(),
    year = factor(year, levels = levels(obs_sdm2$year)),
    week_global = as.integer(week_global)
  )

chunk_size <- 10L
idx_all <- seq_len(nrow(times_present))
chunks <- split(idx_all, ceiling(idx_all / chunk_size))

sum_prob <- numeric(nrow(pred_grid))
n_done   <- 0L

for (idx in chunks) {
  times_chunk <- times_present[idx, , drop = FALSE] %>% arrange(time_id)
  nd_chunk <- tidyr::crossing(pred_grid, times_chunk) %>% arrange(cell_id, time_id)

  pr <- predict(FIT, newdata = nd_chunk, type = "response")
  prob_mat <- matrix(pr$est, nrow = nrow(pred_grid), ncol = length(idx), byrow = TRUE)

  sum_prob <- sum_prob + rowSums(prob_mat, na.rm = TRUE)
  n_done   <- n_done + length(idx)

  rm(nd_chunk, pr, prob_mat)
  gc()
}

pred_mean_present <- pred_grid %>%
  mutate(prob_mean = sum_prob / n_done) %>%
  select(x, y, prob_mean)

p_present <- ggplot() +
  geom_raster(data = pred_mean_present, aes(x = x, y = y, fill = prob_mean)) +
  geom_sf(data = mz_utm_km, fill = "grey80", color = "black", inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "magma", name = "HW presence\nprobability\n(mean)") +
  coord_sf(
    crs = st_crs(mz_utm_km),
    default_crs = st_crs(mz_utm_km),
    expand = FALSE,
    datum = NA
  ) +
  labs(
    x = "Easting (km)",
    y = "Northing (km)",
    title = "Present-day mean predicted humpback whale occurrence\n(full conditional; averaged over observed weeks)"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(fig_dir, "Fig_present_day_mean_occurrence.png"), p_present,
       width = 7.5, height = 5.5, dpi = 300)


# ==============================================================================
# 10) Scenario projections (SST-driven change; conditional on the spatial field)
# ==============================================================================
park_coords <- tibble::tribble(
  ~Longitude, ~Latitude,
  35.54167, -21.45833,
  35.31694, -21.45833,
  35.31694, -21.71306,
  35.36972, -21.92611,
  35.40028, -22.04861,
  35.54167, -22.04861,
  35.54167, -21.92611,
  35.54167, -21.45833
)

sanctuary_coords <- tibble::tribble(
  ~Longitude, ~Latitude,
  35.40028, -22.04861,
  35.54167, -22.04861,
  35.55861, -22.17461,
  35.55933, -22.34078,
  35.54389, -22.34153,
  35.42411, -22.33653,
  35.41489, -22.19300,
  35.40028, -22.04861
)

make_poly_utm <- function(df_ll, crs_in = 4326, crs_out = 32736) {
  m <- as.matrix(df_ll[, c("Longitude", "Latitude")])
  if (!all(m[1, ] == m[nrow(m), ])) m <- rbind(m, m[1, ])
  sf::st_sfc(sf::st_polygon(list(m)), crs = crs_in) %>% sf::st_transform(crs_out)
}

park_utm      <- make_poly_utm(park_coords)
sanctuary_utm <- make_poly_utm(sanctuary_coords)

mu_SST   <- mean(obs_sdm2$SST, na.rm = TRUE)
sd_SST   <- sd(obs_sdm2$SST, na.rm = TRUE)
mu_depth <- mean(obs_sdm2$depth, na.rm = TRUE)
sd_depth <- sd(obs_sdm2$depth, na.rm = TRUE)
depth_mean_scaled <- mean(obs_sdm2$depth_scaled, na.rm = TRUE)

grid_res_km <- 0.5
margin_km   <- 2

bbox_km <- list(
  xmin = min(obs_sdm2$x, na.rm = TRUE) - margin_km,
  xmax = max(obs_sdm2$x, na.rm = TRUE) + margin_km,
  ymin = min(obs_sdm2$y, na.rm = TRUE) - margin_km,
  ymax = max(obs_sdm2$y, na.rm = TRUE) + margin_km
)

r_template <- terra::rast(
  xmin = bbox_km$xmin * 1000, xmax = bbox_km$xmax * 1000,
  ymin = bbox_km$ymin * 1000, ymax = bbox_km$ymax * 1000,
  resolution = grid_res_km * 1000,
  crs = "EPSG:32736"
)

# Land mask in metres
land_sf_m <- mz_utm_m %>%
  st_make_valid() %>%
  st_union() %>%
  st_buffer(100) %>%
  st_simplify(dTolerance = 50)

land_r <- terra::rasterize(terra::vect(land_sf_m), r_template, field = 1)

# Load, average, project, and resample SST rasters
load_sst_mean_on_grid <- function(nc_path, start_year, end_year, template) {
  stopifnot(file.exists(nc_path))
  s <- terra::rast(nc_path)
  try(terra::crs(s) <- "EPSG:4326", silent = TRUE)

  tt <- try(terra::time(s), silent = TRUE)
  if (!inherits(tt, "try-error") && !all(is.na(tt))) {
    idx <- which(lubridate::year(tt) >= start_year & lubridate::year(tt) <= end_year)
    if (length(idx)) s <- s[[idx]]
  }

  s_mean <- terra::app(s, mean, na.rm = TRUE)

  mn <- as.numeric(terra::global(s_mean, "min", na.rm = TRUE)[1, 1])
  if (!is.na(mn) && mn > 100) s_mean <- s_mean - 273.15

  s_utm <- terra::project(s_mean, terra::crs(template), method = "bilinear")
  terra::resample(s_utm, template, method = "bilinear")
}

# Fill coastal NAs introduced by reprojection/resampling
fill_na_nearestish <- function(r, max_w = 21) {
  out <- r
  for (w in seq(3, max_w, by = 2)) {
    if (terra::global(is.na(out), "sum", na.rm = TRUE)[1, 1] == 0) break
    sm <- terra::focal(out, w = matrix(1, w, w), fun = mean, na.rm = TRUE, expand = TRUE)
    out <- terra::cover(out, sm)
  }
  out
}

# Scenario maps use fixed effects plus the fitted spatial field only
predict_from_sst_r <- function(sst_r, label) {
  xy_m <- as.data.frame(terra::xyFromCell(r_template, 1:terra::ncell(r_template)))
  vals <- as.numeric(terra::values(sst_r))

  base_df <- tibble(
    cell = 1:terra::ncell(r_template),
    x    = xy_m$x / 1000,
    y    = xy_m$y / 1000,
    SST  = vals
  ) %>%
    filter(is.finite(SST))

  base_df$depth_scaled <- depth_mean_scaled
  base_df$SST_scaled   <- (base_df$SST - mu_SST) / sd_SST
  base_df$week_global  <- as.integer(round(stats::median(obs_sdm2$week_global, na.rm = TRUE)))
  base_df$year         <- factor(levels(obs_sdm2$year)[1], levels = levels(obs_sdm2$year))

  needed <- intersect(c("SST_scaled", "depth_scaled", "year", "week_global", "x", "y"), names(base_df))

  pr <- predict(
    FIT,
    newdata = base_df[, needed, drop = FALSE],
    type = "link",
    re_form_iid = NA
  )

  eta <- pr$est_non_rf + pr$omega_s
  prob <- plogis(eta)

  r_prob <- r_template
  v <- rep(NA_real_, terra::ncell(r_template))
  v[base_df$cell] <- prob
  terra::values(r_prob) <- v

  list(label = label, raster = r_prob)
}

r_baseline <- fill_na_nearestish(load_sst_mean_on_grid(nc_path_baseline, 2000, 2019, r_template))
r_245_40   <- fill_na_nearestish(load_sst_mean_on_grid(nc_path_ssp245_2040, 2040, 2050, r_template))
r_585_40   <- fill_na_nearestish(load_sst_mean_on_grid(nc_path_ssp585_2040, 2040, 2050, r_template))
r_245_90   <- fill_na_nearestish(load_sst_mean_on_grid(nc_path_ssp245_2090, 2090, 2100, r_template))
r_585_90   <- fill_na_nearestish(load_sst_mean_on_grid(nc_path_ssp585_2090, 2090, 2100, r_template))

# Interpolate observed SST to the projection grid for a present-day reference
build_measured_sst_r <- function() {
  pts <- tibble(x = obs_sdm2$x, y = obs_sdm2$y, SST = obs_sdm2$SST) %>%
    filter(is.finite(x), is.finite(y), is.finite(SST))

  if (nrow(pts) < 10) return(NULL)

  m <- mgcv::gam(SST ~ s(x, y, bs = "tp"), data = pts, method = "REML")

  xy_m  <- as.data.frame(terra::xyFromCell(r_template, 1:terra::ncell(r_template)))
  xy_km <- transmute(xy_m, x = x / 1000, y = y / 1000)

  pred <- as.numeric(predict(m, newdata = xy_km))
  terra::setValues(r_template, pred)
}

r_measured <- build_measured_sst_r()

p_measured <- if (!is.null(r_measured)) {
  predict_from_sst_r(r_measured, "Obs-based SST (interpolated)")
} else {
  NULL
}

p_baseline <- predict_from_sst_r(r_baseline, "Modeled SST 2000–2019")
p_245_40   <- predict_from_sst_r(r_245_40,   "SSP2-4.5 2040–2050")
p_585_40   <- predict_from_sst_r(r_585_40,   "SSP5-8.5 2040–2050")
p_245_90   <- predict_from_sst_r(r_245_90,   "SSP2-4.5 2090–2100")
p_585_90   <- predict_from_sst_r(r_585_90,   "SSP5-8.5 2090–2100")

pred_list <- list(p_measured, p_baseline, p_245_40, p_585_40, p_245_90, p_585_90) %>%
  purrr::compact()

walk(pred_list, function(p) {
  fn <- file.path(pred_dir, paste0("HW_prob_", gsub("[^0-9A-Za-z]+", "_", p$label), ".tif"))
  terra::writeRaster(p$raster, filename = fn, overwrite = TRUE)
})

pred_list <- lapply(pred_list, function(p) {
  p$raster <- terra::mask(p$raster, land_r, inverse = TRUE)
  p
})

plot_df <- purrr::map_dfr(pred_list, function(p) {
  xy  <- as.data.frame(terra::xyFromCell(p$raster, 1:terra::ncell(p$raster)))
  val <- as.numeric(terra::values(p$raster))
  tibble(x = xy$x, y = xy$y, probability = val, Scenario = p$label) %>%
    filter(is.finite(probability))
})

plot_df$Scenario <- recode(
  plot_df$Scenario,
  "Modeled SST 2000–2019"        = "Modeled SST\n2000–2019",
  "SSP2-4.5 2040–2050"           = "SSP2-4.5\n2040–2050",
  "SSP2-4.5 2090–2100"           = "SSP2-4.5\n2090–2100",
  "SSP5-8.5 2040–2050"           = "SSP5-8.5\n2040–2050",
  "SSP5-8.5 2090–2100"           = "SSP5-8.5\n2090–2100",
  "Obs-based SST (interpolated)" = "Obs-based SST\n(interpolated)"
)

plot_df$Scenario <- factor(
  plot_df$Scenario,
  levels = c(
    "Obs-based SST\n(interpolated)",
    "Modeled SST\n2000–2019",
    "SSP2-4.5\n2040–2050",
    "SSP5-8.5\n2040–2050",
    "SSP2-4.5\n2090–2100",
    "SSP5-8.5\n2090–2100"
  )
)

xlim_m <- c(bbox_km$xmin, bbox_km$xmax) * 1000
ylim_m <- c(bbox_km$ymin, bbox_km$ymax) * 1000

p_scen <- ggplot() +
  geom_raster(data = plot_df, aes(x = x, y = y, fill = probability), interpolate = TRUE) +
  geom_sf(data = land_sf_m, fill = "grey70", color = "black", linewidth = 0.2, inherit.aes = FALSE) +
  geom_sf(data = park_utm, fill = NA, color = "#39b54a", linewidth = 0.7, inherit.aes = FALSE) +
  geom_sf(data = sanctuary_utm, fill = NA, color = "#ffd54f", linewidth = 0.7, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_viridis_c(
    option = "magma",
    name = "HW presence\nprobability",
    limits = c(0, 0.8),
    oob = scales::squish
  ) +
  coord_sf(
    xlim = xlim_m,
    ylim = ylim_m,
    crs = 32736,
    default_crs = st_crs(32736),
    expand = FALSE,
    datum = NA
  ) +
  scale_x_continuous("Easting (km)", labels = ~ .x / 1000) +
  scale_y_continuous("Northing (km)", labels = ~ .x / 1000) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold")) +
  facet_wrap(~ Scenario, nrow = 1)

ggsave(file.path(fig_dir, "HW_occurrence_scenarios.png"), p_scen, width = 14, height = 6, dpi = 300)


# ==============================================================================
# 11) Projection uncertainty: area above thresholds
# ==============================================================================
nsim_uncert <- 100L
thresh_vec  <- c(0.5, 0.1)
cell_km2    <- grid_res_km^2

scenario_rasters <- list(
  "Measured SST (current)" = r_measured,
  "Modeled SST 2000–2019"  = r_baseline,
  "SSP2-4.5 2040–2050"     = r_245_40,
  "SSP5-8.5 2040–2050"     = r_585_40,
  "SSP2-4.5 2090–2100"     = r_245_90,
  "SSP5-8.5 2090–2100"     = r_585_90
) %>%
  purrr::compact()

make_nd_from_sst <- function(sst_r) {
  xy_m <- as.data.frame(terra::xyFromCell(r_template, 1:terra::ncell(r_template)))
  sst  <- as.numeric(terra::values(sst_r))

  nd <- tibble(
    cell = 1:terra::ncell(r_template),
    x    = xy_m$x / 1000,
    y    = xy_m$y / 1000,
    SST  = sst
  ) %>%
    filter(is.finite(SST))

  nd$SST_scaled   <- (nd$SST - mu_SST) / sd_SST
  nd$depth_scaled <- depth_mean_scaled
  nd$week_global  <- as.integer(round(median(obs_sdm2$week_global, na.rm = TRUE)))
  nd$year         <- factor(levels(obs_sdm2$year)[1], levels = levels(obs_sdm2$year))

  nd
}

area_quantiles_from_scenario <- function(label, sst_r, model, nsim, thresholds, seed = 42L) {
  keep_ocean <- is.na(terra::values(land_r))
  if (!any(keep_ocean)) {
    return(tibble(Scenario = label, Threshold = thresholds, `2.5%` = 0, `50%` = 0, `97.5%` = 0))
  }

  nd <- make_nd_from_sst(sst_r)
  nd <- nd[keep_ocean[nd$cell], , drop = FALSE]

  if (!nrow(nd)) {
    return(tibble(Scenario = label, Threshold = thresholds, `2.5%` = 0, `50%` = 0, `97.5%` = 0))
  }

  needed <- c("SST_scaled", "depth_scaled", "year", "week_global", "x", "y")
  pr_parts <- predict(model, newdata = nd[, needed, drop = FALSE], type = "link", re_form_iid = NA)
  omega_hat <- pr_parts$omega_s

  set.seed(seed)
  sim_fe <- predict(
    model,
    newdata = nd[, needed, drop = FALSE],
    type = "link",
    nsim = nsim,
    re_form = NA,
    re_form_iid = NA
  )

  sim_mat <- extract_sims_matrix(sim_fe, nsim)
  eta_sim <- sim_mat + omega_hat
  p_sim   <- plogis(eta_sim)

  purrr::map_dfr(thresholds, function(th) {
    n_cells <- apply(p_sim, 2, function(col) sum(col > th, na.rm = TRUE))
    areas   <- n_cells * cell_km2
    qs <- quantile(areas, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    tibble(
      Scenario = label,
      Threshold = th,
      `2.5%` = unname(qs[[1]]),
      `50%`  = unname(qs[[2]]),
      `97.5%` = unname(qs[[3]])
    )
  })
}

uncertainty_tbl <- purrr::imap_dfr(
  scenario_rasters,
  ~ area_quantiles_from_scenario(.y, .x, FIT, nsim_uncert, thresh_vec)
)

write_csv(uncertainty_tbl, file.path(tab_dir, "HW_area_uncertainty.csv"))


# ==============================================================================
# 12) Mesh cutoff sensitivity (1, 2, and 5 km)
# ==============================================================================
cutoffs_km <- c(5, 2, 1)

build_barrier_mesh <- function(cutoff_km) {
  m0 <- sdmTMB::make_mesh(
    data    = obs_sdm2,
    xy_cols = c("x", "y"),
    cutoff  = cutoff_km
  )
  sdmTMBextra::add_barrier_mesh(
    m0, mz_utm_km,
    range_fraction = 0.1,
    proj_scaling   = 1000,
    plot           = FALSE
  )
}

fit_on_mesh <- function(mesh_obj) {
  sdmTMB(
    observed ~ s(depth_scaled, bs = "ts") + s(SST_scaled, bs = "ts") + (1 | year),
    data           = obs_sdm2,
    mesh           = mesh_obj,
    time           = "week_global",
    extra_time     = all_weeks,
    family         = binomial(link = "logit"),
    spatial        = "on",
    spatiotemporal = "ar1",
    reml           = TRUE,
    control        = sdmTMB::sdmTMBcontrol(eval.max = 4000L, iter.max = 2000L)
  )
}

mesh_list <- map(cutoffs_km, build_barrier_mesh)
names(mesh_list) <- paste0("cut", cutoffs_km, "km")

fits <- list()
for (nm in names(mesh_list)) {
  message("Fitting ", nm, " ...")
  fits[[nm]] <- tryCatch(
    fit_on_mesh(mesh_list[[nm]]),
    error = function(e) {
      message("  FAILED: ", conditionMessage(e))
      NULL
    }
  )
  gc()
}
fits <- purrr::compact(fits)

if (length(fits) > 0) {
  pred_probs <- imap(fits, function(fit, nm) {
    tryCatch(
      predict(fit, newdata = obs_sdm2, type = "response")$est,
      error = function(e) {
        message("Predict failed for ", nm, ": ", conditionMessage(e))
        NULL
      }
    )
  }) %>%
    purrr::compact()

  auc_tbl <- imap_dfr(pred_probs, function(p, nm) {
    a <- pROC::auc(obs_sdm2$observed, p)
    tibble(mesh = nm, AUC = as.numeric(a))
  })

  get_ran_pars <- function(fit) {
    rp <- sdmTMB::tidy(fit, effects = "ran_pars") %>%
      mutate(
        term = recode(
          term,
          "sigma_Omega"   = "sigma_O",
          "sigma_Epsilon" = "sigma_E",
          "ar1_phi"       = "rho"
        )
      )
    out <- rp %>%
      select(term, estimate) %>%
      tidyr::pivot_wider(names_from = term, values_from = estimate)

    for (nm in c("range", "sigma_O", "sigma_E", "rho")) {
      if (!nm %in% names(out)) out[[nm]] <- NA_real_
    }
    out
  }

  get_n_knots <- function(mesh_obj) {
    val <- NA_integer_
    try(val <- mesh_obj$mesh$n, silent = TRUE)
    if (is.na(val)) try(val <- nrow(mesh_obj$mesh$loc), silent = TRUE)
    val
  }

  param_tbl <- imap_dfr(fits, function(fit, nm) {
    rp <- get_ran_pars(fit)
    tibble(
      mesh      = nm,
      AIC       = AIC(fit),
      n_knots   = get_n_knots(mesh_list[[nm]]),
      range_km  = rp$range,
      sigma_O   = rp$sigma_O,
      sigma_E   = rp$sigma_E,
      rho       = rp$rho
    )
  }) %>%
    left_join(auc_tbl, by = "mesh") %>%
    mutate(cutoff_km = as.numeric(gsub("^cut|km$", "", mesh))) %>%
    arrange(cutoff_km)

  write_csv(param_tbl, file.path(tab_dir, "mesh_cutoff_sensitivity.csv"))
}

save.image(file = WS_FILE)
cat("Workspace path:", normalizePath(WS_FILE, winslash = "/", mustWork = FALSE), "\n")
