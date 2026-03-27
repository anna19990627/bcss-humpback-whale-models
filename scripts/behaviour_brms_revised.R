
################################################################################
# Behavioural state ~ environment (brms categorical, multinomial)
# BCSS humpback whale surveys, 2020–2024
#
# Aim:
#   Quantify how behavioural state (R, S, T) varies with key environmental and
#   social covariates using a Bayesian multinomial logistic regression.
#
# Workflow:
#   A) Screen social variables on the largest available subsets
#   B) Rank environmental models on a consistent complete-case dataset
#   C) Check spatial sensitivity with a coarse 10-km cell random effect
#   D) Save final-model diagnostics
#   E) Plot predicted probabilities along unscaled gradients
#   F) Export supplement tables
#
# Outputs:
#   figs/behaviour/
#   tables/behaviour/
#   models_behaviour/
################################################################################


# ==============================================================================
# 0) Packages and global settings
# ==============================================================================
library(tidyverse)
library(sf)
library(lubridate)
library(brms)
library(loo)
library(posterior)
library(scales)
library(knitr)
library(readr)
library(patchwork)
library(spdep)

set.seed(123)
options(mc.cores = parallel::detectCores())

REFIT_MODELS <- FALSE


# ==============================================================================
# 1) Paths and output folders
# ==============================================================================
data_path <- "clean_dfs/hw.rds"

fig_dir   <- file.path("figs",   "behaviour")
table_dir <- file.path("tables", "behaviour")
model_dir <- "models_behaviour"

dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

WS_FILE <- file.path(model_dir, "behaviour_HW.RData")
if (file.exists(WS_FILE)) load(WS_FILE)


# ==============================================================================
# 2) Load and harmonize data
# ==============================================================================
hw_sf <- readRDS(data_path)

# Harmonize known column-name typos
hw_sf <- hw_sf %>%
  rename_with(~ "current_speed",   .cols = any_of("currend_speed")) %>%
  rename_with(~ "max_wave_height", .cols = any_of("max_wave_hight"))

# Non-spatial working data frame
hw_df <- hw_sf %>%
  st_drop_geometry() %>%
  mutate(
    b_state        = factor(b_state),
    group_cohesion = factor(group_cohesion),
    year           = factor(year)
  )


# ==============================================================================
# 3) brms fitting helper
# ==============================================================================
ctrl_behav <- list(adapt_delta = 0.99)

fit_cat <- function(formula, data,
                    knots = NULL,
                    file = NULL,
                    refit = FALSE,
                    chains = 4,
                    cores  = 4,
                    iter   = 2000,
                    warmup = 1000) {

  if (!is.null(file) && file.exists(file) && !refit) {
    message("Reading cached model: ", file)
    return(readRDS(file))
  }

  fit <- brm(
    formula   = formula,
    family    = categorical(),
    data      = data,
    knots     = knots,
    chains    = chains,
    cores     = cores,
    iter      = iter,
    warmup    = warmup,
    control   = ctrl_behav,
    seed      = 123,
    save_pars = save_pars(all = TRUE)
  )

  if (!is.null(file)) saveRDS(fit, file)
  fit
}

ensure_loo <- function(fit) {
  if (inherits(fit, "brmsfit") && is.null(fit$criteria$loo)) add_criterion(fit, "loo") else fit
}


# ==============================================================================
# A) Social predictors: do they justify losing observations?
# ==============================================================================
# A1) Group cohesion
dat_gc <- hw_df %>%
  filter(!is.na(b_state), !is.na(year), !is.na(group_cohesion))

m_behav_gc <- fit_cat(
  b_state ~ group_cohesion + (1 | year),
  data  = dat_gc,
  file  = file.path(model_dir, "A1_m_behav_groupcohesion.rds"),
  refit = REFIT_MODELS
)

m_behav_gc_null <- fit_cat(
  b_state ~ 1 + (1 | year),
  data  = dat_gc,
  file  = file.path(model_dir, "A1_m_behav_groupcohesion_null.rds"),
  refit = REFIT_MODELS
)

loo_gc <- loo(m_behav_gc, m_behav_gc_null)
print(loo_gc)

# A2) Group size (QTY)
dat_qty <- hw_df %>%
  filter(!is.na(b_state), !is.na(year), !is.na(QTY))

m_behav_qty <- fit_cat(
  b_state ~ QTY + (1 | year),
  data  = dat_qty,
  file  = file.path(model_dir, "A2_m_behav_QTY_socialcheck.rds"),
  refit = REFIT_MODELS
)

m_behav_qty_null <- fit_cat(
  b_state ~ 1 + (1 | year),
  data  = dat_qty,
  file  = file.path(model_dir, "A2_m_behav_QTY_socialcheck_null.rds"),
  refit = REFIT_MODELS
)

loo_qty <- loo(m_behav_qty, m_behav_qty_null)
print(loo_qty)


# ==============================================================================
# B) Environmental model selection on a consistent complete-case dataset
# ==============================================================================
# B1) Build the analysis dataset
hw_behav <- hw_df %>%
  filter(
    !is.na(b_state),
    !is.na(year),
    !is.na(Date_Time),

    !is.na(QTY),
    !is.na(depth),
    !is.na(dist_m),
    !is.na(SST),

    !is.na(weeknr),
    !is.na(wind_speed),
    !is.na(max_wave_height),
    !is.na(current_speed)
  ) %>%
  mutate(
    ToD = hour(Date_Time) + minute(Date_Time) / 60 + second(Date_Time) / 3600,
    SST_s             = as.numeric(scale(SST)),
    depth_s           = as.numeric(scale(depth)),
    dist_m_s          = as.numeric(scale(dist_m)),
    wind_speed_s      = as.numeric(scale(wind_speed)),
    max_wave_height_s = as.numeric(scale(max_wave_height)),
    current_speed_s   = as.numeric(scale(current_speed)),
    weeknr_s          = as.numeric(scale(weeknr))
  )

behav_meta <- hw_behav %>%
  summarise(
    N        = n(),
    n_R      = sum(b_state == "R"),
    n_S      = sum(b_state == "S"),
    n_T      = sum(b_state == "T"),
    mean_SST = mean(SST, na.rm = TRUE),
    sd_SST   = sd(SST, na.rm = TRUE),
    mean_QTY = mean(QTY, na.rm = TRUE),
    mean_mwh = mean(max_wave_height, na.rm = TRUE),
    sd_mwh   = sd(max_wave_height, na.rm = TRUE),
    n_year   = n_distinct(year)
  )
print(behav_meta)
write_csv(behav_meta, file.path(table_dir, "behaviour_dataset_summary.csv"))

knots_behav <- list(ToD = c(0, 24))

# B2) Candidate models
m_sst <- fit_cat(
  b_state ~ SST_s + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B2_m_sst.rds"),
  refit = REFIT_MODELS
)

m_qty <- fit_cat(
  b_state ~ QTY + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B2_m_qty.rds"),
  refit = REFIT_MODELS
)

m_mwh <- fit_cat(
  b_state ~ max_wave_height_s + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B2_m_mwh.rds"),
  refit = REFIT_MODELS
)

m_week <- fit_cat(
  b_state ~ weeknr + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B2_m_week.rds"),
  refit = REFIT_MODELS
)

m_tod_lin <- fit_cat(
  b_state ~ ToD + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B2_m_tod_lin.rds"),
  refit = REFIT_MODELS
)

m_tod_cc <- fit_cat(
  b_state ~ s(ToD, bs = "cc", k = 8) + (1 | year),
  data  = hw_behav,
  knots = knots_behav,
  file  = file.path(model_dir, "B2_m_tod_cc.rds"),
  refit = REFIT_MODELS
)

m_wind <- fit_cat(
  b_state ~ wind_speed_s + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B2_m_wind.rds"),
  refit = REFIT_MODELS
)

m_current <- fit_cat(
  b_state ~ current_speed_s + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B2_m_current.rds"),
  refit = REFIT_MODELS
)

m_depth <- fit_cat(
  b_state ~ depth_s + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B2_m_depth.rds"),
  refit = REFIT_MODELS
)

m_dist <- fit_cat(
  b_state ~ dist_m_s + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B2_m_dist.rds"),
  refit = REFIT_MODELS
)

# B3) Hypothesis-driven combined models
m_best1 <- m_sst

m_best2 <- fit_cat(
  b_state ~ SST_s + QTY + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B3_m_best2.rds"),
  refit = REFIT_MODELS
)

m_best3 <- fit_cat(
  b_state ~ SST_s + QTY + max_wave_height_s + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B3_m_best3.rds"),
  refit = REFIT_MODELS
)

m_best4 <- fit_cat(
  b_state ~ SST_s + QTY + max_wave_height_s + weeknr + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B3_m_best4.rds"),
  refit = REFIT_MODELS
)

m_best5 <- fit_cat(
  b_state ~ s(ToD, bs = "cc", k = 8) + SST_s + QTY + max_wave_height_s + weeknr + (1 | year),
  data  = hw_behav,
  knots = knots_behav,
  file  = file.path(model_dir, "B3_m_best5.rds"),
  refit = REFIT_MODELS
)

m_weather <- fit_cat(
  b_state ~ max_wave_height_s + wind_speed_s + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B3_m_weather.rds"),
  refit = REFIT_MODELS
)

m_location <- fit_cat(
  b_state ~ dist_m_s + depth_s + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B3_m_location.rds"),
  refit = REFIT_MODELS
)

m_full <- fit_cat(
  b_state ~ SST_s + QTY + max_wave_height_s + weeknr +
    wind_speed_s + depth_s + dist_m_s + current_speed_s + (1 | year),
  data  = hw_behav,
  file  = file.path(model_dir, "B3_m_full.rds"),
  refit = REFIT_MODELS
)

# B4) PSIS-LOO model-selection table
mod_list_behav <- list(
  m_best3    = m_best3,
  m_best5    = m_best5,
  m_best4    = m_best4,
  m_best2    = m_best2,
  m_best1    = m_best1,

  m_full     = m_full,
  m_weather  = m_weather,
  m_location = m_location,

  m_qty      = m_qty,
  m_mwh      = m_mwh,
  m_week     = m_week,
  m_tod_cc   = m_tod_cc,
  m_tod_lin  = m_tod_lin,
  m_wind     = m_wind,
  m_depth    = m_depth,
  m_dist     = m_dist,
  m_current  = m_current
)

mod_list_behav <- purrr::map(mod_list_behav, ensure_loo)

loo_tbl <- imap_dfr(mod_list_behav, \(fit, nm) {
  est <- fit$criteria$loo$estimates
  tibble(
    model    = nm,
    elpd_loo = est["elpd_loo", "Estimate"],
    se_elpd  = est["elpd_loo", "SE"],
    LOOIC    = est["looic", "Estimate"],
    se_LOOIC = est["looic", "SE"],
    p_loo    = est["p_loo", "Estimate"]
  )
})

min_LOOIC <- min(loo_tbl$LOOIC)
loo_tbl <- loo_tbl %>%
  mutate(
    dLOOIC  = LOOIC - min_LOOIC,
    rel_lik = exp(-0.5 * dLOOIC),
    w_LOO   = rel_lik / sum(rel_lik)
  ) %>%
  arrange(dLOOIC)

model_specs <- tribble(
  ~model,      ~Model_label,                         ~Predictors,
  "m_best3",   "SST + QTY + max wave (final)",      "SST, group size (QTY), maximum wave height",
  "m_best5",   "m_best3 + week + s(ToD)",           "SST, QTY, max wave height, week of year, time of day",
  "m_best4",   "m_best3 + week",                    "SST, QTY, max wave height, week of year",
  "m_best2",   "SST + QTY",                         "SST, QTY",
  "m_best1",   "SST only",                          "SST",

  "m_full",    "Full env. (no tide)",               "SST, QTY, max wave height, week, wind, depth, distance, current",
  "m_weather", "Max wave + wind",                   "Max wave height, wind speed",
  "m_location","Depth + distance",                  "Depth, distance to shore",

  "m_qty",     "QTY only",                          "QTY",
  "m_mwh",     "Max wave height only",              "Max wave height",
  "m_week",    "Week of year only",                 "Week of year",
  "m_tod_cc",  "Time of day only (cyclic)",         "Time of day",
  "m_tod_lin", "Time of day only (linear)",         "Time of day",
  "m_wind",    "Wind speed only",                   "Wind speed",
  "m_depth",   "Depth only",                        "Depth",
  "m_dist",    "Distance to shore only",            "Distance to shore",
  "m_current", "Current speed only",                "Current speed"
)

tab_models <- loo_tbl %>%
  left_join(model_specs, by = "model") %>%
  mutate(
    Model_label = ifelse(is.na(Model_label), model, Model_label),
    Predictors  = ifelse(is.na(Predictors),  "—", Predictors)
  ) %>%
  transmute(
    Model      = Model_label,
    Predictors = Predictors,
    LOOIC      = round(LOOIC, 2),
    dLOOIC     = round(dLOOIC, 2),
    w_LOO      = round(w_LOO, 2),
    p_loo      = round(p_loo, 2)
  )

kable(
  tab_models,
  booktabs = TRUE,
  caption = "Table Sx. Behavioural-state model selection using PSIS-LOO (ranked by increasing LOOIC). All models include a year random intercept."
)
write_csv(tab_models, file.path(table_dir, "TableSx_behaviour_model_selection.csv"))

m_behav_final <- mod_list_behav[["m_best3"]]


# ==============================================================================
# C) Spatial sensitivity check: 10-km grid-cell random effect
# ==============================================================================
behav_cell_km <- 10

# Rebuild the final-model dataset with UTM coordinates and a 10-km cell id
hw_behav_sp <- hw_sf %>%
  st_transform(32736) %>%
  mutate(
    x_utm_km = st_coordinates(.)[, 1] / 1000,
    y_utm_km = st_coordinates(.)[, 2] / 1000
  ) %>%
  st_drop_geometry() %>%
  mutate(
    b_state = factor(b_state),
    year    = factor(year),
    cell_x  = floor(x_utm_km / behav_cell_km),
    cell_y  = floor(y_utm_km / behav_cell_km),
    cell    = factor(paste(cell_x, cell_y, sep = "_"))
  ) %>%
  filter(
    !is.na(b_state),
    !is.na(year),
    !is.na(Date_Time),

    !is.na(QTY),
    !is.na(depth),
    !is.na(dist_m),
    !is.na(SST),

    !is.na(weeknr),
    !is.na(wind_speed),
    !is.na(max_wave_height),
    !is.na(current_speed)
  ) %>%
  mutate(
    SST_s             = (SST - behav_meta$mean_SST) / behav_meta$sd_SST,
    max_wave_height_s = (max_wave_height - behav_meta$mean_mwh) / behav_meta$sd_mwh
  )

m_best3_cell <- fit_cat(
  b_state ~ SST_s + QTY + max_wave_height_s + (1 | year) + (1 | cell),
  data  = hw_behav_sp,
  file  = file.path(model_dir, "C_m_best3_cell_10kmRE.rds"),
  refit = REFIT_MODELS
)

m_best3_loo      <- ensure_loo(m_behav_final)
m_best3_cell_loo <- ensure_loo(m_best3_cell)

cat("\nSpatial sensitivity (cell RE) LOO compare:\n")
print(loo_compare(m_best3_loo$criteria$loo, m_best3_cell_loo$criteria$loo))
summary(m_best3_cell)

# Moran's I on class-wise residuals (population-level predictions)
ep <- posterior_epred(m_behav_final, newdata = hw_behav_sp, re_formula = NA)
p_mean <- apply(ep, c(2, 3), mean)
classes <- dimnames(ep)[[3]]
colnames(p_mean) <- classes

coords <- as.matrix(hw_behav_sp[, c("x_utm_km", "y_utm_km")])

# Small jitter (1 m = 0.001 km) avoids identical-coordinate warnings
set.seed(1)
coords_j <- coords + matrix(rnorm(length(coords), sd = 0.001), ncol = 2)

knn_k <- 8
nb <- knn2nb(knearneigh(coords_j, k = knn_k))
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_results <- map_dfr(classes, function(cl) {
  obs01 <- as.integer(hw_behav_sp$b_state == cl)
  resid <- obs01 - p_mean[, cl]
  mt <- moran.test(resid, lw, zero.policy = TRUE)
  tibble(
    class   = cl,
    Moran_I = unname(mt$estimate[["Moran I statistic"]]),
    p_value = mt$p.value
  )
})

print(moran_results)
write_csv(moran_results, file.path(table_dir, "TableSx_behaviour_moranI.csv"))


# ==============================================================================
# D) Diagnostics for the final model
# ==============================================================================
loo_final <- loo(m_behav_final)
print(loo_final)

png(file.path(fig_dir, "diag_behaviour_loo.png"), width = 1400, height = 900, res = 180)
plot(loo_final)
dev.off()

np <- nuts_params(m_behav_final)
n_div <- sum(np$Parameter == "divergent__" & np$Value == 1)
n_treedepth <- sum(np$Parameter == "treedepth__" & np$Value >= 10)

cat("\nDivergences:", n_div, "\n")
cat("Max treedepth hits:", n_treedepth, "\n")

diag_tbl <- tibble(
  metric = c("divergences", "max_treedepth_hits"),
  value  = c(n_div, n_treedepth)
)
write_csv(diag_tbl, file.path(table_dir, "behaviour_sampler_diagnostics.csv"))

if (!is.null(loo_final$diagnostics$pareto_k)) {
  pareto_k_tbl <- tibble(
    observation = seq_along(loo_final$diagnostics$pareto_k),
    pareto_k    = loo_final$diagnostics$pareto_k
  )
  write_csv(pareto_k_tbl, file.path(table_dir, "behaviour_pareto_k.csv"))
}

pp_behav <- pp_check(m_behav_final, type = "bars", ndraws = 100)
ggsave(file.path(fig_dir, "diag_behaviour_ppcheck_bars.png"), pp_behav,
       width = 7.5, height = 4.8, dpi = 300)


# ==============================================================================
# E) Predicted probabilities along unscaled gradients
# ==============================================================================
behav_pal <- c(
  "R" = "#4F5D75",
  "S" = "#EF8354",
  "T" = "#BFC0C0"
)

plot_behav_prob_curve <- function(
    fit, data,
    var = c("SST", "QTY", "max_wave_height"),
    n = 200, ci = 0.95, re_formula = NA,
    trim = c(0.05, 0.95),
    xlab = NULL, title = NULL,
    ndraws = NULL
) {
  var <- match.arg(var)

  x_raw <- data[[var]]
  x_seq <- seq(
    quantile(x_raw, trim[1], na.rm = TRUE),
    quantile(x_raw, trim[2], na.rm = TRUE),
    length.out = n
  )

  mu_SST <- mean(data$SST, na.rm = TRUE)
  sd_SST <- sd(data$SST, na.rm = TRUE)
  mu_QTY <- mean(data$QTY, na.rm = TRUE)
  mu_MWH <- mean(data$max_wave_height, na.rm = TRUE)
  sd_MWH <- sd(data$max_wave_height, na.rm = TRUE)

  nd <- tibble(
    SST             = mu_SST,
    QTY             = mu_QTY,
    max_wave_height = mu_MWH,
    year            = factor(levels(data$year)[1], levels = levels(data$year))
  ) %>%
    slice(rep(1, n))

  nd[[var]] <- x_seq

  nd <- nd %>%
    mutate(
      SST_s             = (SST - mu_SST) / sd_SST,
      max_wave_height_s = (max_wave_height - mu_MWH) / sd_MWH
    )

  ep <- posterior_epred(fit, newdata = nd, re_formula = re_formula, ndraws = ndraws)

  alpha <- (1 - ci) / 2
  probs <- c(alpha, 0.5, 1 - alpha)
  summ  <- apply(ep, c(2, 3), \(v) quantile(v, probs = probs))

  class_names <- dimnames(ep)[[3]]
  if (is.null(class_names)) class_names <- levels(data$b_state)

  out <- bind_rows(lapply(seq_along(class_names), function(k) {
    tibble(
      x     = x_seq,
      state = factor(class_names[k], levels = class_names),
      lo    = summ[1, , k],
      mid   = summ[2, , k],
      hi    = summ[3, , k]
    )
  }))

  behav_lab <- c(
    "R" = "Resting (R)",
    "S" = "Social (S)",
    "T" = "Travelling (T)"
  )
  labs_map <- setNames(levels(out$state), levels(out$state))
  labs_map[names(behav_lab)] <- behav_lab

  if (is.null(xlab)) xlab <- var
  if (is.null(title)) title <- paste0("Behaviour state along a ", gsub("_", "-", var), " gradient")

  ggplot(out, aes(x = x, y = mid, ymin = lo, ymax = hi, colour = state, fill = state)) +
    geom_ribbon(alpha = 0.25, linewidth = 0) +
    geom_line(linewidth = 1.1) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_colour_manual(
      values = behav_pal,
      breaks = levels(out$state),
      labels = labs_map[levels(out$state)],
      drop   = FALSE
    ) +
    scale_fill_manual(
      values = behav_pal,
      breaks = levels(out$state),
      labels = labs_map[levels(out$state)],
      drop   = FALSE
    ) +
    labs(
      x = xlab,
      y = "Predicted probability",
      title = title,
      subtitle = paste0("Median and ", ci * 100, "% credible intervals (population-level)"),
      colour = "Behaviour state",
      fill   = "Behaviour state"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

p_sst <- plot_behav_prob_curve(
  m_behav_final, hw_behav, var = "SST",
  xlab  = "Sea surface temperature (°C)",
  title = "Behaviour state along an SST gradient"
)

p_qty <- plot_behav_prob_curve(
  m_behav_final, hw_behav, var = "QTY",
  xlab  = "Group size (QTY)",
  title = "Behaviour state along a group-size gradient"
)

p_mwh <- plot_behav_prob_curve(
  m_behav_final, hw_behav, var = "max_wave_height",
  xlab  = "Maximum wave height (m)",
  title = "Behaviour state along a wave-height gradient"
)

ggsave(file.path(fig_dir, "fig_behav_sst_prob.png"), p_sst, width = 7.5, height = 4.8, dpi = 400)
ggsave(file.path(fig_dir, "fig_behav_qty_prob.png"), p_qty, width = 7.5, height = 4.8, dpi = 400)
ggsave(file.path(fig_dir, "fig_behav_mwh_prob.png"), p_mwh, width = 7.5, height = 4.8, dpi = 400)

p_all <- (p_sst / p_qty / p_mwh) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "fig_behav_gradients_combined.png"), p_all,
       width = 8.5, height = 13.0, dpi = 400)


# ==============================================================================
# F) Supplement tables
# ==============================================================================
# F1) Odds ratios on original units
state_levels <- levels(m_behav_final$data$b_state)
ref_state    <- state_levels[1]

state_lab <- c(R = "Resting", S = "Social", T = "Travelling")
ref_name  <- ifelse(ref_state %in% names(state_lab), state_lab[[ref_state]], ref_state)

draws <- as_draws_df(m_behav_final)

get_beta <- function(nm) {
  if (!nm %in% names(draws)) stop("Missing coefficient: ", nm)
  draws[[nm]]
}

summ_or <- function(beta_draws, digits = 2) {
  OR <- exp(beta_draws)
  tibble(
    OR_median = round(median(OR), digits),
    lo95      = round(unname(quantile(OR, 0.025)), digits),
    hi95      = round(unname(quantile(OR, 0.975)), digits)
  )
}

sd_SST <- behav_meta$sd_SST
sd_MWH <- behav_meta$sd_mwh

pred_spec <- tribble(
  ~Predictor,                      ~coef_suffix,             ~Unit,           ~conv,
  "Sea-surface temperature (SST)", "SST_s",                  "per 1 °C",      1 / sd_SST,
  "Maximum wave height",           "max_wave_height_s",      "per 1 m",       1 / sd_MWH,
  "Group size (QTY)",              "QTY",                    "per +1 whale",  1
)

nonref_states <- state_levels[-1]

or_tbl <- map_dfr(nonref_states, function(st) {
  st_name  <- ifelse(st %in% names(state_lab), state_lab[[st]], st)
  contrast <- paste0(st_name, " vs ", ref_name)
  prefix   <- paste0("b_mu", st, "_")

  map_dfr(seq_len(nrow(pred_spec)), function(i) {
    p <- pred_spec[i, ]
    beta_unit <- get_beta(paste0(prefix, p$coef_suffix)) * p$conv
    s <- summ_or(beta_unit, digits = 2)
    tibble(
      Contrast = contrast,
      Predictor = p$Predictor,
      Unit = p$Unit,
      `Odds ratio (median)` = s$OR_median,
      `95% CrI` = paste0(s$lo95, "–", s$hi95)
    )
  })
})

caption_or <- paste0(
  "Table Sx. Behavioural-state odds ratios from a Bayesian multinomial logistic regression ",
  "(brms; family = categorical; N = ", behav_meta$N, "; year random intercept with ", behav_meta$n_year, " levels). ",
  ref_name, " was the reference category. Odds ratios are expressed on original units ",
  "(per 1 °C SST, per 1 m maximum wave height, per +1 whale in group size). ",
  "SST and maximum wave height were standardized in the fitted model and back-transformed using SDs from the analysis dataset ",
  "(SD_SST = ", round(sd_SST, 3), " °C; SD_wave = ", round(sd_MWH, 3), " m)."
)

kable(or_tbl, booktabs = TRUE, caption = caption_or)
write_csv(or_tbl, file.path(table_dir, "TableSx_behaviour_OR.csv"))

# F2) Predicted probabilities under baseline and simple perturbations
mu_SST <- behav_meta$mean_SST
mu_QTY <- behav_meta$mean_QTY
mu_MWH <- behav_meta$mean_mwh
ref_year <- factor(levels(hw_behav$year)[1], levels = levels(hw_behav$year))

nd_scen <- tibble(
  Scenario = c("Baseline (means)", "SST +1°C", "QTY +1", "Max wave +0.5 m"),
  SST = c(mu_SST, mu_SST + 1, mu_SST, mu_SST),
  QTY = c(mu_QTY, mu_QTY, mu_QTY + 1, mu_QTY),
  max_wave_height = c(mu_MWH, mu_MWH, mu_MWH, mu_MWH + 0.5),
  year = ref_year
) %>%
  mutate(
    SST_s             = (SST - mu_SST) / sd_SST,
    max_wave_height_s = (max_wave_height - mu_MWH) / sd_MWH
  )

ep_scen <- posterior_epred(m_behav_final, newdata = nd_scen, re_formula = NA)
classes <- dimnames(ep_scen)[[3]]
if (is.null(classes)) classes <- levels(hw_behav$b_state)

summ_ci <- function(x, digits = 3) {
  tibble(
    lo95 = round(unname(quantile(x, 0.025)), digits),
    med  = round(unname(quantile(x, 0.50)),  digits),
    hi95 = round(unname(quantile(x, 0.975)), digits)
  )
}

abs_tbl <- bind_rows(lapply(seq_len(nrow(nd_scen)), function(i) {
  tmp <- bind_rows(lapply(seq_along(classes), function(k) {
    s <- summ_ci(ep_scen[, i, k], digits = 3)
    tibble(State = classes[k], lo95 = s$lo95, med = s$med, hi95 = s$hi95)
  }))
  tmp$Scenario <- nd_scen$Scenario[i]
  tmp
})) %>%
  select(Scenario, State, med, lo95, hi95)

abs_tbl_wide <- abs_tbl %>%
  mutate(val = paste0(med, " [", lo95, ", ", hi95, "]")) %>%
  select(Scenario, State, val) %>%
  pivot_wider(names_from = State, values_from = val)

caption_prob <- paste0(
  "Table Sy. Predicted behavioural-state probabilities (median [95% CrI]) at typical conditions and under simple perturbations. ",
  "Baseline uses predictor means from the analysis dataset (SST = ", round(mu_SST, 2), " °C; QTY = ", round(mu_QTY, 2),
  "; max wave height = ", round(mu_MWH, 2), " m). Perturbations change one predictor at a time ",
  "(SST +1 °C; QTY +1; max wave +0.5 m) holding other predictors at baseline values. ",
  "Predictions are population-level (year random effects excluded)."
)

kable(abs_tbl_wide, booktabs = TRUE, caption = caption_prob)
write_csv(abs_tbl_wide, file.path(table_dir, "TableSy_behaviour_probabilities.csv"))

message("DONE. Behaviour outputs written to: ", fig_dir, " and ", table_dir, ". Cached models in: ", model_dir)

save.image(file = WS_FILE)
cat("Workspace path:", normalizePath(WS_FILE, winslash = "/", mustWork = FALSE), "\n")
