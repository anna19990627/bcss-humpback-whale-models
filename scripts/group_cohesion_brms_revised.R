
################################################################################
# Group cohesion ~ environment (brms categorical, multinomial)
# BCSS humpback whale surveys, 2020–2024
#
# Aim:
#   Quantify how humpback whale group-cohesion classes vary with environmental
#   gradients using a Bayesian multinomial model.
#
# Workflow:
#   1) Prepare the analysis dataset and pool sparse classes if needed
#   2) Fit candidate models and rank them with PSIS-LOO
#   3) Check spatial sensitivity with a 10-km cell random effect
#   4) Save final-model diagnostics
#   5) Plot predicted probabilities along depth and SST gradients
#   6) Export supplement tables
#
# Outputs:
#   figs/group_cohesion/
#   tables/group_cohesion/
#   models_gc/
################################################################################


# ==============================================================================
# 0) Packages and global settings
# ==============================================================================
library(tidyverse)   # dplyr, ggplot2, purrr, tibble, tidyr, forcats
library(sf)
library(lubridate)
library(brms)
library(loo)
library(knitr)
library(scales)
library(posterior)
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

fig_dir   <- file.path("figs",   "group_cohesion")
table_dir <- file.path("tables", "group_cohesion")
model_dir <- "models_gc"

dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

WS_FILE <- file.path(model_dir, "group_cohesion_HW.RData")
if (file.exists(WS_FILE)) load(WS_FILE)


# ==============================================================================
# 2) Data preparation
# ==============================================================================
hw_sf <- readRDS(data_path)

# Harmonize known column-name typos
hw_sf <- hw_sf %>%
  rename_with(~ "current_speed",   .cols = any_of("currend_speed")) %>%
  rename_with(~ "max_wave_height", .cols = any_of("max_wave_hight"))

# Base modelling data set with UTM coordinates in km for diagnostics
hw_gc_base <- hw_sf %>%
  st_transform(32736) %>%
  mutate(
    x_utm_km = st_coordinates(.)[, 1] / 1000,
    y_utm_km = st_coordinates(.)[, 2] / 1000
  ) %>%
  st_drop_geometry() %>%
  mutate(
    group_cohesion = factor(group_cohesion),
    year           = factor(year),
    ToD            = hour(Date_Time) + minute(Date_Time) / 60 + second(Date_Time) / 3600
  ) %>%
  filter(
    !is.na(group_cohesion), !is.na(year),
    !is.na(depth), depth > 0,
    !is.na(dist_m),
    !is.na(SST),
    !is.na(wind_speed),
    !is.na(max_wave_height),
    !is.na(current_speed),
    !is.na(weeknr),
    !is.na(ToD)
  ) %>%
  mutate(
    depth_s           = as.numeric(scale(depth)),
    dist_m_s          = as.numeric(scale(dist_m)),
    SST_s             = as.numeric(scale(SST)),
    wind_speed_s      = as.numeric(scale(wind_speed)),
    max_wave_height_s = as.numeric(scale(max_wave_height)),
    current_speed_s   = as.numeric(scale(current_speed)),
    weeknr_s          = as.numeric(scale(weeknr))
  )

# Pool sparse categories only when needed
min_n <- 20
keep_levels <- names(which(table(hw_gc_base$group_cohesion) >= min_n))

hw_gc <- hw_gc_base %>%
  mutate(
    group_cohesion = forcats::fct_other(
      group_cohesion,
      keep        = keep_levels,
      other_level = "rare"
    )
  ) %>%
  droplevels()

knots_gc <- list(ToD = c(0, 24))

# 10-km grid-cell id for the spatial sensitivity model
gc_cell_km <- 10
hw_gc <- hw_gc %>%
  mutate(
    cell_x = floor(x_utm_km / gc_cell_km),
    cell_y = floor(y_utm_km / gc_cell_km),
    cell   = factor(paste(cell_x, cell_y, sep = "_"))
  )

gc_counts <- hw_gc %>% count(group_cohesion, name = "n") %>% arrange(desc(n))
gc_meta <- hw_gc %>%
  summarise(
    N          = n(),
    mean_depth = mean(depth, na.rm = TRUE),
    sd_depth   = sd(depth,   na.rm = TRUE),
    mean_SST   = mean(SST,   na.rm = TRUE),
    sd_SST     = sd(SST,     na.rm = TRUE),
    n_year     = n_distinct(year),
    n_cell     = n_distinct(cell)
  )

print(gc_counts)
print(gc_meta)
write_csv(gc_counts, file.path(table_dir, "groupcohesion_class_counts.csv"))
write_csv(gc_meta,   file.path(table_dir, "groupcohesion_dataset_summary.csv"))


# ==============================================================================
# 3) brms fitting helper
# ==============================================================================
ctrl_sp <- list(adapt_delta = 0.999, max_treedepth = 14)

fit_cat <- function(formula, data,
                    knots = NULL,
                    file = NULL,
                    refit = FALSE,
                    chains = 4,
                    cores  = 4,
                    iter   = 4000,
                    warmup = 2000) {

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
    control   = ctrl_sp,
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
# 4) Candidate models
# ==============================================================================
m0_gc <- fit_cat(
  group_cohesion ~ 1 + (1 | year),
  data  = hw_gc,
  file  = file.path(model_dir, "m0_gc.rds"),
  refit = REFIT_MODELS
)

m_gc_depth <- fit_cat(
  group_cohesion ~ depth_s + (1 | year),
  data  = hw_gc,
  file  = file.path(model_dir, "m_gc_depth.rds"),
  refit = REFIT_MODELS
)

m_gc_dist <- fit_cat(
  group_cohesion ~ dist_m_s + (1 | year),
  data  = hw_gc,
  file  = file.path(model_dir, "m_gc_dist.rds"),
  refit = REFIT_MODELS
)

m_gc_sst <- fit_cat(
  group_cohesion ~ SST_s + (1 | year),
  data  = hw_gc,
  file  = file.path(model_dir, "m_gc_sst.rds"),
  refit = REFIT_MODELS
)

m_gc_wind_speed <- fit_cat(
  group_cohesion ~ wind_speed_s + (1 | year),
  data  = hw_gc,
  file  = file.path(model_dir, "m_gc_wind_speed.rds"),
  refit = REFIT_MODELS
)

m_gc_max_wave_height <- fit_cat(
  group_cohesion ~ max_wave_height_s + (1 | year),
  data  = hw_gc,
  file  = file.path(model_dir, "m_gc_max_wave_height.rds"),
  refit = REFIT_MODELS
)

m_gc_current_speed <- fit_cat(
  group_cohesion ~ current_speed_s + (1 | year),
  data  = hw_gc,
  file  = file.path(model_dir, "m_gc_current_speed.rds"),
  refit = REFIT_MODELS
)

m_gc_week <- fit_cat(
  group_cohesion ~ weeknr_s + (1 | year),
  data  = hw_gc,
  file  = file.path(model_dir, "m_gc_week.rds"),
  refit = REFIT_MODELS
)

m_gc_tod <- fit_cat(
  group_cohesion ~ s(ToD, bs = "cc", k = 8) + (1 | year),
  data  = hw_gc,
  knots = knots_gc,
  file  = file.path(model_dir, "m_gc_tod.rds"),
  refit = REFIT_MODELS
)

# Combined hypothesis-driven models
m_gc_env_nodist <- fit_cat(
  group_cohesion ~ depth_s + SST_s + (1 | year),
  data  = hw_gc,
  file  = file.path(model_dir, "m_gc_env_nodist.rds"),
  refit = REFIT_MODELS
)

m_gc_env <- fit_cat(
  group_cohesion ~ depth_s + dist_m_s + SST_s + (1 | year),
  data  = hw_gc,
  file  = file.path(model_dir, "m_gc_env.rds"),
  refit = REFIT_MODELS
)

m_gc_env_nosst <- fit_cat(
  group_cohesion ~ depth_s + dist_m_s + (1 | year),
  data  = hw_gc,
  file  = file.path(model_dir, "m_gc_env_nosst.rds"),
  refit = REFIT_MODELS
)

m_gc_env_big <- fit_cat(
  group_cohesion ~ depth_s + dist_m_s + SST_s +
    wind_speed_s + max_wave_height_s + current_speed_s + (1 | year),
  data  = hw_gc,
  file  = file.path(model_dir, "m_gc_env_big.rds"),
  refit = REFIT_MODELS
)

m_gc_env_time <- fit_cat(
  group_cohesion ~ depth_s + dist_m_s + SST_s + weeknr_s +
    s(ToD, bs = "cc", k = 8) + (1 | year),
  data  = hw_gc,
  knots = knots_gc,
  file  = file.path(model_dir, "m_gc_env_time.rds"),
  refit = REFIT_MODELS
)

mod_list_gc <- list(
  m_gc_env_nodist      = m_gc_env_nodist,
  m_gc_env             = m_gc_env,
  m_gc_env_nosst       = m_gc_env_nosst,
  m_gc_env_time        = m_gc_env_time,
  m_gc_env_big         = m_gc_env_big,

  m_gc_depth           = m_gc_depth,
  m_gc_dist            = m_gc_dist,
  m_gc_sst             = m_gc_sst,
  m_gc_wind_speed      = m_gc_wind_speed,
  m_gc_max_wave_height = m_gc_max_wave_height,
  m_gc_current_speed   = m_gc_current_speed,
  m_gc_week            = m_gc_week,
  m_gc_tod             = m_gc_tod,

  m0_gc                = m0_gc
)

mod_list_gc <- purrr::map(mod_list_gc, ensure_loo)


# ==============================================================================
# 5) PSIS-LOO model selection
# ==============================================================================
loo_tbl <- imap_dfr(mod_list_gc, \(fit, nm) {
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
  ~model,                 ~Model_label,                                     ~Predictors,
  "m_gc_env_nodist",      "Depth + SST (final)",                             "Depth, SST",
  "m_gc_env",             "Depth + SST + distance",                          "Depth, SST, distance to shore",
  "m_gc_env_nosst",       "Depth + distance",                                "Depth, distance to shore",
  "m_gc_env_time",        "Depth + SST + distance + week + time-of-day",     "Depth, SST, distance to shore, week of year, time of day",
  "m_gc_env_big",         "Depth + SST + distance + wind + wave + current",  "Depth, SST, distance to shore, wind speed, max wave height, current speed",

  "m_gc_depth",           "Depth only",                                      "Depth",
  "m_gc_dist",            "Distance only",                                   "Distance to shore",
  "m_gc_sst",             "SST only",                                        "SST",
  "m_gc_wind_speed",      "Wind speed only",                                 "Wind speed",
  "m_gc_max_wave_height", "Max wave height only",                            "Max wave height",
  "m_gc_current_speed",   "Current speed only",                              "Current speed",
  "m_gc_week",            "Week of year only",                               "Week of year",
  "m_gc_tod",             "Time of day only (cyclic)",                       "Time of day",
  "m0_gc",                "Null (year only)",                                "Intercept only"
)

tab_models_readable <- loo_tbl %>%
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
  tab_models_readable,
  booktabs = TRUE,
  caption = paste0(
    "Table Sx. Model selection for group cohesion using PSIS-LOO. ",
    "Models ranked by increasing LOOIC (lower is better). ",
    "All candidate models include a random intercept for year."
  )
)

write_csv(tab_models_readable, file.path(table_dir, "TableSx_groupcohesion_model_selection.csv"))

m_gc_final <- mod_list_gc[["m_gc_env_nodist"]]


# ==============================================================================
# 6) Spatial sensitivity check: 10-km cell random effect
# ==============================================================================
m_gc_final_cell <- fit_cat(
  group_cohesion ~ depth_s + SST_s + (1 | year) + (1 | cell),
  data  = hw_gc,
  file  = file.path(model_dir, "m_gc_env_nodist_cell_10kmRE.rds"),
  refit = REFIT_MODELS
)
m_gc_final_cell <- ensure_loo(m_gc_final_cell)

cat("\nSpatial sensitivity (cell RE) LOO compare:\n")
print(loo_compare(m_gc_final$criteria$loo, m_gc_final_cell$criteria$loo))

# Moran's I on class-wise residuals (population-level predictions)
ep <- posterior_epred(m_gc_final, newdata = hw_gc, re_formula = NA)
p_mean <- apply(ep, c(2, 3), mean)

classes <- dimnames(ep)[[3]]
if (is.null(classes)) classes <- levels(hw_gc$group_cohesion)
colnames(p_mean) <- classes

coords <- as.matrix(hw_gc[, c("x_utm_km", "y_utm_km")])

# Small jitter (1 m = 0.001 km) avoids identical-point warnings
set.seed(1)
coords_j <- coords + matrix(rnorm(length(coords), sd = 0.001), ncol = 2)

knn_k <- 8
nb <- knn2nb(knearneigh(coords_j, k = knn_k))
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_results_gc <- map_dfr(classes, function(cl) {
  obs01 <- as.integer(hw_gc$group_cohesion == cl)
  resid <- obs01 - p_mean[, cl]
  mt <- moran.test(resid, lw, zero.policy = TRUE)
  tibble(
    class   = cl,
    Moran_I = unname(mt$estimate[["Moran I statistic"]]),
    p_value = mt$p.value
  )
})

print(moran_results_gc)
write_csv(moran_results_gc, file.path(table_dir, "TableSx_groupcohesion_moranI.csv"))


# ==============================================================================
# 7) Diagnostics for the final model
# ==============================================================================
loo_final <- loo(m_gc_final)
print(loo_final)

png(file.path(fig_dir, "diag_groupcohesion_loo.png"), width = 1400, height = 900, res = 180)
plot(loo_final)
dev.off()

np <- nuts_params(m_gc_final)
n_div <- sum(np$Parameter == "divergent__" & np$Value == 1)
n_treedepth <- sum(np$Parameter == "treedepth__" & np$Value >= ctrl_sp$max_treedepth)

cat("\nDivergences:", n_div, "\n")
cat("Max treedepth hits:", n_treedepth, "\n")

diag_tbl_gc <- tibble(
  metric = c("divergences", "max_treedepth_hits"),
  value  = c(n_div, n_treedepth)
)
write_csv(diag_tbl_gc, file.path(table_dir, "groupcohesion_sampler_diagnostics.csv"))

if (!is.null(loo_final$diagnostics$pareto_k)) {
  pareto_k_tbl_gc <- tibble(
    observation = seq_along(loo_final$diagnostics$pareto_k),
    pareto_k    = loo_final$diagnostics$pareto_k
  )
  write_csv(pareto_k_tbl_gc, file.path(table_dir, "groupcohesion_pareto_k.csv"))
}

pp_gc <- pp_check(m_gc_final, type = "bars", ndraws = 100)
ggsave(file.path(fig_dir, "diag_groupcohesion_ppcheck_bars.png"), pp_gc,
       width = 7.5, height = 4.8, dpi = 300)


# ==============================================================================
# 8) Predicted probabilities along unscaled gradients
# ==============================================================================
gc_pal <- c(
  "AA"  = "#4F5D75",
  "FC"  = "#EF8354",
  "FCE" = "#BFC0C0"
)

plot_gc_prob_curve <- function(
    fit,
    data,
    var = c("depth", "SST"),
    n = 200,
    ci = 0.95,
    re_formula = NA,
    trim = c(0.05, 0.95),
    xlab = NULL,
    title = NULL,
    ndraws = NULL
) {
  var <- match.arg(var)

  x_raw <- data[[var]]
  x_seq <- seq(
    quantile(x_raw, trim[1], na.rm = TRUE),
    quantile(x_raw, trim[2], na.rm = TRUE),
    length.out = n
  )

  mu_depth <- mean(data$depth, na.rm = TRUE)
  sd_depth <- sd(data$depth, na.rm = TRUE)
  mu_SST   <- mean(data$SST,   na.rm = TRUE)
  sd_SST   <- sd(data$SST,     na.rm = TRUE)

  nd <- tibble(
    depth = mu_depth,
    SST   = mu_SST,
    year  = factor(levels(data$year)[1], levels = levels(data$year))
  ) %>%
    slice(rep(1, n))

  nd[[var]] <- x_seq

  nd <- nd %>%
    mutate(
      depth_s = (depth - mu_depth) / sd_depth,
      SST_s   = (SST   - mu_SST)   / sd_SST
    )

  ep <- posterior_epred(fit, newdata = nd, re_formula = re_formula, ndraws = ndraws)

  alpha <- (1 - ci) / 2
  probs <- c(alpha, 0.5, 1 - alpha)
  summ  <- apply(ep, c(2, 3), \(v) quantile(v, probs = probs))

  class_names <- dimnames(ep)[[3]]
  if (is.null(class_names)) class_names <- levels(data$group_cohesion)

  out <- bind_rows(lapply(seq_along(class_names), function(k) {
    tibble(
      x     = x_seq,
      class = factor(class_names[k], levels = class_names),
      lo    = summ[1, , k],
      mid   = summ[2, , k],
      hi    = summ[3, , k]
    )
  }))

  gc_lab <- c(
    "AA"   = "Adults only (AA)",
    "FC"   = "Female–calf (FC)",
    "FCE"  = "Female–calf escorted (FCE)",
    "rare" = "Rare / other"
  )
  labs_map <- setNames(levels(out$class), levels(out$class))
  labs_map[names(gc_lab)] <- gc_lab

  if (is.null(xlab)) xlab <- ifelse(var == "depth", "Depth (m)", "Sea surface temperature (°C)")
  if (is.null(title)) title <- paste0("Group cohesion along a ", ifelse(var == "depth", "depth", "SST"), " gradient")

  ggplot(out, aes(x = x, y = mid, ymin = lo, ymax = hi, colour = class, fill = class)) +
    geom_ribbon(alpha = 0.25, linewidth = 0) +
    geom_line(linewidth = 1.1) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_colour_manual(
      values = gc_pal,
      breaks = levels(out$class),
      labels = labs_map[levels(out$class)],
      drop   = FALSE
    ) +
    scale_fill_manual(
      values = gc_pal,
      breaks = levels(out$class),
      labels = labs_map[levels(out$class)],
      drop   = FALSE
    ) +
    labs(
      x = xlab,
      y = "Predicted probability",
      title = title,
      subtitle = paste0("Median and ", ci * 100, "% credible intervals (population-level)"),
      colour = "Group cohesion",
      fill   = "Group cohesion"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

p_gc_depth <- plot_gc_prob_curve(m_gc_final, hw_gc, var = "depth", xlab = "Depth (m)")
p_gc_sst   <- plot_gc_prob_curve(m_gc_final, hw_gc, var = "SST",   xlab = "Sea surface temperature (°C)")

ggsave(file.path(fig_dir, "fig_gc_depth_prob.png"), p_gc_depth, width = 7.5, height = 4.8, dpi = 400)
ggsave(file.path(fig_dir, "fig_gc_sst_prob.png"),   p_gc_sst,   width = 7.5, height = 4.8, dpi = 400)

p_gc_both <- (p_gc_depth / p_gc_sst) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "fig_gc_gradients_combined.png"), p_gc_both,
       width = 8.5, height = 9.5, dpi = 400)


# ==============================================================================
# 9) Supplement table: odds ratios on original units
# ==============================================================================
draws_final <- as_draws_df(m_gc_final)

get_beta <- function(name) {
  if (!name %in% names(draws_final)) stop("Missing coefficient: ", name)
  draws_final[[name]]
}

summ_or <- function(beta_draws, digits = 2) {
  OR <- exp(beta_draws)
  tibble(
    OR_median = round(median(OR), digits),
    lo95      = round(unname(quantile(OR, 0.025)), digits),
    hi95      = round(unname(quantile(OR, 0.975)), digits)
  )
}

gc_levels <- levels(m_gc_final$data$group_cohesion)
gc_ref    <- gc_levels[1]

gc_lab <- c(
  "AA"   = "Adults only (AA)",
  "FC"   = "Female–calf (FC)",
  "FCE"  = "Female–calf escorted (FCE)",
  "rare" = "Rare / other"
)
ref_name_gc <- ifelse(gc_ref %in% names(gc_lab), gc_lab[[gc_ref]], gc_ref)

sd_depth <- gc_meta$sd_depth
sd_SST   <- gc_meta$sd_SST

depth_step_m <- 10
sst_step_C   <- 1

pred_spec_gc <- tribble(
  ~Predictor, ~coef_suffix, ~Unit,                              ~mult,
  "Depth",    "depth_s",    paste0("per ", depth_step_m, " m"),  (depth_step_m / sd_depth),
  "SST",      "SST_s",      paste0("per ", sst_step_C, " °C"),   (sst_step_C   / sd_SST)
)

nonref_levels <- gc_levels[-1]

or_tbl_gc <- map_dfr(nonref_levels, function(st) {
  st_name  <- ifelse(st %in% names(gc_lab), gc_lab[[st]], st)
  contrast <- paste0(st_name, " vs ", ref_name_gc)
  prefix   <- paste0("b_mu", st, "_")

  map_dfr(seq_len(nrow(pred_spec_gc)), function(i) {
    p <- pred_spec_gc[i, ]
    beta_unit <- get_beta(paste0(prefix, p$coef_suffix)) * p$mult
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

caption_or_gc <- paste0(
  "Table Sx. Odds ratios from the final group-cohesion model (brms; family = categorical; N = ", gc_meta$N,
  "; year random intercept with ", gc_meta$n_year, " levels). ", ref_name_gc, " was the reference category. ",
  "Odds ratios are expressed on original units (per ", depth_step_m, " m depth; per ", sst_step_C, " °C SST). ",
  "Depth and SST were standardized in the fitted model and back-transformed using SDs from the analysis dataset ",
  "(SD_depth = ", round(sd_depth, 1), " m; SD_SST = ", round(sd_SST, 2), " °C)."
)

kable(or_tbl_gc, booktabs = TRUE, caption = caption_or_gc)
write_csv(or_tbl_gc, file.path(table_dir, "TableSx_groupcohesion_OR.csv"))


# ==============================================================================
# 10) Supplement table: predicted probabilities under simple perturbations
# ==============================================================================
mu_depth <- gc_meta$mean_depth
mu_SST   <- gc_meta$mean_SST
ref_year <- factor(levels(hw_gc$year)[1], levels = levels(hw_gc$year))

nd_gc <- tibble(
  Scenario = c("Baseline (means)", paste0("Depth +", depth_step_m, " m"), paste0("SST +", sst_step_C, " °C")),
  depth = c(mu_depth, mu_depth + depth_step_m, mu_depth),
  SST   = c(mu_SST,   mu_SST,                  mu_SST + sst_step_C),
  year  = ref_year
) %>%
  mutate(
    depth_s = (depth - mu_depth) / sd_depth,
    SST_s   = (SST   - mu_SST)   / sd_SST
  )

ep_gc <- posterior_epred(m_gc_final, newdata = nd_gc, re_formula = NA)
classes <- dimnames(ep_gc)[[3]]
if (is.null(classes)) classes <- levels(hw_gc$group_cohesion)

summ_ci <- function(x, digits = 3) {
  tibble(
    lo95 = round(unname(quantile(x, 0.025)), digits),
    med  = round(unname(quantile(x, 0.50)),  digits),
    hi95 = round(unname(quantile(x, 0.975)), digits)
  )
}

abs_tbl_gc <- bind_rows(lapply(seq_len(nrow(nd_gc)), function(i) {
  tmp <- bind_rows(lapply(seq_along(classes), function(k) {
    s <- summ_ci(ep_gc[, i, k], digits = 3)
    tibble(State = classes[k], lo95 = s$lo95, med = s$med, hi95 = s$hi95)
  }))
  tmp$Scenario <- nd_gc$Scenario[i]
  tmp
})) %>%
  select(Scenario, State, med, lo95, hi95)

abs_tbl_gc_wide <- abs_tbl_gc %>%
  mutate(val = paste0(med, " [", lo95, ", ", hi95, "]")) %>%
  select(Scenario, State, val) %>%
  pivot_wider(names_from = State, values_from = val)

caption_prob_gc <- paste0(
  "Table Sy. Predicted group-cohesion probabilities (median [95% CrI]) at typical conditions and under simple perturbations ",
  "(population-level predictions; year random effects excluded). Baseline uses predictor means from the analysis dataset ",
  "(depth = ", round(mu_depth, 1), " m; SST = ", round(mu_SST, 2), " °C). Perturbations change one predictor at a time ",
  "(depth +", depth_step_m, " m; SST +", sst_step_C, " °C), holding the other predictor at baseline."
)

kable(abs_tbl_gc_wide, booktabs = TRUE, caption = caption_prob_gc)
write_csv(abs_tbl_gc_wide, file.path(table_dir, "TableSy_groupcohesion_probabilities.csv"))

message("DONE. Outputs written to: ", fig_dir, " and ", table_dir, ". Models cached in: ", model_dir)

save.image(file = WS_FILE)
cat("Workspace path:", normalizePath(WS_FILE, winslash = "/", mustWork = FALSE), "\n")
