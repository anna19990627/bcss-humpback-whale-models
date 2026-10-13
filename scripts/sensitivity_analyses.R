################################################################################
# REVIEWER-REQUESTED SENSITIVITY ANALYSES — BAZARUTO HUMPBACK WHALES
#
# Run from the project root AFTER the two primary model scripts have been run.
# This file runs the three additional sensitivity analyses:
#   S2: Within-season adjustment of the SST-occurrence response
#   S3: Leave-one-taxonomic-group-out pseudoabsence analysis (SST and depth)
#   S4: Behavioural state ~ group type versus group size (Tables S14/S15)
#
# Required existing project files:
#   models/sdm_env/obs_sdm_prepped.rds
#   models/sdm_env/m_depth_sst_ar1gap.rds
#   clean_dfs/hw.rds
# Outputs are saved under figs/, tables/, and models/reviewer_sensitivity/.
################################################################################

library(tidyverse)
library(lubridate)
library(sf)
library(sdmTMB)
library(sdmTMBextra)
library(fmesher)
library(geodata)
library(brms)
library(loo)
library(scales)
library(patchwork)

set.seed(123)
REFIT <- FALSE  # TRUE forces all seven sensitivity models to be refitted

out_fig <- "figs/sensitivity"
out_tab <- "tables/sensitivity"
out_mod <- "models/sensitivity"
for (p in c(out_fig, out_tab, out_mod)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing project file: ", path,
                               "\nRun the corresponding primary script first.")
  readRDS(path)
}

fit_cached <- function(path, fit_function) {
  if (file.exists(path) && !REFIT) {
    message("Loading: ", path)
    return(readRDS(path))
  }
  message("Fitting: ", basename(path))
  model <- fit_function()
  saveRDS(model, path)
  model
}

# sdmTMB's predict(..., nsim=) has returned different object shapes by version.
simulation_matrix <- function(pred, n_points) {
  sims <- if (is.matrix(pred)) pred else attr(pred, "sims")
  if (is.null(sims) && is.list(pred)) sims <- pred$sims
  if (is.null(sims) && is.data.frame(pred)) {
    cols <- grep("^sim", names(pred), value = TRUE)
    if (length(cols)) sims <- as.matrix(pred[, cols, drop = FALSE])
  }
  if (is.null(sims)) stop("No simulation matrix found in sdmTMB predictions.")
  if (nrow(sims) != n_points && ncol(sims) == n_points) sims <- t(sims)
  if (nrow(sims) != n_points) stop("Unexpected prediction matrix dimensions.")
  sims
}

# ============================================================================
# 1) Occurrence data and original model
# ============================================================================
occ_all <- read_required("models/sdm_env/obs_sdm_prepped.rds")
occ_original <- read_required("models/sdm_env/m_depth_sst_ar1gap.rds")

origin <- as.Date(min(occ_all$Date_Time, na.rm = TRUE))
occ <- occ_all %>%
  mutate(week_global = as.integer(floor(as.numeric(as.Date(Date_Time) - origin) / 7) + 1L)) %>%
  filter(!is.na(week_global), !is.na(x), !is.na(y), !is.na(observed),
         !is.na(depth_scaled), !is.na(SST_scaled), !is.na(year))
all_weeks <- seq.int(min(occ$week_global), max(occ$week_global))

mu_sst <- mean(occ$SST)
sd_sst <- sd(occ$SST)
mu_depth <- mean(occ$depth)
sd_depth <- sd(occ$depth)
ref_year <- factor(levels(occ$year)[1], levels = levels(occ$year))
ref_week <- as.integer(round(median(occ$week_global)))
ref_x <- mean(occ$x)
ref_y <- mean(occ$y)

sst_grid <- seq(quantile(occ$SST, 0.05), quantile(occ$SST, 0.95), length.out = 200)
depth_grid <- seq(quantile(occ$depth, 0.05), quantile(occ$depth, 0.95), length.out = 200)

new_gradient <- function(variable = c("SST", "depth")) {
  variable <- match.arg(variable)
  sst <- if (variable == "SST") sst_grid else rep(mu_sst, length(depth_grid))
  depth <- if (variable == "depth") depth_grid else rep(mu_depth, length(sst_grid))
  tibble(x = ref_x, y = ref_y, year = ref_year, week_global = ref_week,
         week_in_season_c = 0, SST = sst, depth = depth,
         SST_scaled = (sst - mu_sst) / sd_sst,
         depth_scaled = (depth - mu_depth) / sd_depth)
}

predict_gradient <- function(model, grid, variable) {
  set.seed(123)
  pr <- predict(model, newdata = grid, type = "link", nsim = 500,
                re_form = NA, re_form_iid = NA)
  prob <- plogis(simulation_matrix(pr, nrow(grid)))
  tibble(value = grid[[variable]],
         mid = apply(prob, 1, median),
         lo = apply(prob, 1, quantile, probs = 0.025),
         hi = apply(prob, 1, quantile, probs = 0.975))
}

# Reconstruct the same 2-km triangulation/coastline as the primary script.
# A reduced pseudoabsence dataset needs a new observation-to-mesh mapping.
mesh_cache <- new.env(parent = emptyenv())
mesh_for <- function(data = NULL) {
  if (is.null(mesh_cache$triangulation)) {
    mesh_cache$triangulation <- sdmTMB::make_mesh(occ_all, c("x", "y"), cutoff = 2)
    coast <- sf::st_as_sf(geodata::gadm("Mozambique", level = 1, path = "."))
    bbox <- sf::st_bbox(c(xmin = 35, xmax = 36, ymin = -22.5, ymax = -21),
                        crs = sf::st_crs(4326))
    coast <- sf::st_crop(coast, bbox) %>% sf::st_transform(32736)
    mesh_cache$coast <- sf::st_transform(
      coast, crs = fmesher::fm_crs_set_lengthunit(
        fmesher::fm_CRS(SRS_string = "EPSG:32736"), "km"))
  }
  m <- if (is.null(data)) mesh_cache$triangulation else
    sdmTMB::make_mesh(data, xy_cols = c("x", "y"),
                      mesh = mesh_cache$triangulation$mesh)
  sdmTMBextra::add_barrier_mesh(m, mesh_cache$coast, range_fraction = 0.1,
                                proj_scaling = 1000, plot = FALSE)
}

# ============================================================================
# 2) Within-season timing sensitivity — Figure S2
# ============================================================================
occ_season <- occ %>%
  mutate(date = as.Date(Date_Time),
         season_start = as.Date(paste0(year(date), "-05-01")),
         week_in_season = as.numeric(date - season_start) / 7,
         week_in_season_c = week_in_season - mean(week_in_season, na.rm = TRUE))

m_season <- fit_cached(file.path(out_mod, "occurrence_seasonal.rds"), function() {
  sdmTMB(observed ~ s(depth_scaled, bs = "ts") + s(SST_scaled, bs = "ts") +
           s(week_in_season_c, bs = "ts") + (1 | year),
         data = occ_season, mesh = mesh_for(), time = "week_global",
         extra_time = all_weeks, family = binomial(link = "logit"),
         spatial = "on", spatiotemporal = "ar1", reml = TRUE)
})

sst_nd <- new_gradient("SST")
season_curves <- bind_rows(
  predict_gradient(occ_original, sst_nd, "SST") %>% mutate(Model = "Original model"),
  predict_gradient(m_season, sst_nd, "SST") %>% mutate(Model = "Season-adjusted model")
)

season_compare <- season_curves %>%
  select(value, mid, Model) %>%
  pivot_wider(names_from = Model, values_from = mid) %>%
  mutate(abs_difference = abs(`Original model` - `Season-adjusted model`))
season_results <- tibble(AIC_original = as.numeric(AIC(occ_original)),
                         AIC_season = as.numeric(AIC(m_season)),
                         delta_AIC = AIC_season - AIC_original,
                         mean_abs_difference = mean(season_compare$abs_difference),
                         max_abs_difference = max(season_compare$abs_difference))
write_csv(season_results, file.path(out_tab, "S2_seasonality_results.csv"))
print(season_results)

p_s2 <- ggplot(season_curves, aes(value, mid, group = Model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = Model), alpha = 0.12, colour = NA) +
  geom_line(aes(linetype = Model), linewidth = 1) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(x = "Sea-surface temperature (°C)", y = "Predicted relative occurrence",
       linetype = NULL, fill = NULL) + theme_classic()
ggsave(file.path(out_fig, "Figure_S2_seasonality.png"), p_s2,
       width = 7.5, height = 4.8, dpi = 300)

# ============================================================================
# 3) Target-group pseudoabsence sensitivity — Figure S3
# ============================================================================
pa <- occ %>%
  mutate(species_lower = str_to_lower(Species),
         pseudo_group = case_when(
           observed == 1 ~ "Humpback whale",
           str_detect(species_lower, "dolphin|pilot whale|killer whale|minke whale") ~ "Other cetaceans",
           str_detect(species_lower, "turtle") ~ "Turtles",
           str_detect(species_lower, "shark|hammerhead|guitarfish|ray|manta|whipray|stingray") ~ "Sharks and rays",
           TRUE ~ "Other marine taxa"))

pa_counts <- pa %>% filter(observed == 0) %>% count(pseudo_group)
write_csv(pa_counts, file.path(out_tab, "S3_pseudoabsence_counts.csv"))
print(pa_counts)

excluded <- c("Other cetaceans", "Turtles", "Sharks and rays")
pa_models <- list("Original model" = occ_original)
for (taxon in excluded) {
  label <- paste("No", tolower(taxon))
  dat <- pa %>% filter(observed == 1 | pseudo_group != taxon)
  path <- file.path(out_mod, paste0("occurrence_without_", gsub(" ", "_", tolower(taxon)), ".rds"))
  pa_models[[label]] <- fit_cached(path, local({
    d <- dat
    function() sdmTMB(observed ~ s(depth_scaled, bs = "ts") +
                        s(SST_scaled, bs = "ts") + (1 | year),
                      data = d, mesh = mesh_for(d), time = "week_global",
                      extra_time = all_weeks, family = binomial(link = "logit"),
                      spatial = "on", spatiotemporal = "ar1", reml = TRUE)
  }))
}

pa_curves <- map_dfr(names(pa_models), function(label) {
  fit <- pa_models[[label]]
  bind_rows(predict_gradient(fit, sst_nd, "SST") %>% mutate(gradient = "SST"),
            predict_gradient(fit, new_gradient("depth"), "depth") %>% mutate(gradient = "depth")) %>%
    mutate(Model = label)
})

curve_metrics <- function(d) {
  d <- arrange(d, value)
  peak <- which.max(d$mid)
  band <- d$value[d$mid >= 0.95 * max(d$mid)]
  after_peak <- seq_len(nrow(d)) >= peak
  half <- which(d$mid <= 0.5 * d$mid[peak] & after_peak)
  tibble(peak = d$value[peak], band_low = min(band), band_high = max(band),
         half_max_after_peak = if (length(half)) d$value[half[1]] else NA_real_)
}

pa_summary <- pa_curves %>%
  group_by(Model, gradient) %>%
  group_modify(~ {
    original <- pa_curves %>%
      filter(Model == "Original model", gradient == .y$gradient) %>% arrange(value)
    current <- arrange(.x, value)
    bind_cols(curve_metrics(current),
              tibble(mean_abs_difference = mean(abs(current$mid - original$mid)),
                     max_abs_difference = max(abs(current$mid - original$mid))))
  }) %>% ungroup()
write_csv(pa_summary, file.path(out_tab, "S3_pseudoabsence_results.csv"))
print(pa_summary)

plot_pa <- function(gradient, xlab) {
  d <- filter(pa_curves, .data$gradient == .env$gradient)
  ggplot(d, aes(value, mid, group = Model)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Model), alpha = 0.08, colour = NA) +
    geom_line(aes(linetype = Model), linewidth = 0.9) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(x = xlab, y = "Predicted relative occurrence", linetype = NULL, fill = NULL) +
    theme_classic()
}
p_s3 <- (plot_pa("SST", "(A) Sea-surface temperature (°C)") /
         plot_pa("depth", "(B) Water depth (m)")) + plot_layout(guides = "collect")
ggsave(file.path(out_fig, "Figure_S3_pseudoabsence.png"), p_s3,
       width = 7.5, height = 9, dpi = 300)

# ============================================================================
# 4) Group-type sensitivity in the behavioural model — Figure S4, Tables S14/S15
# ============================================================================
hw <- read_required("clean_dfs/hw.rds") %>%
  rename_with(~ "current_speed", any_of("currend_speed")) %>%
  rename_with(~ "max_wave_height", any_of("max_wave_hight")) %>%
  sf::st_drop_geometry() %>%
  mutate(b_state = factor(b_state), group_cohesion = factor(group_cohesion),
         year = factor(year))

# Exactly the main behavioural-model complete cases (N = 376) BEFORE filtering
# for group type. This retains the ORIGINAL SST/wave standardization.
behav_main <- hw %>%
  filter(!is.na(b_state), !is.na(year), !is.na(Date_Time), !is.na(QTY),
         !is.na(depth), !is.na(dist_m), !is.na(SST), !is.na(weeknr),
         !is.na(wind_speed), !is.na(max_wave_height), !is.na(current_speed)) %>%
  mutate(SST_s = as.numeric(scale(SST)),
         max_wave_height_s = as.numeric(scale(max_wave_height)))

behav_common <- behav_main %>%
  filter(!is.na(group_cohesion)) %>%
  mutate(group_cohesion = factor(group_cohesion, levels = c("AA", "FC", "FCE"))) %>%
  droplevels()
message("Behaviour: main N = ", nrow(behav_main), "; common N = ", nrow(behav_common))
print(table(behav_common$group_cohesion))
if (nrow(behav_main) != 376 || nrow(behav_common) != 332) {
  warning("Behaviour sample differs from manuscript (376/332); review your data.")
}

fit_behav <- function(formula, filename) {
  fit_cached(file.path(out_mod, filename), function() {
    brms::brm(formula, data = behav_common, family = brms::categorical(),
              chains = 4, cores = 4, iter = 2000, warmup = 1000,
              control = list(adapt_delta = 0.999), seed = 123,
              save_pars = brms::save_pars(all = TRUE))
  })
}

m_base <- fit_behav(b_state ~ SST_s + max_wave_height_s + (1 | year),
                    "behaviour_environment.rds")
m_qty <- fit_behav(b_state ~ SST_s + QTY + max_wave_height_s + (1 | year),
                   "behaviour_group_size.rds")
m_type <- fit_behav(b_state ~ SST_s + group_cohesion + max_wave_height_s + (1 | year),
                    "behaviour_group_type.rds")

loo_models <- list(Base = loo::loo(m_base), QTY = loo::loo(m_qty),
                   Group_type = loo::loo(m_type))
loo_comparison <- loo::loo_compare(loo_models)
tab_s14 <- as.data.frame(loo_comparison) %>%
  rownames_to_column("Model") %>%
  transmute(Model, delta_ELPD = elpd_diff, SE = se_diff)
write_csv(tab_s14, file.path(out_tab, "Table_S14_group_type_LOO.csv"))
print(tab_s14)

# Eight fixed-effect terms shown in Table S15 (excluding intercepts).
coef_names <- expand_grid(state = c("S", "T"),
                          predictor = c("SST_s", "group_cohesionFC",
                                        "group_cohesionFCE", "max_wave_height_s")) %>%
  mutate(term = paste0("mu", state, "_", predictor),
         Contrast = if_else(state == "S", "Social vs Resting", "Travelling vs Resting"),
         Predictor = recode(predictor, SST_s = "SST",
                            group_cohesionFC = "FC vs AA",
                            group_cohesionFCE = "FCE vs AA",
                            max_wave_height_s = "Wave"))
fixed_terms <- as.data.frame(brms::fixef(m_type)) %>% rownames_to_column("term")
if (!all(coef_names$term %in% fixed_terms$term)) {
  stop("Unexpected brms coefficient names. Found: ", paste(fixed_terms$term, collapse = ", "))
}
tab_s15 <- coef_names %>%
  left_join(fixed_terms, by = "term") %>%
  transmute(Contrast, Predictor, Estimate,
            SD = Est.Error, CrI_low = Q2.5, CrI_high = Q97.5)
write_csv(tab_s15, file.path(out_tab, "Table_S15_group_type_coefficients.csv"))
print(tab_s15)

nd <- tibble(group_cohesion = factor(c("AA", "FC", "FCE"),
                                     levels = c("AA", "FC", "FCE")),
             SST_s = 0, max_wave_height_s = 0,
             year = factor(levels(behav_common$year)[1], levels = levels(behav_common$year)))
posterior_prob <- brms::posterior_epred(m_type, newdata = nd, re_formula = NA)
states <- dimnames(posterior_prob)[[3]]
if (is.null(states)) states <- levels(behav_common$b_state)
probs <- map_dfr(seq_len(nrow(nd)), function(i) {
  map_dfr(seq_along(states), function(k) {
    z <- posterior_prob[, i, k]
    tibble(Group = as.character(nd$group_cohesion[i]), State = states[k],
           median = median(z), lo = unname(quantile(z, 0.025)),
           hi = unname(quantile(z, 0.975)))
  })
}) %>%
  mutate(Group = factor(Group, levels = c("AA", "FC", "FCE")),
         State = factor(State, levels = c("R", "S", "T")))
write_csv(probs, file.path(out_tab, "S4_group_type_probabilities.csv"))

palette <- c(R = "#4F5D75", S = "#EF8354", T = "#BFC0C0")
p_s4 <- ggplot(probs, aes(Group, median, colour = State)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.08,
                position = position_dodge(width = 0.3)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  scale_colour_manual(values = palette,
                      labels = c(R = "Resting (R)", S = "Social (S)", T = "Travelling (T)")) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(x = "Group type", y = "Predicted probability", colour = "Behavioural state") +
  theme_classic()
ggsave(file.path(out_fig, "Figure_S4_group_type.png"), p_s4,
       width = 8, height = 5.2, dpi = 400)

message("Done. Figures: ", out_fig, "; tables: ", out_tab, "; fitted models: ", out_mod)
