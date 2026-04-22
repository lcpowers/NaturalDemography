## ipm_functions.R
## Site * aspect * year vital rate models (microtopo-first then climate).
## Source this file at the top of IPM analysis scripts.
## Assumes working directory is the project root (natural_demog/).

library(tidyverse)
library(lme4)

##### 1. Data loading #####
sites <- c("guanella","niwot","jones")
aspects <- c("north","south","top")
years <- c(2023:2025)

mod_data <- readRDS("data/demog/vr_modData/mod_data.rds")
demog_data <- readRDS("data/demog/vr_modData/mod_demog_data.rds")

clim <- read_csv("data/climate/for_analyses/clim_summary.csv", show_col_types = FALSE)
clim_vars <- names(clim)[!names(clim) %in% c("site","pop","year")]

microtopo <- read_csv("data/output/model_microtopo_25cm.csv", show_col_types = FALSE)
mt_vars <- names(microtopo)[!names(microtopo) %in% c("site","pop","plantID","elv_pt")]

scale_params <- read_csv("data/demog/vr_modData/scale_params.csv", show_col_types = FALSE)

# Germination rates by site * pop; replace zeros; derive aspect.cat from pop name
germ_rates_pop <- demog_data$germ_rates %>%
  mutate(germ_rate = ifelse(germ_rate == 0,
                            min(germ_rate[germ_rate > 0], na.rm = TRUE),
                            germ_rate),
         aspect.cat = case_when(
           str_starts(pop, "North") ~ "north",
           str_starts(pop, "South") ~ "south",
           str_starts(pop, "Top") ~ "top"
         ))

# Climate means at site * pop * year level
clim_means_pop <- mod_data$growth_mod_data %>%
  group_by(site, pop, aspect.cat, year) %>%
  summarise(across(all_of(clim_vars), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(year = as.character(year))

# Microtopo means at site * pop level
mt_means_pop <- mod_data$growth_mod_data %>%
  mutate(aspect.cat = case_when(
    str_starts(pop, "North") ~ "north",
    str_starts(pop, "South") ~ "south",
    str_starts(pop, "Top") ~ "top"
  )) %>%
  group_by(site, pop, aspect.cat) %>%
  summarise(across(all_of(mt_vars), ~mean(.x, na.rm = TRUE)), .groups = "drop")
##### End 1 #####

##### 2. Load models #####
models <- readRDS("r/modelFitResults/best_mt2clim_mods_2026-04-12.RDS")
list2env(models, envir = .GlobalEnv)

# clim_only_models <- readRDS(tail(sort(list.files("r/modelFitResults",
#                                                    pattern = "best_climOnly_mods.*\\.RDS",
#                                                    full.names = TRUE)), 1))
# list2env(clim_only_models, envir = .GlobalEnv)

clim_models <- readRDS("r/modelFitResults/best_clim_mods_2026-04-15.RDS")
list2env(clim_models, envir = .GlobalEnv)
##### End 2 #####

##### 3. IPM bin setup #####
log_sizes_all <- c(mod_data$growth_mod_data$log_size_t0,
                   mod_data$growth_mod_data$log_size_t0 + mod_data$growth_mod_data$growth)

minsize <- min(log_sizes_all, na.rm = TRUE)
maxsize <- max(log_sizes_all, na.rm = TRUE) + 0.1

n.bin <- 100
vec.bin <- seq(minsize, maxsize, length.out = n.bin + 1)
h <- vec.bin[2] - vec.bin[1]
binmids <- vec.bin[1:n.bin] + h / 2

recruit.mean.log <- demog_data$recruit_log_size
recruit.sd.log <- 0.25
sdlgszcdf <- pnorm(q = vec.bin, mean = recruit.mean.log, sd = recruit.sd.log)
sdlgszprobs <- diff(sdlgszcdf)
sdlgszprobs <- sdlgszprobs / sum(sdlgszprobs)

nutlets_per_flower <- 1.2
##### End 3 #####

##### 4. Kernel functions #####

## predict_lambda(df)
## df must have columns: site, pop, aspect.cat, year.
## Additional raw predictor columns override pop-level means.
## Returns df with lambda column appended.
predict_lambda <- function(df) {
  log_size_center <- scale_params$center[scale_params$var == "log_size_t0"]
  log_size_scale <- scale_params$scale[scale_params$var == "log_size_t0"]

  df$lambda <- NA_real_

  for (i in seq_len(nrow(df))) {
    row <- df[i,]
    site_i <- row$site; pop_i <- row$pop; aspect_i <- row$aspect.cat; year_i <- as.character(row$year)

    clim_row <- clim_means_pop %>% filter(site == site_i, pop == pop_i, year == year_i)
    mt_row <- mt_means_pop %>% filter(site == site_i, pop == pop_i)

    pred_data <- data.frame(log_size_t0 = (binmids - log_size_center) / log_size_scale,
                            site = site_i, aspect.cat = aspect_i,
                            year = factor(year_i, levels = as.character(years)))
    for (v in clim_vars) {
      raw_v <- if (v %in% names(row) && !is.na(row[[v]])) row[[v]] else clim_row[[v]]
      pred_data[[v]] <- (raw_v - scale_params$center[scale_params$var == v]) /
                         scale_params$scale[scale_params$var == v]
    }
    for (v in mt_vars) {
      raw_v <- if (v %in% names(row) && !is.na(row[[v]])) row[[v]] else mt_row[[v]]
      pred_data[[v]] <- (raw_v - scale_params$center[scale_params$var == v]) /
                         scale_params$scale[scale_params$var == v]
    }

    surv_prob <- predict(surv_mt2clim_best_mod, newdata = pred_data, type = "response", re.form = NA)
    growth_mean <- binmids + predict(growth_mt2clim_best_mod, newdata = pred_data, re.form = NA)
    growth_var <- pmax(predict(growth_var_mt2clim_best_mod, newdata = pred_data, re.form = NA), 0)
    flower_prob <- predict(repro_mt2clim_best_mod, newdata = pred_data, type = "response", re.form = NA)
    flower_count <- predict(nflrs_mt2clim_best_mod, newdata = pred_data, type = "response", re.form = NA)

    pop_germ <- germ_rates_pop$germ_rate[germ_rates_pop$site == site_i & germ_rates_pop$pop == pop_i]

    growth_mx <- matrix(0, n.bin, n.bin)
    for (ss in 1:n.bin) {
      growcdf <- pnorm(vec.bin, growth_mean[ss], sqrt(growth_var[ss]))
      grows <- diff(growcdf)
      if (sum(grows) > 0) growth_mx[, ss] <- grows / sum(grows)
    }

    P_kernel <- t(t(growth_mx) * surv_prob)
    F_kernel <- outer(sdlgszprobs, flower_prob * flower_count * nutlets_per_flower * pop_germ)
    df$lambda[i] <- Re(eigen(P_kernel + F_kernel, only.values = TRUE)$values[1])
  }
  df
}

## predict_indiv_lambda(df)
## Expects scaled predictor values (e.g. from growth_mod_data_scaled).
## Requires site, pop, aspect.cat, year, log_size_t0, site.pop, all clim_vars, all mt_vars.
## Returns df with lambda and per-plant vital rate predictions appended.
predict_indiv_lambda <- function(df,use.re=FALSE) {
  log_size_center <- scale_params$center[scale_params$var == "log_size_t0"]
  log_size_scale  <- scale_params$scale[scale_params$var == "log_size_t0"]
  scaled_binmids  <- (binmids - log_size_center) / log_size_scale

  df$lambda <- NA_real_
  df$surv_pred <- NA_real_
  df$growth_pred <- NA_real_
  df$growth_var_pred <- NA_real_
  df$flower_prob_pred <- NA_real_
  df$flower_count_pred <- NA_real_

  re_form <- if(use.re) {~(1|site.pop)} else NA

  for (i in seq_len(nrow(df))) {
    row <- df[i,]
    site_i <- row$site; pop_i <- row$pop; aspect_i <- row$aspect.cat; year_i <- as.character(row$year)

    pred_data <- data.frame(log_size_t0 = scaled_binmids,
                            site = site_i, aspect.cat = aspect_i,
                            year = factor(year_i, levels = as.character(years)),
                            site.pop = paste(site_i, pop_i, sep = "."))
    for (v in c(clim_vars, mt_vars)) {
      pred_data[[v]] <- row[[v]]
    }

    surv_prob <- predict(surv_mt2clim_best_mod, newdata = pred_data, type = "response", re.form = re_form)
    growth_mean <- binmids + predict(growth_mt2clim_best_mod, newdata = pred_data, re.form = re_form)
    growth_var <- pmax(predict(growth_var_mt2clim_best_mod, newdata = pred_data, re.form = re_form), 0)
    flower_prob <- predict(repro_mt2clim_best_mod, newdata = pred_data, type = "response", re.form = re_form)
    flower_count <- predict(nflrs_mt2clim_best_mod, newdata = pred_data, type = "response", re.form = re_form)

    pop_germ <- germ_rates_pop$germ_rate[germ_rates_pop$pop == pop_i & germ_rates_pop$site == site_i]

    growth_mx <- matrix(0, n.bin, n.bin)
    for (ss in 1:n.bin) {
      growcdf <- pnorm(vec.bin, growth_mean[ss], sqrt(growth_var[ss]))
      grows <- diff(growcdf)
      if (sum(grows) > 0) growth_mx[, ss] <- grows / sum(grows)
    }

    P_kernel <- t(t(growth_mx) * surv_prob)
    F_kernel <- outer(sdlgszprobs, flower_prob * flower_count * nutlets_per_flower * pop_germ)
    df$lambda[i] <- Re(eigen(P_kernel + F_kernel, only.values = TRUE)$values[1])

    pred_data_1 <- pred_data[1,]
    pred_data_1$log_size_t0 <- row$log_size_t0
    df$surv_pred[i] <- predict(surv_mt2clim_best_mod, newdata = pred_data_1, type = "response", re.form = re_form)
    df$growth_pred[i] <- predict(growth_mt2clim_best_mod, newdata = pred_data_1, re.form = re_form)
    df$growth_var_pred[i] <- predict(growth_var_mt2clim_best_mod, newdata = pred_data_1, re.form = re_form)
    df$flower_prob_pred[i] <- predict(repro_mt2clim_best_mod, newdata = pred_data_1, type = "response", re.form = re_form)
    df$flower_count_pred[i] <- predict(nflrs_mt2clim_best_mod, newdata = pred_data_1, type = "response", re.form = re_form)
  }
  df
}

## pop_kernel(site, pop, year)
## Returns K = P + F using pop-level covariate means and pop-specific germ rate.
## Includes site.pop random effects. aspect.cat is derived from pop name.
pop_kernel <- function(site, pop, year,use.re=F) {
  aspect.cat <- case_when(
    str_starts(pop, "North") ~ "north",
    str_starts(pop, "South") ~ "south",
    str_starts(pop, "Top") ~ "top"
  )
  re_form <- if(use.re) {~(1|site.pop)} else NA
  
  log_size_center <- scale_params$center[scale_params$var == "log_size_t0"]
  log_size_scale <- scale_params$scale[scale_params$var == "log_size_t0"]

  clim_row <- clim_means_pop %>% filter(site == !!site, pop == !!pop, year == as.character(!!year))
  mt_row <- mt_means_pop %>% filter(site == !!site, pop == !!pop)

  pred_data <- data.frame(log_size_t0 = (binmids - log_size_center) / log_size_scale,
                          site = site, aspect.cat = aspect.cat,
                          year = factor(as.character(year), levels = as.character(years)),
                          site.pop = paste(site, pop, sep = "."))
  for (v in clim_vars) {
    pred_data[[v]] <- (clim_row[[v]] - scale_params$center[scale_params$var == v]) /
                       scale_params$scale[scale_params$var == v]
  }
  for (v in mt_vars) {
    pred_data[[v]] <- (mt_row[[v]] - scale_params$center[scale_params$var == v]) /
                       scale_params$scale[scale_params$var == v]
  }

  surv_prob <- predict(surv_mt2clim_best_mod, newdata = pred_data, type = "response", re.form = re_form)
  growth_mean <- binmids + predict(growth_mt2clim_best_mod, newdata = pred_data, re.form = re_form)
  growth_var <- pmax(predict(growth_var_mt2clim_best_mod, newdata = pred_data, re.form = re_form), 0)
  flower_prob <- predict(repro_mt2clim_best_mod, newdata = pred_data, type = "response", re.form = re_form)
  flower_count <- predict(nflrs_mt2clim_best_mod, newdata = pred_data, type = "response", re.form = re_form)

  pop_germ <- germ_rates_pop$germ_rate[germ_rates_pop$pop == pop & germ_rates_pop$site == site]

  growth_mx <- matrix(0, n.bin, n.bin)
  for (ss in 1:n.bin) {
    growcdf <- pnorm(vec.bin, growth_mean[ss], sqrt(growth_var[ss]))
    grows <- diff(growcdf)
    if (sum(grows) > 0) growth_mx[, ss] <- grows / sum(grows)
  }

  P_kernel <- t(t(growth_mx) * surv_prob)
  F_kernel <- outer(sdlgszprobs, flower_prob * flower_count * nutlets_per_flower * pop_germ)
  P_kernel + F_kernel
}

## pop_kernel_climOnly(site, pop, year)
## Like pop_kernel but uses climate-only models — no MT predictors.
## Returns K = P + F using pop-level climate means and pop-specific germ rate.
pop_kernel_climOnly <- function(site, pop, year, use.re = FALSE) {
  aspect.cat <- case_when(
    str_starts(pop, "North") ~ "north",
    str_starts(pop, "South") ~ "south",
    str_starts(pop, "Top") ~ "top"
  )
  re_form <- if (use.re) {~(1|site.pop)} else NA

  log_size_center <- scale_params$center[scale_params$var == "log_size_t0"]
  log_size_scale <- scale_params$scale[scale_params$var == "log_size_t0"]

  clim_row <- clim_means_pop %>% filter(site == !!site, pop == !!pop, year == as.character(!!year))

  pred_data <- data.frame(log_size_t0 = (binmids - log_size_center) / log_size_scale,
                          site = site, aspect.cat = aspect.cat,
                          year = factor(as.character(year), levels = as.character(years)),
                          site.pop = paste(site, pop, sep = "."))
  for (v in clim_vars) {
    pred_data[[v]] <- (clim_row[[v]] - scale_params$center[scale_params$var == v]) /
                       scale_params$scale[scale_params$var == v]
  }

  surv_prob <- predict(surv_climOnly_best_mod, newdata = pred_data, type = "response", re.form = re_form)
  growth_mean <- binmids + predict(growth_climOnly_best_mod, newdata = pred_data, re.form = re_form)
  growth_var <- pmax(predict(growth_var_climOnly_best_mod, newdata = pred_data, re.form = re_form), 0)
  flower_prob <- predict(repro_climOnly_best_mod, newdata = pred_data, type = "response", re.form = re_form)
  flower_count <- predict(nflrs_climOnly_best_mod, newdata = pred_data, type = "response", re.form = re_form)
  # surv_prob <- predict(surv_clim_best_mod, newdata = pred_data, type = "response", re.form = re_form)
  # growth_mean <- binmids + predict(growth_clim_best_mod, newdata = pred_data, re.form = re_form)
  # growth_var <- pmax(predict(growth_var_clim_best_mod, newdata = pred_data, re.form = re_form), 0)
  # flower_prob <- predict(repro_clim_best_mod, newdata = pred_data, type = "response", re.form = re_form)
  # flower_count <- predict(nflrs_clim_best_mod, newdata = pred_data, type = "response", re.form = re_form)

  pop_germ <- germ_rates_pop$germ_rate[germ_rates_pop$pop == pop & germ_rates_pop$site == site]

  growth_mx <- matrix(0, n.bin, n.bin)
  for (ss in 1:n.bin) {
    growcdf <- pnorm(vec.bin, growth_mean[ss], sqrt(growth_var[ss]))
    grows <- diff(growcdf)
    if (sum(grows) > 0) growth_mx[, ss] <- grows / sum(grows)
  }

  P_kernel <- t(t(growth_mx) * surv_prob)
  F_kernel <- outer(sdlgszprobs, flower_prob * flower_count * nutlets_per_flower * pop_germ)
  P_kernel + F_kernel
}

## vr_surface_df(vr, var1, var2, ...)
## Returns a data frame of vital rate predictions over a grid of two scaled predictors.
## vr: one of "surv", "growth", "growth_var", "repro", "nflrs"
## var1, var2: column names from growth_mod_data_scaled
## Non-focal continuous predictors held at pop-level means (averaged across years).
## log_size_t0 expanded at Q1/median/Q3. year fixed at 2024 (averaged out).
## n_grid: number of values per axis (default 20)
## pop_filter: subset of populations to include
vr_surface_df <- function(vr, var1, var2, n_grid = 20,
                           pop_filter = unique(mt_means_pop$pop)) {
  mod <- switch(vr,
    surv = surv_mt2clim_best_mod,
    growth = growth_mt2clim_best_mod,
    growth_var = growth_var_mt2clim_best_mod,
    repro = repro_mt2clim_best_mod,
    nflrs = nflrs_mt2clim_best_mod
  )

  scaled_data <- mod_data$growth_mod_data_scaled
  cont_vars <- c(clim_vars, mt_vars)

  var1_seq <- seq(min(scaled_data[[var1]], na.rm = TRUE),
                  max(scaled_data[[var1]], na.rm = TRUE), length.out = n_grid)
  var2_seq <- seq(min(scaled_data[[var2]], na.rm = TRUE),
                  max(scaled_data[[var2]], na.rm = TRUE), length.out = n_grid)
  size_q <- setNames(quantile(scaled_data$log_size_t0, c(0.25, 0.5, 0.75), na.rm = TRUE),
                     c("small", "medium", "large"))

  # Pop-level means averaged across years; log_size_t0 excluded (handled via size_q)
  pop_means <- scaled_data %>%
    group_by(site, pop, aspect.cat) %>%
    summarise(across(all_of(cont_vars), ~mean(., na.rm = TRUE)),
              site.pop = first(site.pop), .groups = "drop") %>%
    filter(pop %in% pop_filter) %>%
    mutate(year = factor("2024", levels = as.character(years)),
           pop_id = row_number())

  # Grid: pop * size * var1 * var2; other vars from pop means
  out <- expand_grid(pop_id = pop_means$pop_id,
                     size_cat = factor(names(size_q), levels = c("small","medium","large")),
                     !!var1 := var1_seq,
                     !!var2 := var2_seq) %>%
    mutate(log_size_t0 = size_q[as.character(size_cat)]) %>%
    left_join(pop_means %>% dplyr::select(-all_of(c(var1, var2))), by = "pop_id") %>%
    dplyr::select(-pop_id)

  if (inherits(mod, "lmerMod")) {
    out$pred <- predict(mod, newdata = out, re.form = NA)
  } else {
    out$pred <- predict(mod, newdata = out, type = "response", re.form = NA)
  }
  out
}

## lambda_surface_plot(var1, var2, ...)
## Contour plot of lambda over a grid of two continuous predictors (raw/unscaled units).
## var1, var2: column names from growth_mod_data (clim or microtopo variables)
## ref_year: year to hold fixed (default 2023)
## n_surf: grid resolution per axis (default 20)
## site_filter, pop_filter: subsets of sites/pops to include
lambda_surface_plot <- function(var1, var2, ref_year = 2023, n_surf = 20,
                                site_filter = sites, pop_filter = unique(mt_means_pop$pop)) {
  all_data <- mod_data$growth_mod_data

  var1_seq <- seq(min(all_data[[var1]], na.rm = TRUE),
                  max(all_data[[var1]], na.rm = TRUE),
                  length.out = n_surf)
  var2_seq <- seq(min(all_data[[var2]], na.rm = TRUE),
                  max(all_data[[var2]], na.rm = TRUE),
                  length.out = n_surf)

  pop_info <- mt_means_pop %>%
    filter(site %in% site_filter, pop %in% pop_filter) %>%
    dplyr::select(site, pop, aspect.cat)

  grid_df <- expand_grid(pop_info, year = ref_year,
                         !!var1 := var1_seq, !!var2 := var2_seq) %>%
    predict_lambda()

  ggplot(grid_df, aes(x = .data[[var1]], y = .data[[var2]], fill = lambda)) +
    geom_raster() +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#D6604D", midpoint = 1) +
    facet_grid(pop ~ site) +
    labs(x = var1, y = var2, fill = expression(lambda),
         title = paste0("Lambda surface: ", var1, " \u00d7 ", var2, " (year = ", ref_year, ")")) +
    theme_bw()
}
##### End 4 #####

