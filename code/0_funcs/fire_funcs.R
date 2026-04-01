# Shared simulation functions for Carbon at Risk (CaR) fire risk analysis
# Adapted from code/funcs/fire_funcs.R for standalone  pipeline

# HELPER FUNCTIONS------------------------------------------------------------

#' Fetch fire data from EFFIS API
#' @param geo_id The geographic ID (e.g., "US.5_1")
#' @param cache_dir Optional directory to cache API responses as CSV
#' @return A data frame with annual burned area data
fetch_fire_data <- function(geo_id, cache_dir = NULL) {
  # Check cache first
  if (!is.null(cache_dir)) {
    cache_file <- file.path(cache_dir, sprintf("effis_fire_%s.csv", gsub("\\.", "_", geo_id)))
    if (file.exists(cache_file)) {
      cat(sprintf("  Loading cached fire data for %s\n", geo_id))
      return(read.csv(cache_file, stringsAsFactors = FALSE))
    }
  }

  fires_url <- sprintf(
    "https://cprof.effis.emergency.copernicus.eu/api/v3/banf?level=ADM1&value=%s&year=2023&yearFrom=2002&yearTo=2023&env=PROD",
    geo_id
  )

  response <- GET(fires_url)
  if (status_code(response) != 200) {
    warning(sprintf("Failed to fetch fire data for %s", geo_id))
    return(NULL)
  }

  fires_json <- content(response, as = "parsed", type = "application/json")
  fires_df <- bind_rows(fires_json$banfyear)

  # Save to cache
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(fires_df, cache_file, row.names = FALSE)
    cat(sprintf("  Cached fire data to %s\n", cache_file))
  }

  return(fires_df)
}

#' Fetch forest indicators from EFFIS API
#' @param geo_id The geographic ID
#' @param cache_dir Optional directory to cache API responses
#' @return A list with forest indicators
fetch_forest_indicators <- function(geo_id, cache_dir = NULL) {
  # Check cache first
  if (!is.null(cache_dir)) {
    cache_file <- file.path(cache_dir, sprintf("effis_forest_%s.rds", gsub("\\.", "_", geo_id)))
    if (file.exists(cache_file)) {
      cat(sprintf("  Loading cached forest indicators for %s\n", geo_id))
      return(readRDS(cache_file))
    }
  }

  forest_url <- sprintf(
    "https://cprof.effis.emergency.copernicus.eu/api/v3/indicators?level=ADM1&value=%s&env=PROD",
    geo_id
  )

  response <- GET(forest_url)
  if (status_code(response) != 200) {
    warning(sprintf("Failed to fetch forest indicators for %s", geo_id))
    return(NULL)
  }

  forest_json <- content(response, as = "parsed", type = "application/json")

  # Save to cache
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(forest_json, cache_file)
    cat(sprintf("  Cached forest indicators to %s\n", cache_file))
  }

  return(forest_json)
}

#' Generate correlated burn fractions using Gaussian copula
#' @param n_simulations Number of Monte Carlo draws
#' @param n_projects Number of projects
#' @param empirical_burn_fractions Vector of historical burn fractions
#' @param correlation Common pairwise correlation (0 to 1)
#' @return Matrix (n_simulations x n_projects) of correlated burn fractions
sample_burn_copula <- function(n_simulations, n_projects,
                               empirical_burn_fractions, correlation) {
  # Fast equicorrelated normals: X = sqrt(rho)*Z_common + sqrt(1-rho)*Z_indep
  z_common <- rnorm(n_simulations)
  z_indep <- matrix(rnorm(n_simulations * n_projects), n_simulations, n_projects)
  Z <- sqrt(correlation) * z_common + sqrt(1 - correlation) * z_indep

  # Transform to uniform [0,1] via normal CDF
  U <- pnorm(Z)

  # Fast vectorized quantile mapping using linear interpolation
  sorted_vals <- sort(empirical_burn_fractions)
  n <- length(sorted_vals)
  probs <- (1:n - 0.5) / n

  u_flat <- as.vector(U)
  result <- approx(probs, sorted_vals, xout = u_flat, rule = 2)$y
  matrix(result, nrow = n_simulations, ncol = n_projects)
}

#' Simulate burn dynamics for K projects over a time horizon
#' VECTORIZED: runs all n_simulations at once using matrix operations
#' @param empirical_burn_fractions Vector of historical burn fractions
#' @param regrowth_rate Annual regrowth rate
#' @param climate_rate Annual climate-driven increase in burn rate
#' @param max_horizon Maximum number of years to simulate
#' @param n_projects Number of projects (default 1)
#' @param n_simulations Number of Monte Carlo simulations
#' @param correlation Pairwise correlation between projects (default 0)
#' @return Matrix (n_simulations x max_horizon+1): cumulative burned fraction
simulate_burn_dynamics <- function(empirical_burn_fractions,
                                   regrowth_rate,
                                   climate_rate,
                                   max_horizon,
                                   n_projects = 1,
                                   n_simulations = 1000,
                                   correlation = 0) {
  result <- matrix(0.0, nrow = n_simulations, ncol = max_horizon + 1)

  if (max_horizon == 0) return(result)

  unburned <- matrix(1.0, nrow = n_simulations, ncol = n_projects)
  burned   <- matrix(0.0, nrow = n_simulations, ncol = n_projects)

  for (year in 1:max_horizon) {
    if (correlation == 0) {
      base_burn <- matrix(
        sample(empirical_burn_fractions, n_simulations * n_projects, replace = TRUE),
        nrow = n_simulations, ncol = n_projects
      )
    } else if (correlation == 1) {
      base_burn <- matrix(
        rep(sample(empirical_burn_fractions, n_simulations, replace = TRUE), n_projects),
        nrow = n_simulations, ncol = n_projects
      )
    } else {
      base_burn <- sample_burn_copula(
        n_simulations, n_projects, empirical_burn_fractions, correlation
      )
    }
    burn_now <- pmin(base_burn * (1 + climate_rate * (year - 1)), 1.0)

    new_burn <- unburned * burn_now
    regrowth <- burned * regrowth_rate

    unburned <- unburned - new_burn + regrowth
    burned   <- burned + new_burn - regrowth

    unburned <- pmin(pmax(unburned, 0), 1)
    burned   <- pmin(pmax(burned, 0), 1)

    result[, year + 1] <- rowMeans(1 - unburned)
  }

  return(result)
}

#' Run Monte Carlo simulation for Carbon at Risk
#' @param empirical_burn_fractions Vector of historical burn fractions
#' @param annual_regrowth_rate Annual regrowth rate
#' @param time_horizons Vector of time horizons to evaluate
#' @param n_simulations Number of Monte Carlo simulations
#' @param climate_rate Annual climate-driven increase in burn rate
#' @param n_projects Number of projects (default 1)
#' @param correlation Pairwise correlation (default 0)
#' @return Data frame with time horizon, mean, and 95th percentile CaR
run_car_simulation <- function(empirical_burn_fractions,
                               annual_regrowth_rate,
                               time_horizons = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
                                                 110, 120, 130, 140, 150, 160, 170, 180, 200),
                               n_simulations = 1000,
                               climate_rate = 0.005,
                               n_projects = 1,
                               correlation = 0) {

  max_horizon <- max(time_horizons)

  all_results <- simulate_burn_dynamics(
    empirical_burn_fractions,
    annual_regrowth_rate,
    climate_rate,
    max_horizon = max_horizon,
    n_projects = n_projects,
    n_simulations = n_simulations,
    correlation = correlation
  )

  mean_car <- colMeans(all_results)[time_horizons + 1]
  car_95 <- apply(all_results, 2, quantile, 0.95)[time_horizons + 1]

  return(data.frame(
    time_horizon = time_horizons,
    mean_car = mean_car,
    car_95 = as.numeric(car_95)
  ))
}

#' Calculate burn fraction from fire data
#' @param fires_df Fire data frame
#' @param forest_lc1 Total forest area from indicators
#' @param project_area Project area in hectares (default 100000)
#' @param rescale_firesize Whether to apply fire-size rescaling (default TRUE)
#' @return Vector of burn fractions
calculate_burn_fractions <- function(fires_df, forest_lc1, project_area = 100000,
                                     rescale_firesize = TRUE) {
  if (rescale_firesize) {
    fires_df <- fires_df %>%
      mutate(
        ratio_stretched = {
          ratio <- firesize / project_area
          1 + (ratio - min(ratio)) * (3 - 1) / (max(ratio) - min(ratio) + 1e-10)
        },
        burn_fraction = (lc1 / forest_lc1) * ratio_stretched
      )
  } else {
    fires_df <- fires_df %>%
      mutate(burn_fraction = lc1 / forest_lc1)
  }

  return(fires_df$burn_fraction)
}

# MAIN PROCESSING FUNCTIONS -----------------

#' Process selected regions: fetch data, compute burn fractions, run CaR
#' @param gdf_gadm GeoDataFrame with admin boundaries and GPP
#' @param selected_regions List of list(country, subcountry)
#' @param CLIMATE_RATE Climate rate (default 0.005)
#' @param rescale_firesize Whether to rescale fire size (default TRUE)
#' @param project_area Project area in hectares (default 100000)
#' @param regrowth_rates Named vector of annual regrowth rates (from get_regrowth_rates())
#' @param n_simulations Number of MC simulations (default 1000)
#' @param time_horizons Vector of time horizons
#' @param cache_dir Optional cache directory for EFFIS API responses
#' @return List with car_results, plots, geometries, burn_fractions, regrowth
process_selected_geometries <- function(gdf_gadm, selected_regions,
                                         regrowth_rates,
                                         CLIMATE_RATE = 0.005,
                                         rescale_firesize = TRUE,
                                         project_area = 100000,
                                         n_simulations = 1000,
                                         time_horizons = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
                                                           110, 120, 130, 140, 150, 160, 170, 180, 200),
                                         cache_dir = NULL) {

  selected_geometries <- list()

  for (region in selected_regions) {
    selected_geom <- gdf_gadm %>%
      filter(NAME_0 == region$country, NAME_1 == region$subcountry)

    if (nrow(selected_geom) > 0) {
      selected_geometries[[length(selected_geometries) + 1]] <- selected_geom[1, ]
    } else {
      warning(sprintf("Region not found: %s - %s", region$country, region$subcountry))
    }
  }

  if (length(selected_geometries) == 0) {
    stop("No valid regions selected")
  }

  car_results <- list()
  all_plots <- list()
  all_burn_fractions <- list()
  all_regrowth <- list()

  cat("Processing", length(selected_geometries), "regions\n")

  for (i in seq_along(selected_geometries)) {
    row <- selected_geometries[[i]]

    geo_id    <- paste0(row$country_code, ".", row$ID_1, "_1")
    geo_label <- paste(row$NAME_0, "-", row$NAME_1)

    cat(sprintf("----- %s (%s) -----\n", geo_label, geo_id))

    region_name <- row$NAME_1
    if (!region_name %in% names(regrowth_rates)) {
      stop(sprintf("No regrowth rate defined for region '%s'", region_name))
    }
    annual_regrowth_rate <- regrowth_rates[[region_name]]
    cat(sprintf("  Regrowth rate (Zang): %.2f%%\n", annual_regrowth_rate * 100))
    all_regrowth[[geo_label]] <- annual_regrowth_rate

    cat("  Fetching fire data from EFFIS API...\n")
    fires_df    <- fetch_fire_data(geo_id, cache_dir = cache_dir)
    forest_json <- fetch_forest_indicators(geo_id, cache_dir = cache_dir)

    if (is.null(fires_df) || is.null(forest_json)) {
      warning(sprintf("Skipping %s due to API errors", geo_label))
      next
    }

    fires_df <- fires_df %>% arrange(year)

    p <- plot_burned_area(fires_df, geo_label)
    all_plots[[geo_label]] <- p

    empirical_burn_fractions <- calculate_burn_fractions(
      fires_df,
      forest_json$lc1,
      project_area,
      rescale_firesize = rescale_firesize
    )
    all_burn_fractions[[geo_label]] <- empirical_burn_fractions

    cat(sprintf("  Burn fraction range: %.4f - %.4f\n",
                min(empirical_burn_fractions), max(empirical_burn_fractions)))

    cat("  Running Monte Carlo simulation...\n")
    car_df <- run_car_simulation(
      empirical_burn_fractions,
      annual_regrowth_rate,
      time_horizons,
      n_simulations,
      CLIMATE_RATE
    )

    car_df$geo_label <- geo_label
    car_df$regrowth_rate <- annual_regrowth_rate

    car_results[[geo_label]] <- car_df

    cat(sprintf("  CaR at 1yr: %.0f kg/tonne\n", car_df$car_95[car_df$time_horizon == 1] * 1000))
    cat(sprintf("  CaR at 100yr: %.0f kg/tonne\n", car_df$car_95[car_df$time_horizon == 100] * 1000))
    cat("\n")
  }

  car_combined <- bind_rows(car_results)

  all_regrowth <- bind_rows(
    lapply(names(all_regrowth), function(name) {
      tibble(
        geo_label = name,
        regrowth_rate = all_regrowth[[name]]
      )
    })
  )

  all_burn_fractions <- bind_rows(
    lapply(names(all_burn_fractions), function(name) {
      tibble(
        year = 1:length(all_burn_fractions[[name]]),
        geo_label = name,
        burn_fraction = all_burn_fractions[[name]]
      )
    })
  )

  return(list(
    car_results = car_combined,
    plots = all_plots,
    geometries = selected_geometries,
    empirical_burn_fractions = all_burn_fractions,
    regrowth = all_regrowth
  ))
}


#' Run diversification Monte Carlo simulation
#' @param empirical_burn_fractions Vector of historical burn fractions
#' @param regrowth_rate Annual regrowth rate
#' @param climate_rate Annual climate-driven increase in burn rate
#' @param years Vector of time horizons
#' @param n_simulations Number of Monte Carlo simulations
#' @param n_projects Number of projects
#' @param estate_area Total estate area in hectares
#' @param correlation Pairwise correlation (default 0)
#' @param return_raw Whether to return raw MC results (default FALSE)
#' @return List with mean_loss and car_95 vectors (in hectares)
simulate_diversified_car <- function(empirical_burn_fractions,
                                     regrowth_rate,
                                     climate_rate,
                                     years,
                                     n_simulations,
                                     n_projects,
                                     estate_area,
                                     correlation = 0,
                                     return_raw = FALSE) {

  max_year <- max(years)

  all_results <- simulate_burn_dynamics(
    empirical_burn_fractions,
    regrowth_rate,
    climate_rate,
    max_horizon = max_year,
    n_projects = n_projects,
    n_simulations = n_simulations,
    correlation = correlation
  )

  mc_hectares <- all_results * estate_area

  mean_loss <- colMeans(mc_hectares)[years + 1]
  car_95 <- apply(mc_hectares, 2, quantile, 0.95)[years + 1]

  result <- list(mean_loss = as.numeric(mean_loss), car_95 = as.numeric(car_95))
  if (return_raw) {
    result$mc_hectares <- mc_hectares
  }
  return(result)
}


#' Run full diversification analysis for a region
#' @param gdf_gadm GeoDataFrame with admin boundaries
#' @param country Country name (default "United States")
#' @param subcountry Subcountry name (default "California")
#' @param estate_area Total estate area in hectares (default 1000)
#' @param n_simulations Number of Monte Carlo simulations (default 500)
#' @param max_years Maximum time horizon (default 200)
#' @param correlation Correlation between projects (default 0)
#' @param return_raw Whether to return raw MC results (default FALSE)
#' @param regrowth_rates Named vector of annual regrowth rates (from get_regrowth_rates())
#' @param climate_rate Annual climate increase rate (default 0.005)
#' @param rescale_firesize Whether to rescale fire size (default TRUE)
#' @param cache_dir Optional cache directory for EFFIS API responses
#' @return List with plot, years, results, region, estate_area, summary
run_diversification_analysis <- function(gdf_gadm,
                                         regrowth_rates,
                                         country = "United States",
                                         subcountry = "California",
                                         estate_area = 1000,
                                         n_simulations = 500,
                                         max_years = 200,
                                         correlation = 0,
                                         return_raw = FALSE,
                                         climate_rate = 0.005,
                                         rescale_firesize = TRUE,
                                         cache_dir = NULL) {

  row <- gdf_gadm %>%
    filter(NAME_0 == country, NAME_1 == subcountry) %>%
    slice(1)

  geo_id    <- paste0(row$country_code, ".", row$ID_1, "_1")
  geo_label <- paste(country, "-", subcountry)
  fires_df  <- fetch_fire_data(geo_id, cache_dir = cache_dir)
  forest_json <- fetch_forest_indicators(geo_id, cache_dir = cache_dir)

  if (is.null(fires_df) || is.null(forest_json)) stop("Failed to fetch fire data from API")

  empirical_burn_fractions <- calculate_burn_fractions(
    fires_df, forest_json$lc1, project_area = estate_area,
    rescale_firesize = rescale_firesize
  )

  if (!subcountry %in% names(regrowth_rates)) {
    stop(sprintf("No regrowth rate defined for region '%s'", subcountry))
  }
  regrowth_rate <- regrowth_rates[[subcountry]]
  cat(sprintf("Regrowth rate (Zang): %.2f%%\n", regrowth_rate * 100))

  time_horizons <- c(1, seq(10, max_years, by = 10))
  years <- 0:max_years

  results <- purrr::map(
    purrr::set_names(c(1, 10, 100)),
    function(n_projects) {
      cat(sprintf("Running simulations for K=%d projects...\n", n_projects))
      simulate_diversified_car(
        empirical_burn_fractions, regrowth_rate, climate_rate,
        time_horizons, n_simulations, n_projects, estate_area,
        correlation = correlation,
        return_raw = return_raw
      )
    }
  )

  cat("Generating diversification plot...\n")
  div_plot <- plot_diversification(
    time_horizons,
    results[['1']]$mean_loss,
    results[['1']]$car_95,
    results[['10']]$car_95,
    results[['100']]$car_95,
    estate_area,
    geo_label
  )

  cat("\n--- Summary at Key Time Horizons ---\n")
  idx_50  <- which(time_horizons == 50)
  idx_100 <- which(time_horizons == 100)
  idx_200 <- which(time_horizons == 200)
  summary_df <- data.frame(
    Horizon  = c(50, 100, 200),
    Expected = round(results[['1']]$mean_loss[c(idx_50, idx_100, idx_200)], 1),
    CaR_K1   = round(results[['1']]$car_95[c(idx_50, idx_100, idx_200)], 1),
    CaR_K10  = round(results[['10']]$car_95[c(idx_50, idx_100, idx_200)], 1),
    CaR_K100 = round(results[['100']]$car_95[c(idx_50, idx_100, idx_200)], 1)
  )
  summary_df$Reduction_K100 <- paste0(
    round((1 - summary_df$CaR_K100 / summary_df$CaR_K1) * 100, 1), "%"
  )

  return(list(
    plot = div_plot,
    years = years,
    results = results,
    region = geo_label,
    estate_area = estate_area,
    summary = summary_df
  ))
}


# PLOTTING FUNCTIONS -----------------

#' Plot annual burned area bar chart
#' @param fires_df Data frame with fire data
#' @param label Region label
#' @return A ggplot object
plot_burned_area <- function(fires_df, label) {
  fires_df %>%
    arrange(year) %>%
    ggplot(aes(x = factor(year), y = ba_area_ha)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(
      title = sprintf("Annual Burned Area in %s", label),
      x = "Year",
      y = "Burned Area (ha)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
}


#' Plot combined CaR curves for all regions
#' @param car_results Combined CaR results data frame
#' @return A ggplot object
plot_car_curves <- function(car_results) {
  car_results %>%
    mutate(
      car_kg_tonne = car_95 * 1000,
      regrowth_pct = round(regrowth_rate * 100)
    ) %>%
    ggplot(aes(x = time_horizon, y = car_kg_tonne,
               color = geo_label, group = geo_label)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_y_continuous(limits = c(0, 1050)) +
    labs(
      title = "Carbon at Risk (95th percentile burned)",
      x = "Time Horizon (years)",
      y = "Carbon at Risk (kg/tonne)",
      color = "Region"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(color = guide_legend(nrow = 2))
}

#' Create CaR summary table
#' @param car_results Combined CaR results data frame
#' @return A data frame formatted as a summary table
create_car_table <- function(car_results) {
  car_results %>%
    filter(time_horizon %in% c(1, 10, 100, 200)) %>%
    mutate(car_kg_tonne = round(car_95 * 1000)) %>%
    dplyr::select(geo_label, time_horizon, car_kg_tonne) %>%
    pivot_wider(
      names_from = time_horizon,
      values_from = car_kg_tonne,
      names_prefix = "CaR_"
    ) %>%
    rename(Project = geo_label) %>%
    rename_with(~ gsub("CaR_", "CaR at ", .x), starts_with("CaR_")) %>%
    rename_with(~ paste0(.x, "yr"), starts_with("CaR at "))
}

#' Plot correlation vs distance for spatial fire correlation analysis
#' @param pairs_df Data frame with columns: distance_km, cor_spearman
#' @param rho0 Fitted intercept (correlation at distance 0)
#' @param lambda_km Fitted decay length in km
#' @param region_label Label for the region (optional)
#' @param show_fit Whether to show the exponential fit curve (default TRUE)
#' @param use_ribbon Whether to show ribbon with quantiles (default FALSE)
#' @param bin_width Distance bin width in km for ribbon (default 100)
#' @param ca_color Color for data (default "#E67E22")
#' @param base_size Base font size (default 22)
#' @return A ggplot object
plot_correlation_vs_distance <- function(pairs_df, rho0, lambda_km,
                                          region_label = NULL,
                                          show_fit = TRUE,
                                          use_ribbon = FALSE,
                                          bin_width = 100,
                                          ca_color = "#E67E22",
                                          base_size = 22) {

  max_dist <- max(pairs_df$distance_km, na.rm = TRUE)

  p <- ggplot() +
    geom_hline(yintercept = 0, color = "gray50", linetype = "solid", linewidth = 0.3) +
    labs(
      x = "Distance (km)",
      y = "Spearman correlation"
    ) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )

  if (use_ribbon) {
    ribbon_df <- pairs_df %>%
      mutate(dist_bin = floor(distance_km / bin_width) * bin_width + bin_width / 2) %>%
      group_by(dist_bin) %>%
      summarise(
        cor_median = median(cor_spearman, na.rm = TRUE),
        cor_q25 = quantile(cor_spearman, 0.25, na.rm = TRUE),
        cor_q75 = quantile(cor_spearman, 0.75, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      filter(n >= 3)

    p <- p +
      geom_ribbon(
        data = ribbon_df,
        aes(x = dist_bin, ymin = cor_q25, ymax = cor_q75),
        fill = ca_color, alpha = 0.3
      ) +
      geom_line(
        data = ribbon_df,
        aes(x = dist_bin, y = cor_median),
        color = ca_color, linewidth = 1
      ) +
      geom_point(
        data = ribbon_df,
        aes(x = dist_bin, y = cor_median),
        color = ca_color, size = 2
      )
  } else {
    p <- p +
      geom_point(
        data = pairs_df,
        aes(x = distance_km, y = cor_spearman),
        alpha = 0.4, size = 1.5, color = ca_color
      )
  }

  if (show_fit && !is.na(rho0) && !is.na(lambda_km)) {
    fit_df <- data.frame(
      distance_km = seq(0, max_dist, length.out = 100)
    )
    fit_df$cor_fit <- rho0 * exp(-fit_df$distance_km / lambda_km)

    p <- p +
      geom_line(
        data = fit_df,
        aes(x = distance_km, y = cor_fit),
        color = "gray30", linewidth = 0.8, linetype = "dashed"
      ) +
      annotate(
        "text",
        x = max_dist * 0.45, y = rho0 * 0.85,
        label = sprintf("rho(d) == %.2f * e^{-d/%.0f~km}", rho0, lambda_km),
        parse = TRUE, hjust = 0, size = 2
      )
  }

  if (!is.null(region_label)) {
    p <- p + ggtitle(sprintf("Spatial Correlation: %s", region_label))
  }

  return(p)
}


#' Plot diversification analysis chart
#' @param years Vector of time horizons
#' @param mean_loss Expected value (mean loss)
#' @param car_1 95% CaR for K=1 project
#' @param car_10 95% CaR for K=10 projects
#' @param car_100 95% CaR for K=100 projects
#' @param estate_area Total estate area
#' @param region_label Label for the region
#' @return A ggplot object
plot_diversification <- function(years, mean_loss, car_1, car_10, car_100,
                                 estate_area, region_label) {

  df <- data.frame(
    year = rep(years, 4),
    value = c(mean_loss, car_1, car_10, car_100),
    series = factor(rep(c("Expected Value", "1 Project",
                          "10 Projects", "100 Projects"),
                        each = length(years)),
                    levels = c("Expected Value", "1 Project",
                               "10 Projects", "100 Projects"))
  )

  ggplot(df, aes(x = year, y = value, color = series)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c("Expected Value" = "gray30",
                                  "1 Project" = "#E67E22",
                                  "10 Projects" = "#7F8C8D",
                                  "100 Projects" = "#95A5A6")) +
    scale_y_continuous(limits = c(0, estate_area)) +
    labs(
      title = NULL,
      x = "Time Horizon (years)",
      y = "CaR",
      color = NULL
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 15),
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.key.size = unit(0.6, "cm"),
      legend.text = element_text(size = 12),
      legend.background = element_rect(fill = "white", color = NA)
    )
}
