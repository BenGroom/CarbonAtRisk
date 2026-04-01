# Portfolio functions for cross-technology CaR analysis
# Bernoulli project model with within- and between-technology correlation
# Adapted from code/funcs/fig4.R for standalone ncc use

library(truncnorm)

#' Compute portfolio statistics for a grid of (n_DACS, n_Forest) combinations
#'
#' @param df Data frame with columns n_DACS, n_Forest
#' @param q Tonnes per project if successful
#' @param p Named vector of survival probabilities (DACS, Forest)
#' @param costs Named vector of project costs (DACS, Forest)
#' @param z z-score for CaR quantile (e.g., qnorm(0.05) for 95% CaR)
#' @param rho_within Named vector of within-technology correlations
#' @param rho_between Between-technology correlation scalar
portfolio_stats <- function(df, q, p, costs, z,
                            rho_within  = c(DACS = 0, Forest = 0),
                            rho_between = 0) {
  df %>%
    mutate(
      n_total = n_DACS + n_Forest,
      mu = q * (n_DACS * p["DACS"] + n_Forest * p["Forest"]),
      v_DACS   = q^2 * p["DACS"]   * (1 - p["DACS"]),
      v_Forest = q^2 * p["Forest"] * (1 - p["Forest"]),
      var_DACS = v_DACS   * (n_DACS   + n_DACS   * pmax(n_DACS   - 1, 0) * rho_within["DACS"]),
      var_F    = v_Forest * (n_Forest + n_Forest * pmax(n_Forest - 1, 0) * rho_within["Forest"]),
      cov_DF   = (n_DACS * n_Forest) * rho_between * sqrt(v_DACS * v_Forest),
      var = var_DACS + var_F + 2 * cov_DF,
      sd  = sqrt(pmax(var, 0)),
      p5  = mu + z * sd,
      CaR95 = n_total * q - p5,
      cost = n_DACS * costs["DACS"] + n_Forest * costs["Forest"],
      cost_per_expected_tonne = cost / mu,
      forest_share = if_else(n_total > 0, n_Forest / n_total, NA_real_),
      lower = 0,
      upper = n_total * q
    ) %>%
    select(-v_DACS, -v_Forest, -var_DACS, -var_F, -cov_DF)
}

#' Find the minimum-cost feasible portfolio
pick_min_cost_feasible <- function(df, target, target_var = "p5") {
  if (target_var == "p5") {
    df %>%
      filter(p5 >= target, n_total > 0) %>%
      slice_min(order_by = cost, with_ties = FALSE)
  } else if (target_var == "CaR95") {
    df %>%
      filter(CaR95 <= target, n_total > 0) %>%
      filter(n_total == max(n_total))
  } else {
    stop("target_var must be 'p5' or 'CaR95'")
  }
}

#' Get optimal bundles (All DACS, All Forest, 50:50 mix) for a given scenario
get_bundles <- function(opts, q, survival_probs, project_costs, z_p5,
                        target, target_var,
                        rho_between = 0,
                        rho_within = c(DACS = 0, Forest = 0)) {
  purrr::map_dfr(opts,
    function(grid) {
      portfolio_stats(grid,
        q = q, p = survival_probs, costs = project_costs,
        z = z_p5,
        rho_between = rho_between, rho_within = rho_within
      ) %>%
        pick_min_cost_feasible(target = target, target_var = target_var)
    },
    .id = "label"
  )
}

#' Generate PDF data for truncated normal densities
get_pdf <- function(df,
                    label_levels = c("All Forest", "50:50 mix", "All DACS")) {
  x_max <- max(df$upper)
  df %>%
    tidyr::crossing(x = seq(0, x_max, length.out = 1000)) %>%
    mutate(
      density = dtruncnorm(x, a = lower, b = upper, mean = mu, sd = sd),
      label = factor(label, levels = label_levels),
      legend = paste0(
        label, " (", n_DACS, "D + ", n_Forest, "F, $",
        scales::comma(cost), ")"
      )
    ) %>%
    mutate(legend = fct_reorder(legend, mu), .by = "scenario")
}
