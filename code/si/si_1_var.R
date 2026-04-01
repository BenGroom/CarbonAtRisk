# SI Figure 1: Value at Risk (VaR) Definition
# Extracted from code/fig1/figure1.R panel_a

library(ggplot2)

# Parameters ---------------------------------------------------------------

VAR_MEAN <- 0
VAR_SD   <- 5
confidence_levels <- c(0.90, 0.95, 0.98)
BLUE_SHADES <- c("0.9" = "#9ecae1", "0.95" = "#4292c6", "0.98" = "#08519c")

# Output directory ---------------------------------------------------------

out_dir <- "outputs/si"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Generate PDF data --------------------------------------------------------

x_var <- seq(VAR_MEAN - 4 * VAR_SD, VAR_MEAN + 4 * VAR_SD, length.out = 1000)
y_var <- dnorm(x_var, VAR_MEAN, VAR_SD)
df_var <- data.frame(x = x_var, y = y_var)

# VaR thresholds (left-tail quantiles)
var_thresholds <- qnorm(1 - confidence_levels, VAR_MEAN, VAR_SD)

# Shading data: shade from left edge to each threshold
shade_var <- do.call(rbind, lapply(seq_along(confidence_levels), function(i) {
  cl <- confidence_levels[i]
  thresh <- var_thresholds[i]
  xs <- seq(min(x_var), thresh, length.out = 300)
  data.frame(
    x = xs,
    y = dnorm(xs, VAR_MEAN, VAR_SD),
    level = as.character(cl)
  )
}))
shade_var$level <- factor(shade_var$level, levels = c("0.9", "0.95", "0.98"))

# Build plot ---------------------------------------------------------------

y_max_a <- max(y_var)

panel_a <- ggplot() +
  geom_area(
    data = shade_var, aes(x = x, y = y, fill = level),
    position = "identity", alpha = 0.7
  ) +
  geom_line(data = df_var, aes(x = x, y = y), linewidth = 0.6) +
  scale_fill_manual(
    values = BLUE_SHADES,
    labels = paste0(
      confidence_levels * 100, "% VaR = ",
      sprintf("%.2f", var_thresholds), "%"
    ),
    name = NULL
  ) +
  # Vertical dashed lines at thresholds
  annotate(
    "segment",
    x = var_thresholds, xend = var_thresholds,
    y = 0, yend = y_max_a * 1.05,
    linetype = "dashed", color = unname(BLUE_SHADES), linewidth = 0.4
  ) +
  # Vertical line at x=0
  annotate("segment", x = 0, xend = 0,
           y = 0, yend = y_max_a * 1.05,
           linetype = "dotted", color = "grey50", linewidth = 0.5) +
  # VaR arrows from 0 to each quantile
  annotate(
    "segment",
    x = rep(0, 3), xend = var_thresholds,
    y = y_max_a * c(0.30, 0.18, 0.06),
    yend = y_max_a * c(0.30, 0.18, 0.06),
    arrow = arrow(length = unit(0.15, "cm")),
    color = unname(BLUE_SHADES), linewidth = 0.8
  ) +
  annotate(
    "text",
    x = rep(var_thresholds[1] / 2, 3),
    y = y_max_a * c(0.34, 0.22, 0.10),
    label = paste0(confidence_levels * 100, "% VaR"),
    size = 5, color = unname(BLUE_SHADES), fontface = "bold",
    hjust = 0.5
  ) +
  # Below-zero: Loss/Profit arrow
  annotate("segment", x = -14, xend = 14,
           y = -y_max_a * 0.15, yend = -y_max_a * 0.15,
           arrow = arrow(length = unit(0.15, "cm"), ends = "both"),
           color = "grey50", linewidth = 0.6) +
  annotate("text", x = -15, y = -y_max_a * 0.15,
           label = "Loss", size = 5, hjust = 1,
           color = "grey50", fontface = "italic") +
  annotate("text", x = 15, y = -y_max_a * 0.15,
           label = "Profit", size = 5, hjust = 0,
           color = "grey50", fontface = "italic") +
  # Callout
  annotate("text", x = 0, y = -y_max_a * 0.28,
           label = 'italic("How much money could I lose?")',
           parse = TRUE, size = 5.5, color = "#4292c6", hjust = 0.5) +
  scale_x_continuous(
    breaks = seq(-20, 20, by = 5),
    labels = function(x) paste0(x, "%")
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_y_continuous(breaks = 0, labels = "0") +
  coord_cartesian(xlim = c(-18, 18),
                  ylim = c(-y_max_a * 0.35, y_max_a * 1.05),
                  clip = "off") +
  labs(
    x = "Return on Investment (%)",
    y = "Probability Density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9),
    legend.background = element_rect(fill = alpha("white", 0.9),
                                     color = "black", linewidth = 0.4),
    plot.margin = margin(5, 10, 5, 5)
  ) +
  guides(fill = guide_legend(override.aes = list(color = "grey40",
                                                  linewidth = 0.3)))

# Save --------------------------------------------------------------------

ggsave(file.path(out_dir, "si_1_var.pdf"), panel_a,
       width = 8, height = 5)

cat("Saved:", file.path(out_dir, "si_1_var.pdf"), "\n")
