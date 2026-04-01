# Figure 1: CaR Definition and Buffer Interpretation
# Panel (a): CaR definition (truncated normal with CaR arrows + decomposition)
# Panel (b): Buffer interpretation diagram (contract bar + delivery distribution)
#
# The buffer interpretation follows notes/buffer_interpretation.tex:
#   CaR = Buffer, by definition, when both computed at contracted scale Q.
#   Contract Q = Q* + CaR_95, so p_5(Q) = Q*.
#   The p_5 calculation is on ALL Q tonnes (target + buffer), so
#   buffer-on-buffer failure is already accounted for.

# Libraries ----------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(truncnorm)
library(grid)

# Output directory ---------------------------------------------------------

out_dir <- "outputs/main"
subfig_dir <- file.path(out_dir, "subfigs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(subfig_dir)) dir.create(subfig_dir, recursive = TRUE)

# Parameters ---------------------------------------------------------------

CAR_MEAN  <- 850
CAR_SD    <- 150
CAR_LOWER <- 0
CAR_UPPER <- Inf   # no upper truncation: delivery can exceed Q
Q_TARGET  <- 1000  # Q: total contracted removal

# Expected value of truncated normal
CAR_MU <- etruncnorm(a = CAR_LOWER, b = CAR_UPPER, mean = CAR_MEAN, sd = CAR_SD)

# Confidence level for main display
ALPHA <- 0.95

# Colors
GREEN_SHADES <- c("0.9" = "#b2df8a", "0.95" = "#33a02c", "0.98" = "#1b7837")
RED_SHADES   <- c("0.9" = "#fcae91", "0.95" = "#fb6a4a", "0.98" = "#cb181d")
GREEN_95 <- "#33a02c"
RED_95   <- "#fb6a4a"
DARK_RED <- "#8b0000"

# Confidence levels
confidence_levels <- c(0.90, 0.95, 0.98)

# Generate distribution data -----------------------------------------------

x_car <- seq(CAR_LOWER, Q_TARGET + 2 * CAR_SD, length.out = 1000)
y_car <- dtruncnorm(x_car, a = CAR_LOWER, b = CAR_UPPER,
                     mean = CAR_MEAN, sd = CAR_SD)
df_car <- data.frame(x = x_car, y = y_car)

# CaR thresholds for all three confidence levels
car_thresholds <- qtruncnorm(1 - confidence_levels,
                              a = CAR_LOWER, b = CAR_UPPER,
                              mean = CAR_MEAN, sd = CAR_SD)
car_values <- Q_TARGET - car_thresholds

# Shorthand for 95%
p5 <- car_thresholds[2]
car_95 <- car_values[2]

y_max <- max(y_car)

# Shading data for all three levels
shade_car <- do.call(rbind, lapply(seq_along(confidence_levels), function(i) {
  cl <- confidence_levels[i]
  thresh <- car_thresholds[i]
  xs <- seq(CAR_LOWER + 0.1, thresh, length.out = 300)
  data.frame(
    x = xs,
    y = dtruncnorm(xs, a = CAR_LOWER, b = CAR_UPPER,
                    mean = CAR_MEAN, sd = CAR_SD),
    level = as.character(cl)
  )
}))
shade_car$level <- factor(shade_car$level, levels = c("0.9", "0.95", "0.98"))

# Shading data: tail region below p5 (for panel b)
shade_tail <- data.frame(
  x = seq(CAR_LOWER + 0.1, p5, length.out = 300),
  y = dtruncnorm(seq(CAR_LOWER + 0.1, p5, length.out = 300),
                  a = CAR_LOWER, b = CAR_UPPER,
                  mean = CAR_MEAN, sd = CAR_SD)
)

# "Panel (a): CaR Definition" ---------------------------------------------

# Arrow y-positions for 90%, 95%, 98% (staggered)
arrow_ys <- y_max * c(0.32, 0.19, 0.06)

panel_a <- ggplot() +
  # Green shading: all three tail regions
  geom_area(
    data = shade_car, aes(x = x, y = y, fill = level),
    position = "identity", alpha = 0.6
  ) +
  # Density curve
  geom_line(data = df_car, aes(x = x, y = y), linewidth = 0.6) +
  scale_fill_manual(
    values = GREEN_SHADES,
    labels = paste0(confidence_levels * 100, "% CaR = ",
                    round(car_values), " kg"),
    name = NULL
  ) +
  # Vertical dashed lines at thresholds
  annotate("segment",
           x = car_thresholds, xend = car_thresholds,
           y = 0, yend = y_max * 1.05,
           linetype = "dashed", color = unname(GREEN_SHADES),
           linewidth = 0.4) +
  # Q target line + label
  geom_vline(xintercept = Q_TARGET, linetype = "dotted",
             color = "grey50", linewidth = 0.5) +
  annotate("text", x = Q_TARGET + 40, y = y_max * 0.70,
           label = 'italic(Q)*": contracted"', parse = TRUE, hjust = 0,
           size = 2, color = "grey40") +
  annotate("segment", x = Q_TARGET + 35, xend = Q_TARGET + 3,
           y = y_max * 0.70, yend = y_max * 0.70,
           arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
           color = "grey40", linewidth = 0.3) +
  # mu tick mark at y=0 level
  annotate("segment", x = CAR_MU, xend = CAR_MU,
           y = -y_max * 0.012, yend = y_max * 0.012,
           color = "grey40", linewidth = 0.4) +
  annotate("text", x = CAR_MU, y = -y_max * 0.04,
           label = expression(mu), size = 2, color = "grey40") +
  # x=0 vertical line
  annotate("segment", x = 0, xend = 0,
           y = 0, yend = y_max * 1.05,
           linetype = "dotted", color = "grey50", linewidth = 0.5) +
  # Green "Carbon Removed" arrows (0 -> each threshold)
  annotate("segment",
           x = rep(0, 3), xend = car_thresholds,
           y = arrow_ys, yend = arrow_ys,
           arrow = arrow(length = unit(0.12, "cm")),
           color = unname(GREEN_SHADES), linewidth = 0.7) +
  annotate("text",
           x = car_thresholds / 2,
           y = arrow_ys + y_max * 0.04,
           label = paste0("Carbon Removed (", confidence_levels * 100, "% sure)"),
           size = 1.8, color = unname(GREEN_SHADES), fontface = "bold",
           hjust = 0.5) +
  # Red "CaR" arrows (Q -> each threshold)
  annotate("segment",
           x = rep(Q_TARGET, 3), xend = car_thresholds,
           y = arrow_ys, yend = arrow_ys,
           arrow = arrow(length = unit(0.12, "cm")),
           color = unname(RED_SHADES), linewidth = 0.7) +
  annotate("text",
           x = (car_thresholds + Q_TARGET) / 2,
           y = arrow_ys + y_max * 0.04,
           label = paste0(confidence_levels * 100, "% CaR"),
           size = 1.8, color = unname(RED_SHADES), fontface = "bold",
           hjust = 0.5) +
  # Below-axis: "More failure / Less failure" arrow
  annotate("segment", x = 200, xend = Q_TARGET,
           y = -y_max * 0.15, yend = -y_max * 0.15,
           arrow = arrow(length = unit(0.12, "cm"), ends = "both"),
           color = "grey50", linewidth = 0.5) +
  annotate("text", x = 180, y = -y_max * 0.15,
           label = "More\nfailure",
           size = 1.7, hjust = 1, color = "grey50", fontface = "italic",
           lineheight = 0.85) +
  annotate("text", x = Q_TARGET + 50, y = -y_max * 0.15,
           label = "Less\nfailure",
           size = 1.7, hjust = 0, color = "grey50", fontface = "italic",
           lineheight = 0.85) +
  # Formatting
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_x_continuous(breaks = seq(0, 1000, by = 200)) +
  scale_y_continuous(breaks = 0, labels = "0") +
  coord_cartesian(xlim = c(0, 1300),
                  ylim = c(-y_max * 0.25, y_max * 1.05),
                  clip = "off") +
  labs(
    x = expression("Carbon removed (kg per tonne CO"[2]*"e)"),
    y = "Probability density"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.margin = margin(3, 3, 3, 3),
    axis.title = element_text(size = 6.5),
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 5.5),
    legend.background = element_rect(fill = alpha("white", 0.9),
                                     color = "black", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(override.aes = list(color = "grey40",
                                                  linewidth = 0.3)))

# "Panel (b): Buffer Interpretation" ---------------------------------------
# Recreates the TikZ diagram from cambridge.tex slide 13 in ggplot2.
# Shows: contract bar (Q* + buffer), delivery distribution, p5 = Q*.
# Consistent with buffer_interpretation.tex analysis.

# Buffer interpretation parameters
# Q_star is the target delivery; Q is total contracted
Q_STAR <- p5        # target delivery = 5th percentile of distribution
Q_TOTAL <- Q_TARGET # total contracted (= Q* + CaR_95)
BUFFER  <- Q_TOTAL - Q_STAR  # = CaR_95

# Vertical layout positions (in data coordinates of the delivery axis)
bar_top    <- y_max * 1.50
bar_bottom <- y_max * 1.20
brace_y    <- y_max * 1.65
arrow_mid  <- y_max * 0.85

panel_b <- ggplot() +
  # === Contract bar: Q* portion (green) ===
  annotate("rect",
           xmin = 0, xmax = Q_STAR,
           ymin = bar_bottom, ymax = bar_top,
           fill = GREEN_95, alpha = 0.25, color = "black", linewidth = 0.4) +
  annotate("text", x = Q_STAR / 2, y = (bar_top + bar_bottom) / 2,
           label = 'italic(Q)*"*: target"', parse = TRUE, size = 2.2) +
  # === Contract bar: Buffer portion (red) ===
  annotate("rect",
           xmin = Q_STAR, xmax = Q_TOTAL,
           ymin = bar_bottom, ymax = bar_top,
           fill = RED_95, alpha = 0.25, color = "black", linewidth = 0.4) +
  annotate("text", x = (Q_STAR + Q_TOTAL) / 2,
           y = (bar_top + bar_bottom) / 2,
           label = expression(CaR[95]),
           size = 2.2) +
  # === Brace above: Q = Q* + CaR_95 ===
  # Bracket line just above the bar
  annotate("segment", x = 0, xend = 0,
           y = bar_top + y_max * 0.02, yend = bar_top + y_max * 0.08,
           color = "grey40", linewidth = 0.3) +
  annotate("segment", x = Q_TOTAL, xend = Q_TOTAL,
           y = bar_top + y_max * 0.02, yend = bar_top + y_max * 0.08,
           color = "grey40", linewidth = 0.3) +
  annotate("segment", x = 0, xend = Q_TOTAL,
           y = bar_top + y_max * 0.08, yend = bar_top + y_max * 0.08,
           color = "grey40", linewidth = 0.3) +
  # Text above the bracket line
  annotate("text", x = Q_TOTAL / 2, y = bar_top + y_max * 0.15,
           label = expression(italic(Q) == italic(Q)*"* + " * CaR[95]),
           size = 2.3) +
  # === Arrow down: "some projects fail" ===
  annotate("segment",
           x = Q_TOTAL / 2, xend = Q_TOTAL / 2,
           y = bar_bottom - y_max * 0.02,
           yend = arrow_mid + y_max * 0.05,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           color = "grey50", linewidth = 0.5) +
  annotate("text", x = Q_TOTAL / 2 - 50, y = (bar_bottom + arrow_mid) / 2,
           label = "some projects fail",
           size = 1.8, color = "grey50", hjust = 1, fontface = "italic") +
  # === Delivery distribution ===
  geom_line(data = df_car, aes(x = x, y = y), linewidth = 0.6,
            color = "grey20") +
  # Red tail shading (below p5)
  geom_area(data = shade_tail, aes(x = x, y = y),
            fill = RED_95, alpha = 0.3) +
  # (5% tail label removed for clarity)
  # === Dashed line at p5 = Q* ===
  annotate("segment", x = Q_STAR, xend = Q_STAR,
           y = -y_max * 0.02, yend = bar_bottom,
           linetype = "dashed", color = RED_95, linewidth = 0.5) +
  annotate("text", x = Q_STAR, y = -y_max * 0.05,
           label = expression(italic(p)[5] == italic(Q)*"*"),
           size = 1.8, color = RED_95) +
  # === Dotted line at Q ===
  annotate("segment", x = Q_TOTAL, xend = Q_TOTAL,
           y = -y_max * 0.02, yend = bar_bottom,
           linetype = "dotted", color = "grey50", linewidth = 0.5) +
  # === "Contract" and "Delivery" labels ===
  annotate("text", x = -80, y = (bar_top + bar_bottom) / 2,
           label = "Contract:", size = 2, fontface = "bold",
           hjust = 1, color = "grey30") +
  annotate("text", x = -80, y = y_max * 0.45,
           label = "Delivery:", size = 2, fontface = "bold",
           hjust = 1, color = "grey30") +
  # === CaR decomposition bar below distribution ===
  # Downside risk bracket: p5 to mu
  annotate("segment", x = p5, xend = CAR_MU,
           y = -y_max * 0.26, yend = -y_max * 0.26,
           linewidth = 1.8, color = RED_95) +
  annotate("segment", x = p5, xend = p5,
           y = -y_max * 0.23, yend = -y_max * 0.29,
           linewidth = 0.6, color = RED_95) +
  annotate("segment", x = CAR_MU, xend = CAR_MU,
           y = -y_max * 0.23, yend = -y_max * 0.29,
           linewidth = 0.6, color = RED_95) +
  annotate("text", x = p5, y = -y_max * 0.33,
           label = "Downside risk", size = 1.75, fontface = "bold",
           hjust = 0, color = RED_95) +
  # Expected delivery gap bracket: mu to Q
  annotate("segment", x = CAR_MU, xend = Q_TOTAL,
           y = -y_max * 0.26, yend = -y_max * 0.26,
           linewidth = 1.8, color = DARK_RED) +
  annotate("segment", x = Q_TOTAL, xend = Q_TOTAL,
           y = -y_max * 0.23, yend = -y_max * 0.29,
           linewidth = 0.6, color = DARK_RED) +
  annotate("text", x = CAR_MU, y = -y_max * 0.33,
           label = "Exp. delivery gap", size = 1.75, fontface = "bold",
           hjust = 0, color = DARK_RED) +
  # mu tick mark at y=0 level
  annotate("segment", x = CAR_MU, xend = CAR_MU,
           y = -y_max * 0.02, yend = y_max * 0.02,
           color = "grey40", linewidth = 0.4) +
  annotate("text", x = CAR_MU, y = -y_max * 0.06,
           label = expression(mu), size = 2, color = "grey40") +
  # Q tick mark at y=0 level (separate from dotted line)
  annotate("text", x = Q_TARGET, y = -y_max * 0.06,
           label = expression(italic(Q)), size = 2, color = "grey40") +
  # Formatting
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_x_continuous(breaks = seq(0, 1000, by = 200)) +
  scale_y_continuous(breaks = NULL) +
  coord_cartesian(xlim = c(-120, 1100),
                  ylim = c(-y_max * 0.40, brace_y + y_max * 0.10),
                  clip = "off") +
  labs(
    x = expression("Carbon removed (kg per tonne CO"[2]*"e)"),
    y = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.margin = margin(3, 3, 3, 3),
    axis.title = element_text(size = 6.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

# Assemble & save ----------------------------------------------------------

fig1 <- panel_a + panel_b +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9))

# 183 mm wide = 7.2 inches; height scaled proportionally
ggsave(file.path(out_dir, "figure1.pdf"), fig1,
       width = 183, height = 75, units = "mm")
ggsave(file.path(subfig_dir, "figure1a.pdf"), panel_a,
       width = 91.5, height = 75, units = "mm")
ggsave(file.path(subfig_dir, "figure1b.pdf"), panel_b,
       width = 91.5, height = 75, units = "mm")

cat("Saved figure1.pdf to", out_dir, "\n")
cat("Saved figure1a.pdf, figure1b.pdf to", subfig_dir, "\n")
