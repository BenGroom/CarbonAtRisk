# SI figure: the CaR decomposition when Q < mu (conservative crediting)
#
# Simplified version of main-text Fig. 1a, showing only the 95% CaR, with the
# contract set below expected delivery. The contracting gap Q - mu is then
# negative and partly offsets the tail component, so CaR is smaller than the
# tail component alone. The two decomposition brackets overlap in x, so they are
# drawn on separate rows.
#
# Run from repo root:
#   Rscript code/si/si_car_negative_gap.R
#
# Output:
#   - outputs/si/si_car_negative_gap.pdf

library(ggplot2)
library(truncnorm)

si_dir <- "outputs/si"
if (!dir.exists(si_dir)) dir.create(si_dir, recursive = TRUE)

# Parameters ---------------------------------------------------------------
# Distribution is identical to main-text Fig. 1; only Q moves, from 1000 to 700.

CAR_MEAN  <- 850
CAR_SD    <- 150
CAR_LOWER <- 0
CAR_UPPER <- Inf
Q_TARGET  <- 700    # contracted volume, deliberately below expected delivery

CAR_MU <- etruncnorm(a = CAR_LOWER, b = CAR_UPPER, mean = CAR_MEAN, sd = CAR_SD)
D5     <- qtruncnorm(0.05, a = CAR_LOWER, b = CAR_UPPER,
                     mean = CAR_MEAN, sd = CAR_SD)

car_95      <- Q_TARGET - D5      # net CaR, positive but small
gap         <- Q_TARGET - CAR_MU  # contracting gap, negative here
tail_comp   <- CAR_MU - D5        # tail component, unchanged in concept

stopifnot(abs((gap + tail_comp) - car_95) < 1e-8)

GREEN_95 <- "#33a02c"
RED_95   <- "#fb6a4a"
DARK_RED <- "#8b0000"
BLUE     <- "#2171B5"

# Distribution -------------------------------------------------------------

x_car <- seq(CAR_LOWER, CAR_MEAN + 3 * CAR_SD, length.out = 1000)
df_car <- data.frame(
  x = x_car,
  y = dtruncnorm(x_car, a = CAR_LOWER, b = CAR_UPPER,
                 mean = CAR_MEAN, sd = CAR_SD)
)
y_max <- max(df_car$y)

shade_tail <- subset(df_car, x <= D5)

# Rows for the decomposition brackets, below the axis
row_tail <- -y_max * 0.20
row_gap  <- -y_max * 0.34

fig <- ggplot() +
  # 5% tail
  geom_area(data = shade_tail, aes(x = x, y = y),
            fill = GREEN_95, alpha = 0.55) +
  geom_line(data = df_car, aes(x = x, y = y), linewidth = 0.6) +

  # Reference verticals
  annotate("segment", x = D5, xend = D5, y = 0, yend = y_max * 1.02,
           linetype = "dashed", color = GREEN_95, linewidth = 0.45) +
  annotate("segment", x = Q_TARGET, xend = Q_TARGET, y = 0, yend = y_max * 1.02,
           linetype = "dotted", color = "grey35", linewidth = 0.5) +
  annotate("segment", x = CAR_MU, xend = CAR_MU, y = 0, yend = y_max * 1.02,
           linetype = "dotted", color = BLUE, linewidth = 0.5) +

  # Top labels
  annotate("text", x = Q_TARGET, y = y_max * 1.09,
           label = 'italic(Q)', parse = TRUE, size = 2.4, color = "grey25") +
  annotate("text", x = CAR_MU, y = y_max * 1.09,
           label = 'mu', parse = TRUE, size = 2.4, color = BLUE) +
  annotate("text", x = D5, y = y_max * 1.09,
           label = 'italic(D)[5]', parse = TRUE, size = 2.4, color = GREEN_95) +

  # The 95% CaR itself, above the axis
  annotate("segment", x = Q_TARGET, xend = D5,
           y = y_max * 0.50, yend = y_max * 0.50,
           arrow = arrow(length = unit(0.09, "cm"), ends = "both",
                         type = "closed"),
           color = RED_95, linewidth = 0.6) +
  annotate("text", x = (Q_TARGET + D5) / 2, y = y_max * 0.585,
           label = sprintf("95%% CaR = %.0f kg", car_95),
           size = 2.2, color = RED_95, fontface = "bold") +

  # Row 1: tail component, D5 -> mu, positive
  annotate("segment", x = D5, xend = CAR_MU, y = row_tail, yend = row_tail,
           linewidth = 1.6, color = RED_95) +
  annotate("segment", x = c(D5, CAR_MU), xend = c(D5, CAR_MU),
           y = row_tail + y_max * 0.03, yend = row_tail - y_max * 0.03,
           linewidth = 0.5, color = RED_95) +
  annotate("text", x = CAR_MU + 30, y = row_tail,
           label = sprintf('"Tail component"~(mu - italic(D)[5]) == "+%.0f kg"',
                           tail_comp),
           parse = TRUE, size = 2, hjust = 0, color = RED_95) +

  # Row 2: contracting gap, mu -> Q, negative. Arrow points left to signal sign.
  annotate("segment", x = CAR_MU, xend = Q_TARGET, y = row_gap, yend = row_gap,
           linewidth = 1.6, color = DARK_RED,
           arrow = arrow(length = unit(0.09, "cm"), type = "closed")) +
  annotate("text", x = CAR_MU + 30, y = row_gap,
           label = sprintf('"Contracting gap"~(italic(Q) - mu) == "%s%.0f kg"',
                           "−", abs(gap)),
           parse = TRUE, size = 2, hjust = 0, color = DARK_RED) +

  # The arithmetic
  annotate("text", x = 20, y = y_max * 0.95,
           label = 'CaR[95] == (italic(Q) - mu) + (mu - italic(D)[5])',
           parse = TRUE, size = 2.2, hjust = 0, color = "grey20") +
  annotate("text", x = 20, y = y_max * 0.845,
           label = sprintf('phantom(CaR[95]) ~ "=  %s%.0f  +  %.0f  =  %.0f kg"',
                           "−", abs(gap), tail_comp, car_95),
           parse = TRUE, size = 2.2, hjust = 0, color = "grey20") +
  annotate("text", x = 20, y = y_max * 0.73,
           label = '"Conservative crediting:"~italic(Q) < mu',
           parse = TRUE, size = 2, hjust = 0, color = "grey45") +

  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_x_continuous(breaks = seq(0, 1200, by = 200)) +
  scale_y_continuous(breaks = 0, labels = "0") +
  coord_cartesian(xlim = c(0, 1420),
                  ylim = c(row_gap - y_max * 0.08, y_max * 1.12),
                  clip = "off") +
  labs(x = expression("Carbon dioxide removed (kg per tonne CO"[2]*"e)"),
       y = "Probability density") +
  theme_classic(base_size = 8) +
  theme(plot.margin = margin(4, 4, 4, 4),
        axis.title = element_text(size = 7))

# Rendered at the width it is displayed at in supplement.tex.
SI_LINEWIDTH_IN <- 500.484 / 72.27
FRAC <- 0.62

ggsave(file.path(si_dir, "si_car_negative_gap.pdf"), fig,
       width = FRAC * SI_LINEWIDTH_IN,
       height = FRAC * SI_LINEWIDTH_IN * 0.72)

cat(sprintf("Q = %.0f, mu = %.1f, D5 = %.1f\n", Q_TARGET, CAR_MU, D5))
cat(sprintf("contracting gap = %.1f, tail component = %.1f, CaR95 = %.1f\n",
            gap, tail_comp, car_95))
cat("Saved", file.path(si_dir, "si_car_negative_gap.pdf"), "\n")
