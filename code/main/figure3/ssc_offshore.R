# Storage Security Calculator — Well-regulated offshore scenario
# From Alcalde et al. (2018), adapted for standalone use in 
# Original code: code/fig3/DACCS200.R

source("code/main/figure3/ssc_common.R")

set.seed(2100)

# Monte Carlo — OFFSHORE parameter distributions
SSCMC <- function(n = 5) {
  offshore_sampler <- function() {
    list(
      CO2target        = 1.2e10,
      injectperWell    = rnorm(1, 0.75e6, 0.00415e6),
      InjectionPeriod  = 30,
      meanPlumeArea    = rlnorm(1, -0.7595, 0.1763),
      ActiveWellFreq   = rlnorm(1, -2.17, 0.6),
      SlowLeakInjector = rnorm(1, 158.5, 5.2),
      MinorBlowFreq    = runif(1, 0.062, 0.0762),
      MinorBlowout     = exp(rlnorm(1, 1.27, 0.0397)),
      MajorBlowFreq    = rnorm(1, 1.48e-4, 3.33e-5),
      CO2MajorBlowout  = exp(rlnorm(1, 2.57, 0.045)),
      KnownWellDensity = runif(1, 0.4, 0.48),
      wellUnderEst     = 1,
      UnPlugWells      = 0,
      DegradWells      = rlnorm(1, -2.17, 0.6),
      IntactHighRate   = runif(1, 0.049, 0.0595),
      IntactDegrad     = 0,
      PlugBlowoutYear  = rlnorm(1, -8.6125, 0.23),
      BlowoutWellYear  = runif(1, 1e-5, 1e-4),
      CO2largeBlowout  = rlnorm(1, 13.4, 0.35),
      CO2degraded      = runif(1, 270, 330),
      CO2intactHigh    = runif(1, 207, 253),
      CO2intactLow     = runif(1, 0.0036, 0.0044),
      NatLeakRate      = rlnorm(1, 0.693, 0.37),
      A                = rtriangle(1, 3, 53, 12),
      B                = runif(1, 0.0143, 0.5),
      res_sat          = rnorm(1, 0.58, 0.0286),
      migration        = 1
    )
  }
  SSCMC_generic(n, offshore_sampler)
}
