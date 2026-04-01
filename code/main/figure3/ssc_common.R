# Storage Security Calculator — Common functions
# Shared by ssc_offshore.R and ssc_onshore.R
# From Alcalde et al. (2018), adapted for standalone use in 

# 1.1 Triangle distribution (courtesy of R. A. Godfrey)
# https://www.rdocumentation.org/packages/ExtDist/versions/0.6-3/topics/Triangular
rtriangle <- function(n = 1, a = 0, b = 1, c = (a + b) / 2) {
  if (length(n) > 1) n <- length(n)
  if (n < 1 | is.na(n)) stop(paste("invalid argument: n =", n))
  n <- floor(n)
  if (any(is.na(c(a, b, c)))) return(rep(NaN, times = n))
  if (a > c | b < c) return(rep(NaN, times = n))
  if (any(is.infinite(c(a, b, c)))) return(rep(NaN, times = n))
  p <- runif(n)
  if (a != c) {
    i <- which((a + sqrt(p * (b - a) * (c - a))) <= c)
    j <- which((b - sqrt((1 - p) * (b - a) * (b - c))) > c)
  } else {
    i <- which((a + sqrt(p * (b - a) * (c - a))) < c)
    j <- which((b - sqrt((1 - p) * (b - a) * (b - c))) >= c)
  }
  if (length(i) != 0) p[i] <- a + sqrt(p[i] * (b - a) * (c - a))
  if (length(j) != 0) p[j] <- b - sqrt((1 - p[j]) * (b - a) * (b - c))
  return(p)
}

# 1.2 Abandoned well leakage — calculates AB1 and AB2 rates
#' @param KnownWellDensity wells/km2
#' @param wellUnderEst well under-estimation factor
#' @param UnPlugWells proportion of unplugged wells
#' @param DegradWells proportion of degraded wells
#' @param IntactDegrad pre-injection remediation fraction
#' @param IntactHighRate proportion intact with high leak rate
#' @param CO2largeBlowout mass CO2 lost during blowout (t)
#' @param CO2degraded mass CO2 lost per degraded well (t)
#' @param CO2intactHigh mass CO2 lost per intact well, high rate (t)
#' @param CO2intactLow mass CO2 lost per intact well, low rate (t)
#' @param InjectionPeriod years of injection
#' @param BlowoutWellYear long-term blowout rate (events/well/year)
#' @param PlugBlowoutYear short-term blowout rate (wells/injection period)
#' @param ActiveWellsAB1 active wells converted to abandoned post-injection
AbSetUp <- function(KnownWellDensity, wellUnderEst, UnPlugWells, DegradWells,
                    IntactDegrad, IntactHighRate, CO2largeBlowout, CO2degraded,
                    CO2intactHigh, CO2intactLow, InjectionPeriod,
                    BlowoutWellYear, PlugBlowoutYear, ActiveWellsAB1) {
  TrueWellDensity <- KnownWellDensity * wellUnderEst
  UnKnownWellDensity <- TrueWellDensity - KnownWellDensity
  KnownUnPlugWells <- KnownWellDensity * UnPlugWells
  UnKnownUnPlugWells <- UnKnownWellDensity * UnPlugWells
  KnownPlugWells <- KnownWellDensity * (1 - UnPlugWells)
  UnKnownPlugWells <- UnKnownWellDensity * (1 - UnPlugWells)
  TotalUnPlugWells <- KnownUnPlugWells + UnKnownUnPlugWells
  TotalPlugWells <- KnownPlugWells + UnKnownPlugWells
  KnownDegPlug <- KnownPlugWells * DegradWells
  UnKnownDegPlug <- UnKnownPlugWells * DegradWells
  KnownIntactPlug <- KnownPlugWells * (1 - DegradWells)
  UnKnownIntactPlug <- UnKnownPlugWells * (1 - DegradWells)
  # AB1
  RemUnPlugWells <- KnownUnPlugWells
  KnownPlugWells <- RemUnPlugWells + KnownPlugWells
  TotalPlugwells <- KnownPlugWells + UnKnownPlugWells
  RemDegPlug <- KnownDegPlug * IntactDegrad
  KnownIntactWells <- RemDegPlug + RemUnPlugWells + KnownIntactPlug
  KnownDegWells <- KnownDegPlug - RemDegPlug
  UnKnownDegWells <- UnKnownPlugWells * DegradWells
  TotalUnKnownPlug <- UnKnownIntactPlug + UnKnownDegWells
  BlowoutUnPlugged <- UnKnownUnPlugWells / InjectionPeriod
  BlowoutPlugYear <- PlugBlowoutYear / InjectionPeriod
  BlowoutKnownKmYear <- BlowoutPlugYear * KnownWellDensity
  BlowoutUnKnownKmYear <- UnKnownPlugWells * BlowoutPlugYear
  BlowoutTotal <- BlowoutKnownKmYear + BlowoutUnKnownKmYear + BlowoutUnPlugged
  CO2blowout <- BlowoutTotal * CO2largeBlowout
  TotalDegWells <- UnKnownDegWells + KnownDegWells
  TotalIntWells <- UnKnownIntactPlug + KnownIntactWells
  IntactWellsHighRate <- TotalIntWells * IntactHighRate
  IntactWellsLowRate <- TotalIntWells - IntactWellsHighRate
  CO2lostDegrad <- TotalDegWells * CO2degraded
  CO2lostIntactHigh <- IntactWellsHighRate * CO2intactHigh
  CO2lostIntactLow <- IntactWellsLowRate * CO2intactLow
  AB1 <- CO2blowout + CO2lostDegrad + CO2lostIntactHigh + CO2lostIntactLow
  # AB2
  UnKnownBlewOut <- BlowoutUnKnownKmYear * InjectionPeriod
  UnKnownUnRemediated <- UnKnownDegWells - UnKnownBlewOut
  NewTotalWells <- TrueWellDensity + ActiveWellsAB1
  KnownWellsAB2 <- KnownIntactWells + ActiveWellsAB1 + UnKnownUnPlugWells + UnKnownBlewOut + KnownDegWells
  UnKnownIntactHighRate <- UnKnownIntactPlug * IntactHighRate
  UnKnownIntactLowRate <- UnKnownIntactPlug * (1 - IntactHighRate)
  TotalIntactLowRate <- KnownWellsAB2 + UnKnownIntactLowRate
  BlowoutKm2Year <- BlowoutWellYear * NewTotalWells
  CO2lostBlowoutAB2 <- BlowoutKm2Year * CO2largeBlowout
  CO2lostDegrad <- UnKnownUnRemediated * CO2degraded
  CO2lostIntactHigh <- UnKnownIntactHighRate * CO2intactHigh
  CO2lostIntactlow <- TotalIntactLowRate * CO2intactLow
  AB2 <- CO2lostBlowoutAB2 + CO2lostDegrad + CO2lostIntactHigh + CO2lostIntactlow
  output <- c(AB1, AB2)
}

# 2. Base Case — computes CO2 leak/trap matrix over 10,000 years
#' @param CO2target total injection target (tonnes)
#' @param injectperWell tonnes per annum per well
#' @param InjectionPeriod years of injection
#' @param meanPlumeArea km2 per Mt
#' @param ActiveWellFreq slow leak frequency
#' @param SlowLeakInjector slow leak mass (t CO2/yr/well)
#' @param MinorBlowFreq minor blowout frequency (events/well/year)
#' @param MinorBlowout minor blowout mass (t CO2/blowout)
#' @param MajorBlowFreq major blowout frequency
#' @param CO2MajorBlowout major blowout mass (t)
#' @param KnownWellDensity wells/km2
#' @param wellUnderEst well under-estimation factor
#' @param UnPlugWells proportion unplugged
#' @param DegradWells proportion degraded
#' @param IntactHighRate proportion intact with high rate
#' @param IntactDegrad pre-remediation fraction
#' @param PlugBlowoutYear short-term blowout rate
#' @param BlowoutWellYear long-term blowout rate
#' @param CO2largeBlowout blowout mass (t)
#' @param CO2degraded degraded well leak mass (t)
#' @param CO2intactHigh intact well leak mass, high (t)
#' @param CO2intactLow intact well leak mass, low (t)
#' @param NatLeakRate natural leakage rate (t/km2/year)
#' @param A long-term leakage rate
#' @param B exponential decay factor
#' @param res_sat residual saturation fraction
#' @param migration number of volumes of rock CO2 sees
SSCBase <- function(CO2target, injectperWell, InjectionPeriod, meanPlumeArea,
                    ActiveWellFreq, SlowLeakInjector, MinorBlowFreq,
                    MinorBlowout, MajorBlowFreq, CO2MajorBlowout,
                    KnownWellDensity, wellUnderEst, UnPlugWells, DegradWells,
                    IntactHighRate, IntactDegrad, PlugBlowoutYear,
                    BlowoutWellYear, CO2largeBlowout, CO2degraded,
                    CO2intactHigh, CO2intactLow, NatLeakRate, A, B,
                    res_sat, migration) {
  time <- 10000
  annualinject <- CO2target / InjectionPeriod
  wells <- annualinject / injectperWell
  leakAreaNP <- CO2target * meanPlumeArea / 1e6
  leakNatPaths <- leakAreaNP * NatLeakRate
  leakActiveWells <- wells * (SlowLeakInjector * ActiveWellFreq +
    MinorBlowFreq * MinorBlowout + CO2MajorBlowout * MajorBlowFreq)
  ActiveWellsAB1 <- wells / CO2target / meanPlumeArea * 1e6
  AB <- AbSetUp(KnownWellDensity, wellUnderEst, UnPlugWells,
                DegradWells, IntactDegrad, IntactHighRate,
                CO2largeBlowout, CO2degraded, CO2intactHigh, CO2intactLow,
                InjectionPeriod, BlowoutWellYear, PlugBlowoutYear, ActiveWellsAB1)
  leakAbWells1 <- CO2target * meanPlumeArea * AB[1] / 1e6
  leakAbWells2 <- CO2target * meanPlumeArea * AB[2] / 1e6
  leak <- matrix(0, ncol = 10, nrow = time, byrow = F)
  dimnames(leak)[[2]] <- c("time", "injCO2", "leakCO2", "cumleakCO2", "A",
                            "B", "mintrap", "soltrap", "res_trap", "free CO2")
  leak[, 1] <- seq(1, time)
  leak[1:30, 2] <- CO2target * leak[1:30, 1] / 30
  leak[31:time, 2] <- CO2target
  leak[1:30, 3] <- (leakNatPaths + leakActiveWells + leakAbWells1) * leak[1:30, 2] / CO2target
  leak[31:time, 3] <- (A + (100 - A) * exp(-B * (leak[31:time, 1] - 31))) *
    (leakNatPaths + leakAbWells2) / 100
  leak[1, 4] <- 0.5 * leak[1, 3]
  for (i in 2:time) {
    leak[i, 4] <- leak[(i - 1), 4] + 0.5 * (leak[i, 3] + leak[i - 1, 3])
  }
  leak[, 7] <- (leak[, 2] - leak[, 4]) *
    (-1.669497e-13 * leak[, 1]^3 + 2.903893e-9 * leak[, 1]^2 + 1.403357e-5 * leak[, 1])
  leak[, 8] <- (leak[, 2] - leak[, 4]) * 0.204 * (leak[, 1])^0.0342144145
  leak[, 9] <- (leak[, 2] - (leak[, 7] + leak[, 8])) * (1 - (1 - res_sat)^migration)
  leak[, 10] <- leak[, 2] - (leak[, 7] + leak[, 8] + leak[, 4] + leak[, 9])
  leak[, 5] <- A
  leak[, 6] <- B
  output <- leak[leak[, 10] > 0, ]
}

# 3. Generic Monte Carlo loop
#' @param n number of MC realisations
#' @param param_sampler a function() that returns a named list of all
#'   parameters accepted by SSCBase, drawn from the appropriate distributions
#' @return list of 17 output matrices at years 1,3,10,30,100,200,500,
#'   1000,2000,...,10000
SSCMC_generic <- function(n = 5, param_sampler) {
  # Time points to extract
  time_points <- c(1, 3, 10, 30, 100, 200, 500,
                   1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
  n_tp <- length(time_points)

  # Initialise output matrices (one per time point)
  outputs <- lapply(seq_len(n_tp), function(x) matrix(rep(0, 10), ncol = 10, nrow = 1))

  for (i in 1:n) {
    params <- param_sampler()
    newLine <- do.call(SSCBase, params)

    nr <- dim(newLine)[[1]]
    for (k in seq_len(n_tp)) {
      tp <- time_points[k]
      if (nr > tp) {
        outputs[[k]] <- rbind(outputs[[k]], newLine[tp, ])
      } else {
        outputs[[k]] <- rbind(outputs[[k]], newLine[nr, ])
      }
    }
  }

  # Remove the dummy initialisation row from each matrix
  for (k in seq_len(n_tp)) {
    outputs[[k]] <- outputs[[k]][outputs[[k]][, 1] != 0, ]
  }

  return(outputs)
}
