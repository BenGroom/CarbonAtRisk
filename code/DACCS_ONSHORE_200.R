# The Storage Security Calculator (SSC): Calculates CO2 geological storage security over 10,000 years.
# Code associated with the manuscript:
#"Assessing the capability of CO2 Storage Security to deliver on climate mitigation."
#By Juan Alcalde, Stephanie Flude,
#Mark Wilkinson, Gareth Johnson, Katriona Edlmann, Clare E. Bond, 
#Vivian Scott, Stuart M. V. Gilfillan, Xènia Ogaya and R. Stuart Haszeldine.

# Instructions:
# Run the entire code. (We recommend running in RStudio)
# For the base case scenario: 
# Open BasicTable to review the results of the base case scenario
# A graphical display of the base case scenario is automatically generated
# View the key and change displayed parameters in Section 4.3 (Line 565)
# For Monte Carlo uncertainty analysis:
# run: temp <- SSCMC(n) - where n is the required number of realisations
# run FigLoss()
# Results table and CO2 loss curves over time will be displayed
# For the sensitivity analysis results
# If required, change the number of standard deviations being assessed in Section 4.2 (line 526)
# Open MinMaxSA to view results

####################################################################################
########################## Section 1: General Set-up ###############################
####################################################################################
set.seed(2100)

############################## 1.1 Define rtriangle ############################
################## Code to define a triangle distribution courtesy of R. A. Godfrey
################## https://www.rdocumentation.org/packages/ExtDist/versions/0.6-3/topics/Triangular

rtriangle <- function(n=1, a=0, b=1, c=(a+b)/2){
  if(length(n)>1) n <- length(n)
  if(n<1 | is.na(n)) stop(paste("invalid argument: n =", n))
  n <- floor(n)
  if(any(is.na(c(a,b,c)))) return(rep(NaN, times=n)) # to match behavior of runif
  if(a > c | b < c) return(rep(NaN, times=n)) # to match behavior of runif
  if(any(is.infinite(c(a,b,c)))) return(rep(NaN, times=n))
  
  p <- runif(n)
  
  if(a != c){
    # if a = c then i is always true
    i <- which((a + sqrt(p * (b - a)*(c - a))) <= c)
    j <- which((b - sqrt((1 - p) * (b - a) * (b - c))) > c)
  } else {
    i <- which((a + sqrt(p * (b - a)*(c - a))) < c)
    j <- which((b - sqrt((1 - p) * (b - a) * (b - c))) >= c)
  }
  if(length(i)!=0)
    p[i] <- a + sqrt(p[i] * (b - a) * (c - a))
  if(length(j)!=0)
    p[j] <- b - sqrt((1 - p[j]) * (b - a) * (b - c))
  
  return(p)
}

##########################################################################
###################### 1.2 Define AbSetUp function #############################
################This calculates Abandoned Well leakage rates AB1 and AB2

AbSetUp <- function(
    
  KnownWellDensity=2.5,# WELL DENSITY - wells / km2 Known abandoned well density: ONSHORE
  wellUnderEst=1.55, # WELL DENSITY - Well under-estimation factor: ONSHORE (Unregulated)
  UnPlugWells=0.3, # WELL STATUS  - Proportion of unplugged wells
  DegradWells=exp(-2.89), # PLUGGED WELL CONDITIONS - Proportion of Degraded wells: ONSHORE
  IntactDegrad=0, # Pre-injection conversions: Known degraded wells remediated to intact. (Function not used here)
  IntactHighRate=0.054, # SLOW LEAKAGE - Plugged & intact wells with higher leak rate
  CO2largeBlowout=exp(13.4), # Mass of CO2 lost during abandoned well blowout
  CO2degraded=300, # Mass of CO2 lost per leaking well during continuous leakage through degraded wells
  CO2intactHigh=230, # Mass of CO2 lost per well for intact wells with the high leak rate
  CO2intactLow=0.004, # Mass CO2 lost per well for intact wells with the low leak rate
  InjectionPeriod=30, # Number of years of injection 
  BlowoutWellYear=0.00002, # BLOWOUTS  - Long-term blowout rate - events per well year over 10 kyr (NOTE: written as 0.00005 in paper)
  PlugBlowoutYear=exp(-8.6125), # BLOWOUTS - Proportion of plugged wells that blowout during injection period
  ActiveWellsAB1=0.05128205 # The number of active wells that convert to abandoned wells post injection
){
  
  ## Initial Conditions
  
  TrueWellDensity <- KnownWellDensity * wellUnderEst
  UnKnownWellDensity <- TrueWellDensity - KnownWellDensity
  
  KnownUnPlugWells <- KnownWellDensity * UnPlugWells
  UnKnownUnPlugWells <- UnKnownWellDensity * UnPlugWells
  KnownPlugWells <- KnownWellDensity * (1-UnPlugWells)
  UnKnownPlugWells <- UnKnownWellDensity * (1-UnPlugWells)
  
  TotalUnPlugWells <- KnownUnPlugWells + UnKnownUnPlugWells
  TotalPlugWells <- KnownPlugWells + UnKnownPlugWells 
  
  ############ Plugged well conditions
  
  KnownDegPlug <- KnownPlugWells * DegradWells
  UnKnownDegPlug <- UnKnownPlugWells * DegradWells
  KnownIntactPlug <- KnownPlugWells * (1-DegradWells)
  UnKnownIntactPlug <- UnKnownPlugWells * (1-DegradWells)
  
  ######### AB1 ################################
  ## Plugged well status
  
  RemUnPlugWells <- KnownUnPlugWells
  KnownPlugWells <- RemUnPlugWells + KnownPlugWells
  #UnKnownPlugWells 
  TotalPlugwells <- KnownPlugWells + UnKnownPlugWells
  
  #UnKnownUnPlugWells
  
  ###Known plugged well conditions
  
  RemDegPlug <- KnownDegPlug * IntactDegrad
  KnownIntactWells <- RemDegPlug + RemUnPlugWells + KnownIntactPlug
  KnownDegWells <- KnownDegPlug - RemDegPlug
  
  ####Unknown well conditions
  
  #UnKnownUnPlugWells
  UnKnownDegWells <- UnKnownPlugWells * DegradWells
  #UnKnownIntactPlug
  TotalUnKnownPlug <- UnKnownIntactPlug + UnKnownDegWells
  
  ################### BLOWOUTS ##########
  
  BlowoutUnPlugged <- UnKnownUnPlugWells / InjectionPeriod
  BlowoutPlugYear <- PlugBlowoutYear / InjectionPeriod
  BlowoutKnownKmYear <- BlowoutPlugYear * KnownWellDensity
  BlowoutUnKnownKmYear <- UnKnownPlugWells * BlowoutPlugYear
  
  BlowoutTotal <- BlowoutKnownKmYear + BlowoutUnKnownKmYear + BlowoutUnPlugged
  
  CO2blowout <- BlowoutTotal * CO2largeBlowout
  
  ################### SLOW LEAKAGE ######
  
  TotalDegWells <- UnKnownDegWells + KnownDegWells
  
  TotalIntWells <- UnKnownIntactPlug + KnownIntactWells
  IntactWellsHighRate <- TotalIntWells * IntactHighRate
  IntactWellsLowRate <- TotalIntWells - IntactWellsHighRate
  
  CO2lostDegrad <- TotalDegWells * CO2degraded
  CO2lostIntactHigh <- IntactWellsHighRate * CO2intactHigh
  CO2lostIntactLow <- IntactWellsLowRate * CO2intactLow
  
  AB1 <- CO2blowout + CO2lostDegrad + CO2lostIntactHigh + CO2lostIntactLow
  
  ############################### AB2 #####################
  ###Well status at start of AB2 (end AB1)
  
  #UnKnownUnPlugWells
  UnKnownBlewOut <- BlowoutUnKnownKmYear * InjectionPeriod
  #KnownDegWells
  UnKnownUnRemediated <- UnKnownDegWells - UnKnownBlewOut
  
  ###### New well conditions
  NewTotalWells <- TrueWellDensity + ActiveWellsAB1
  
  ######  Known Wells
  KnownWellsAB2 <- KnownIntactWells + ActiveWellsAB1 + UnKnownUnPlugWells + UnKnownBlewOut + KnownDegWells
  
  #######  UnknownWells
  #UnKnownUnRemediated
  #UnKnownIntactPlug
  UnKnownIntactHighRate <- UnKnownIntactPlug * IntactHighRate
  UnKnownIntactLowRate <- UnKnownIntactPlug * (1 - IntactHighRate)
  
  #######  Total Wells by Status
  #UnKnownUnRemediated
  #UnKnownIntactHighRate
  TotalIntactLowRate <- KnownWellsAB2 + UnKnownIntactLowRate
  
  ########### BLOWOUTS
  
  BlowoutKm2Year <- BlowoutWellYear * NewTotalWells
  CO2lostBlowoutAB2 <- BlowoutKm2Year * CO2largeBlowout
  
  ###########  SLOW LEAKAGE
  
  CO2lostDegrad <- UnKnownUnRemediated * CO2degraded
  CO2lostIntactHigh <- UnKnownIntactHighRate * CO2intactHigh
  CO2lostIntactlow <- TotalIntactLowRate * CO2intactLow
  
  AB2 <- CO2lostBlowoutAB2 + CO2lostDegrad + CO2lostIntactHigh + CO2lostIntactlow
  
  
  outputAB1 <- c(TrueWellDensity, UnKnownWellDensity, KnownUnPlugWells, 
                 UnKnownUnPlugWells, KnownPlugWells, UnKnownPlugWells, 
                 TotalUnPlugWells, TotalPlugWells, KnownDegPlug, UnKnownDegPlug, KnownIntactPlug,
                 UnKnownIntactPlug, RemUnPlugWells, KnownPlugWells, UnKnownPlugWells, TotalPlugwells,
                 UnKnownUnPlugWells, RemDegPlug, KnownIntactWells, KnownDegWells, UnKnownUnPlugWells,
                 UnKnownDegWells, UnKnownIntactPlug, TotalUnKnownPlug, BlowoutUnPlugged, BlowoutPlugYear,
                 BlowoutKnownKmYear, BlowoutUnKnownKmYear, BlowoutTotal, CO2blowout, TotalDegWells,
                 TotalIntWells, IntactWellsHighRate, IntactWellsLowRate, CO2lostDegrad, CO2lostIntactHigh,
                 CO2lostIntactLow, AB1)
  
  outputAB2 <- c(UnKnownUnPlugWells, UnKnownBlewOut, KnownDegWells, UnKnownUnRemediated,
                 UnKnownUnRemediated, UnKnownIntactPlug, UnKnownIntactHighRate, UnKnownIntactLowRate, 
                 UnKnownUnRemediated, UnKnownIntactHighRate, TotalIntactLowRate, BlowoutKm2Year,
                 CO2lostBlowoutAB2, CO2lostDegrad, CO2lostIntactHigh, CO2lostIntactlow, AB2)
  
  output <- c(AB1,AB2)
  
}



############################################################################
#######################Section 2 - Base Case Function########################
############################################################################
##This code runs a base case scenario, based on the parameters as defined below.
## It calls AbSetUp and replaces the AbSetUp default parameters with the ones defined here.

SSCBase <- function(CO2target = 1.2e10,  # total injection target, tonnes
                    
                    injectperWell=0.75e6, #####  tpa
                    
                    InjectionPeriod=30, ##### years
                    
                    meanPlumeArea=exp(-0.7595), ####  km2 per Mt
                    
                    ####################  ACTIVE WELL LEAKAGE
                    
                    ActiveWellFreq = exp(-2.89),    # Slow leak frequency, proportion active wells leaking (5) : ONSHORE
                    
                    SlowLeakInjector = 158.5,  # slow leak mass - t CO2 per annum per well 
                    
                    MinorBlowFreq = 0.0693,    # Minor blowout frequency, events / well / year
                    
                    MinorBlowout = exp(exp(1.27)),  # Minor blowout (t CO2 / blowout)
                    
                    MajorBlowFreq = 1.35e-4, # Major blowout frequency (t / event / year) : ONSHORE
                    
                    CO2MajorBlowout = exp(exp(2.57)),  # Major blowout mass (t)  (10)
                    
                    ##################  Abandonned well LEAKAGE
                    
                    KnownWellDensity= 2.5,#### well density (wells / km2) : ONSHORE
                    
                    wellUnderEst=1.55, #### well underestimation factor 
                    
                    UnPlugWells=0.3, #### proportion unplugged wells : ONSHORE
                    
                    DegradWells=exp(-2.89), #### proportion degraded wells
                    
                    IntactHighRate=0.054, #### proportion intact with high rate (15)
                    
                    IntactDegrad=0, ### Pre-remediate degraded? (Extra function not used in this scenario)
                    
                    PlugBlowoutYear=exp(-8.6125), #### Short term blowout risk (wells / 30 years)
                    
                    BlowoutWellYear=2e-5,  #### long trm blowout risk (wells / yr)
                    
                    CO2largeBlowout=exp(13.4), #### Blowout Mass (t)
                    
                    CO2degraded=300, #### Degraded well leak mass (t) (20)
                    
                    CO2intactHigh=230,  #### Intact well leak mass - high (t)
                    
                    CO2intactLow=0.004, #### Intact well leak mass - high (t)
                    
                    #########################
                    
                    NatLeakRate = exp(0.693), #### Natural leakage rate (t / km2 / year)
                    
                    ####################### INPUTS for variable leakage decay curves
                    
                    A = 12,  ##  long-term leakage rate, central case
                    
                    B = 0.25715, ## exponential decay factor, central case (25)
                    
                    ##################
                    
                    res_sat = 0.58,  # Res sat fraction
                    
                    migration = 1 #  no of volumes of rock the CO2 'sees'
                    
){
  
  # has res trapping prop to FREE CO2 ( = injected - (min + sol trapped); not total injected CO2)
  
  ######################### input data
  
  time <- 10000  # total run time for model
  
  annualinject <- CO2target / InjectionPeriod   # Mt
  
  wells <- annualinject/injectperWell
  
  #################### natural pathways
  
  leakAreaNP <- CO2target * meanPlumeArea / 1e6 # km2 (1e6 as leak/area in Mt but target in t)
  
  leakNatPaths <- leakAreaNP * NatLeakRate  #  t/year
  
  #################### Active wells  ##  tonnes / year
  
  leakActiveWells <-  wells * (SlowLeakInjector * ActiveWellFreq + MinorBlowFreq * MinorBlowout + CO2MajorBlowout * MajorBlowFreq)
  
  ##################### Abandoned wells
  
  ## calc ActiveWellsAB1, i.e. the Active wells converted to abandoned at the end of injection
  
  ActiveWellsAB1 <- wells / CO2target / meanPlumeArea * 1e6 
  
  AB <- AbSetUp(KnownWellDensity,wellUnderEst,UnPlugWells, 
                DegradWells, IntactDegrad, IntactHighRate, 
                CO2largeBlowout,CO2degraded,CO2intactHigh,CO2intactLow,
                InjectionPeriod,BlowoutWellYear,PlugBlowoutYear, ActiveWellsAB1)
  
  leakAbWells1 <- CO2target * meanPlumeArea * AB[1] / 1e6  
  leakAbWells2 <- CO2target * meanPlumeArea * AB[2] / 1e6  
  ################# OUTPUT MATRIX
  
  leak <- matrix(0,ncol=10,nrow=time,byrow=F)
  dimnames(leak)[[2]] <- c("time","injCO2","leakCO2","cumleakCO2","A", # 1-5
                           "B","mintrap","soltrap" , "res_trap", "free CO2")
  leak[,1] <- seq(1,time)
  
  ############### CO2 injection
  
  leak[1:30,2] <- CO2target * leak[1:30,1]/30
  leak[31:time,2] <- CO2target
  
  ################# LEAKAGE ###################################################
  
  # rows year <=30
  leak[1:30,3] <- (leakNatPaths + leakActiveWells + leakAbWells1)*leak[1:30,2]/CO2target
  
  # rows year > starts at year 31 as 100 % of the max leakage rate 
  leak[31:time,3] <- (A + (100- A)* exp(-B * (leak[31:time,1]-31)))  * (leakNatPaths + leakAbWells2) / 100
  
  ################## cummulative leakage 
  
  leak[1,4] <- 0.5 * leak[1,3]
  
  for (i in 2:time){
    leak[i,4] <- leak[(i-1),4] + 0.5 * (leak[i,3] + leak[i-1,3])
  }
  
  ################## TRAPPING - mineral and sol first, rest goes to residual trapping
  
  # min trap proportional to total UNLEAKED CO2
  
  # mineral trapping ADJUSTED to START at t=0 (was t=31) otherwise 0 (default in 'leak')
  
  leak[,7] <- (leak[,2] - leak[,4]) * (-1.669497e-13 * leak[,1]^3 + 2.903893e-9 * leak[,1]^2 + 1.403357e-5*leak[,1]) #   
  
  
  # solubility trapping ADJUSTED to START at t=0 (was t=31)
  # sol trap proportional to total UNLEAKED CO2
  
  leak[,8] <-  (leak[,2]-leak[,4]) * 0.204 * (leak[,1])^0.0342144145 ## Xu et al. (2003)
  
  # residual trapping
  
  leak[,9] <- (leak[,2]- (leak[,7] + leak[,8])) * (1 - (1-res_sat)^migration) # leaked not removed!!
  
  #######free CO2
  
  leak[,10] <- leak[,2]- (leak[,7] + leak[,8] + leak[,4] + leak[,9])
  
  leak[,5] <- A # returns the variables from the leakage decay model
  leak[,6] <- B
  
  output <- leak[leak[,10]>0,]  # trims the output if free CO2 < 0
  
  #Run the Basic function (Section 4) for final leakage values.
  # Run the BasicTable code (Section 4) for a plot of CO2 partitioning over time.
  
}

#################################################################
###################Section 3 - Monte Carlo Analysis###########################
#################################################################
##This code carries out a Monte Carlo simulation for n realisations.
## Random numbers are selected for each parameter from the distributions defined below.
## Minimum, P05, P50, P95, and maximum values are extracted for selected years.

SSCMC <- function(n=5){
  
  # Calls SSCBase for range of plume areas, does MC for n times
  
  # Outputs data at years time
  
  # output at 1000 year intervals (output 1 is 1000 years etc)
  output0001 <- matrix(rep(0,10),ncol=10,nrow=1) # 1 year
  output0003 <- matrix(rep(0,10),ncol=10,nrow=1) # 3 year
  output0010 <- matrix(rep(0,10),ncol=10,nrow=1) # 10 year
  output0030 <- matrix(rep(0,10),ncol=10,nrow=1) # 30 years
  output0100 <- matrix(rep(0,10),ncol=10,nrow=1) # 100 years
  output0200 <- matrix(rep(0,10),ncol=10,nrow=1) # 200 years
  output0500 <- matrix(rep(0,10),ncol=10,nrow=1) # 500 years
  output1 <- matrix(rep(0,10),ncol=10,nrow=1) # 1000 years
  output2 <- matrix(rep(0,10),ncol=10,nrow=1) # 2000 years
  output3 <- matrix(rep(0,10),ncol=10,nrow=1) # 3000 years
  output4 <- matrix(rep(0,10),ncol=10,nrow=1) # 4000 years
  output5 <- matrix(rep(0,10),ncol=10,nrow=1) # 5000 years
  output6 <- matrix(rep(0,10),ncol=10,nrow=1) # 6000 years
  output7 <- matrix(rep(0,10),ncol=10,nrow=1) # 7000 years
  output8 <- matrix(rep(0,10),ncol=10,nrow=1) # 8000 years
  output9 <- matrix(rep(0,10),ncol=10,nrow=1) # 9000 years
  output10 <- matrix(rep(0,10),ncol=10,nrow=1) # 10,000 years
  
  for (i in 1:n) {
    
    newLine <- SSCBase(CO2target = 1.2e10,  # total injection target, t, FIXED
                       
                       injectperWell=rnorm(1,0.75e6,0.00415e6), #####  tonnes per annum
                       
                       InjectionPeriod=30, ##### years FIXED
                       
                       meanPlumeArea=rlnorm(1,-0.7595,0.1763), ####  km2 per Mt from S North Sea gas fields
                       
                       ####################  ACTIVE WELL LEAKAGE
                       
                       ActiveWellFreq = rlnorm(1,-2.89,0.7),    # Slow leak frequency, proportion active wells leaking (5)
                       
                       SlowLeakInjector = rnorm(1,158.5,5.2),  # slow leak mass - t CO2 per annum per well 
                       
                       MinorBlowFreq = runif(1,0.062,0.0762), ##### 0.0693 +/- 10% Minor blowout frequency, events / well / year
                       
                       MinorBlowout = exp(rlnorm(1,1.27,0.0397)),  # Minor blowout mass (t CO2 / blowout)
                       
                       MajorBlowFreq = rnorm(1,1.35e-4,4.4e-5), ##### value = 4.96e-4 +/- 10%, Major blowout frequency (t / event / year)
                       
                       CO2MajorBlowout = exp(rlnorm(1,2.57,0.045)),  # (Juan suggested +/- 0.045. for the mean In original code it is 0.011) Major blowout mass (t)  (10) exp(6.6)=750 t as requested 
                       
                       ##################  Abandoned well LEAKAGE
                       
                       KnownWellDensity=runif(1,2.25,2.75), ##### 2.5 +/- 10%, well density (wells / km2)
                       
                       wellUnderEst= runif(1,1.1,2), ####(Recalculated by BG: s.e. for log normal. But original is uniform distribution)  well underestimation factor
                       
                       UnPlugWells=0.3, #### proportion unplugged wells
                       
                       DegradWells=rlnorm(1,-2.89,0.7), # proportion degraded wells
                       
                       IntactHighRate=runif(1,0.049,0.0595),  #####0.054 +/- 10%, proportion intact with high rate (15)
                       
                       IntactDegrad=0, ##### Pre-remediate degraded?
                       
                       PlugBlowoutYear=rlnorm(1,-8.6125,0.23),# 1/2000 to 1/9000, Short term blowout risk (wells / 30 years)
                       
                       BlowoutWellYear=runif(1,1e-5,1e-4), #####1e-4 +/-10% long trm blowout risk (wells / yr)
                       
                       CO2largeBlowout=rlnorm(1,13.4,0.35), # Large Blowout Mass (t)
                       
                       CO2degraded=runif(1,270,330), ####300 +/- 10% Degraded well leak mass (t) (20)
                       
                       CO2intactHigh=runif(1,207,253), ####230 +/- 10% Intact well leak mass - high (t)
                       
                       CO2intactLow=runif(1,0.0036,0.0044), ####0.004 +/- 10%,  Intact well leak mass - high (t)
                       
                       #########################
                       
                       NatLeakRate = rlnorm(1,0.693,0.37), #### Natural leakage rate (t / km2 / year)
                       
                       ####################### INPUTS for variable leakage decay curves
                       
                       A = rtriangle(1,3,53,12),  #  long-term leakage rate, central case
                       
                       B = runif(1,0.0143,0.5), # exponential decay factor, central case (25)
                       
                       ##################
                       
                       res_sat = rnorm(1,0.58,0.0286),  # Res sat fraction
                       
                       migration = 1) #  no of volumes of rock the CO2 'sees' NOT VARIED
    
    ####  take data at year X , take last data point if ran out of leakable CO2 before time
    if (dim(newLine)[[1]]>1)    {output0001 <- rbind(output0001,newLine[1,])}    else {output0001 <- rbind(output0001,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>3)    {output0003 <- rbind(output0003,newLine[3,])}    else {output0003 <- rbind(output0003,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>10)   {output0010 <- rbind(output0010,newLine[10,])}    else {output0010 <- rbind(output0010,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>30)   {output0030 <- rbind(output0030,newLine[30,])}   else {output0030 <- rbind(output0030,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>100)  {output0100 <- rbind(output0100,newLine[100,])}  else {output0100 <- rbind(output0100,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>200)  {output0200 <- rbind(output0200,newLine[200,])}  else {output0200 <- rbind(output0200,newLine[dim(newLine)[[1]],])} # New line for 200 years
    if (dim(newLine)[[1]]>500)  {output0500 <- rbind(output0500,newLine[500,])}  else {output0500 <- rbind(output0500,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>1000) {output1 <- rbind(output1,newLine[1000,])} else {output1 <- rbind(output1,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>2000) {output2 <- rbind(output2,newLine[2000,])} else {output2 <- rbind(output2,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>3000) {output3 <- rbind(output3,newLine[3000,])} else {output3 <- rbind(output3,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>4000) {output4 <- rbind(output4,newLine[4000,])} else {output4 <- rbind(output4,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>5000) {output5 <- rbind(output5,newLine[5000,])} else {output5 <- rbind(output5,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>6000) {output6 <- rbind(output6,newLine[6000,])} else {output6 <- rbind(output6,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>7000) {output7 <- rbind(output7,newLine[7000,])} else {output7 <- rbind(output7,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>8000) {output8 <- rbind(output8,newLine[8000,])} else {output8 <- rbind(output8,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>9000) {output9 <- rbind(output9,newLine[9000,])} else {output9 <- rbind(output9,newLine[dim(newLine)[[1]],])}
    if (dim(newLine)[[1]]>10000) {output10 <- rbind(output10,newLine[10000,])} else {output10 <- rbind(output10,newLine[dim(newLine)[[1]],])}
    
  }
  
  output0001 <- output0001[output0001[,1]!=0,] #  removes the dummy row used to make the matrix 'output'
  output0003 <- output0003[output0003[,1]!=0,]
  output0010 <- output0010[output0010[,1]!=0,]
  output0030 <- output0030[output0030[,1]!=0,]
  output0100 <- output0100[output0100[,1]!=0,]
  output0200 <- output0200[output0200[,1]!=0,]
  output0500 <- output0500[output0500[,1]!=0,]
  output1 <- output1[output1[,1]!=0,]
  output2 <- output2[output2[,1]!=0,]
  output3 <- output3[output3[,1]!=0,]
  output4 <- output4[output4[,1]!=0,]
  output5 <- output5[output5[,1]!=0,]
  output6 <- output6[output6[,1]!=0,]
  output7 <- output7[output7[,1]!=0,]
  output8 <- output8[output8[,1]!=0,]
  output9 <- output9[output9[,1]!=0,]
  output10 <- output10[output10[,1]!=0,]

  output <- list(output0001,output0003,output0010,output0030,output0100,output0200,output0500,output1,output2,output3,output4,output5,output6,output7,output8,output9,output10)
  
}
## Run the FigLoss function (FigLoss(),) to produce output table and plot


#############################################################################
#########################Section 4 - Data Interrogation######################
############################################################################
#Sections of code to be run as required

###### 4.1 Basic function returns the total leakage at yr = 10,000 as a % of injected.#####
Basic <- function(...) { 
  Basic <- SSCBase(...)
  Basic <- Basic[dim(Basic)[[1]],] # takes last line
  Basic <- Basic[4]/Basic[2]*100  # % leak at end of run
  return(Basic)
}


########### 4.2 Code for sensitivity Analysis#################
#This calculates the % leak at t=10,000 for the base case, but varying each parameter according to the min and max pairs below
SDN <- 2 # Set the number of standard deviations to test the sensitivity on
MinMaxSA <- as.data.frame(t(as.matrix(data.frame(
  BC = Basic(), # Base case scenario
  ###### Min and Max pairs below reset the min and max values for each parameter
  Min_injectperWell = Basic(injectperWell=0.75e6-(SDN*0.083e6)), # SDN Standard Deviations
  Max_injectperWell = Basic(injectperWell=0.75e6+(SDN*0.083e6)),
  Min_meanPlumeArea = Basic(meanPlumeArea=exp(-0.7595-(SDN*0.8815))), # SDN Standard Deviations
  Max_meanPlumeArea = Basic(meanPlumeArea=exp(-0.7595+(SDN*0.8815))), #
  Min_ActiveWellFreq = Basic(ActiveWellFreq=exp(-2.89-(SDN*0.7))), # SDN Standard Deviations
  Max_ActiveWellFreq = Basic(ActiveWellFreq=exp(-2.89+(SDN*0.7))), # 
  Min_SlowLeakInjector = Basic(SlowLeakInjector=158.5-(SDN*18.83)), # SDN Standard Deviations
  Max_SlowLeakInjector = Basic(SlowLeakInjector=158.5+(SDN*18.83)), #
  Min_MinorBlowout = Basic(MinorBlowout=exp(exp(1.27-(SDN*0.21)))), # SDN Standard Deviations
  Max_MinorBlowout = Basic(MinorBlowout=exp(exp(1.27+(SDN*0.21)))), #
  Min_MajorBlowFreq = Basic(MajorBlowFreq=1.35e-4-(SDN*4.4e-5)), # SDN Standard Deviations
  Max_MajorBlowFreq = Basic(MajorBlowFreq=1.35e-4+(SDN*4.4e-5)), #
  Min_CO2MajorBlowout = Basic(CO2MajorBlowout=exp(exp(2.57-(SDN*0.045)))), # SDN Standard Deviations
  Max_CO2MajorBlowout = Basic(CO2MajorBlowout=exp(exp(2.57+(SDN*0.045)))), #
  Min_KnownWellDensity = Basic(KnownWellDensity=0) ,
  Max_KnownWellDensity = Basic(KnownWellDensity=5) ,
  Min_wellUnderEst = Basic(wellUnderEst=1.1), # Minimum value
  Max_wellUnderEst = Basic(wellUnderEst=2), # Maximum value
  Min_DegradWells = Basic(DegradWells=exp(-2.89-(SDN*0.7))), # SDN Standard Deviations
  Max_DegradWells = Basic(DegradWells=exp(-2.89+(SDN*0.7))), #
  Min_PlugBlowoutYear = Basic(PlugBlowoutYear=exp(-8.6125-(SDN*0.23))), # SDN Standard Deviations
  Max_PlugBlowoutYear = Basic(PlugBlowoutYear=exp(-8.6125+(SDN*0.23))), #
  Min_BlowoutWellYear = Basic(BlowoutWellYear=1e-5), # Minimum value
  Max_BlowoutWellYear = Basic(BlowoutWellYear=1e-4), # Maximum value
  Min_CO2largeBlowout = Basic(CO2largeBlowout=exp(13.4-(SDN*0.35))), # SDN Standard Deviations
  Max_CO2largeBlowout = Basic(CO2largeBlowout=exp(13.4+(SDN*0.35))), #
  Min_NatLeakRate = Basic(NatLeakRate=exp(0.693-(SDN*0.37))), # SDN Standard Deviations
  Max_NatLeakRate = Basic(NatLeakRate=exp(0.693+(SDN*0.37))),
  Min_A = Basic(A=3), # Minimum value
  Max_A = Basic(A=53), # Maximum value
  Min_B = Basic(B=0.0143), # Minimum value
  Max_B = Basic(B=0.5), #Maximum value
  Min_res_sat = Basic(res_sat=0.58-(SDN*0.1897)), # SDN Standard Deviations
  Max_res_sat = Basic(res_sat=0.58+(SDN*0.1897)) #
))))
## Change the SDN as required. Run Code. Call "MnMaxSA" for output.

#################### 4.3 Plot CO2 partitioning over time ######################
## Code creates a plot of CO2 partitioning over 10 kyrs
BasicTable <- SSCBase()
plot(BasicTable[,"time"],BasicTable[,"injCO2"], type="l", xlab = "time (years)", ylab = "CO2 (t)", log = "x") #, ylim = c(1e5, 2e10),
lines(BasicTable[,"time"],BasicTable[,"cumleakCO2"], col="red") # Plots cumulatively leaked
lines(BasicTable[,"time"],BasicTable[,("mintrap")]+BasicTable[,("soltrap")], col="green") # Plots chemically trapped
lines(BasicTable[,"time"],BasicTable[,"res_trap"], col="blue") #Plots residually trapped
lines(BasicTable[,"time"],BasicTable[,"free CO2"], col="purple", lty=2)
lines(BasicTable[,"time"],BasicTable[,("mintrap")]+BasicTable[,("soltrap")]+BasicTable[,("res_trap")]+BasicTable[,("free CO2")], col="dark red", lty=2) # Plots CO2 remaining in reservoir
lines(BasicTable[,"time"],BasicTable[,("mintrap")]+BasicTable[,("soltrap")]+BasicTable[,("res_trap")], col="orange") # Plots immobilised CO2 
#Comment out lines that aren't required

############ 4.4 Plot the different probabilities of leakage over time#############
####Plots the Max, p95, P59, P5, and Min model results for 15 time steps
####To creat plot, run temp <- SSCMC(n) where n is the required iterations
FigLoss <- function(){
  
  # Plots graph of % CO2 loss quantiles vs time
  # Also plots histogram of % CO2 loss in year 10,000
  
  output <- matrix(rep(0, 17 * 6), ncol = 6, nrow = 17)  # col1 = time # changed to 17 to account for the 200 year scenario
  output[,1] <- c(1, 3, 10, 30, 100, 200, 500, seq(1000, 10000, 1000))
  
  print(output)
  
  
  loss_final <- NULL  # to store % loss at year 10,000
  
  for (i in 1:17) {
    loss <- temp[[i]][,4] / temp[[i]][,2] * 100  # % CO2 loss
    output[i, 2:6] <- quantile(loss, c(0, 0.1, 0.5, 0.95, 1))
    
    if (i == 17) {  # capture % loss for year 10,000
      loss_final <- loss
    }
  }
  
  # Plot quantiles over time
  plot(output[,1], output[,6], type = "l",
       xlim = c(1, 10000), ylim = c(0, output[17,6]),
       main = "black = max, blue = P95, red = P50, green = P5, purple = min",
       xlab = "Time (years)", ylab = "% CO₂ Loss", log = "x")
  
  lines(output[,1], output[,2], col = "purple")  # min
  lines(output[,1], output[,3], col = "blue")    # P95
  lines(output[,1], output[,4], col = "red")     # P50
  lines(output[,1], output[,5], col = "green")   # P5
  
  dimnames(output)[[2]] <- c("time", "min", "P95", "P50", "P05", "max")

  return(output)
}

PlotLossHistogram <- function(i = 17, temp, time_points = c(1, 3, 10, 30, 200, 500, seq(1000, 10000, 1000))) {
  loss <- temp[[i]][,4] / temp[[i]][,2] * 100
  hist(loss_final,
       breaks = 30,
       col = "steelblue",
       border = "white",
       main = paste("Histogram of % CO₂ Loss at Year", time_points[i]),
       xlab = "% CO₂ Loss",
       ylab = "Frequency")
  
  PlotLossHistogram(i = 17, temp = temp)
  
  
}


#######################################
########Temp plot for year 200 ########
######## This is output 6 #############
#######################################

temp <- SSCMC(10000)  # or SSCMC(1000) or however many simulations you want

library(ggplot2)
library(dplyr)

# Extract % loss at year 200
loss_200 <- temp[[6]][, 4] / temp[[6]][, 2] * 100
n_sim <- length(loss_200)

# Histogram bin width
hist_data <- hist(loss_200, breaks = 30, plot = FALSE)
bin_width <- hist_data$breaks[2] - hist_data$breaks[1]

# Density, scaled to frequency, truncated at x >= 0
dens <- density(loss_200)
density_df <- data.frame(
  x = dens$x,
  y = dens$y * bin_width * n_sim
) %>% filter(x >= 0)

# Calculate percentiles
p90 <- quantile(loss_200, 0.90)
p95 <- quantile(loss_200, 0.95)
p98 <- quantile(loss_200, 0.98)

# Build custom legend labels with values
label_90_95 <- paste0("90–95%: ≥ ", round(p90, 2), "%")
label_95_98 <- paste0("95–98%: ≥ ", round(p95, 2), "%")
label_98_100 <- paste0("98–100%: ≥ ", round(p98, 2), "%")

# Shaded areas
shade_90_95 <- subset(density_df, x >= p90 & x < p95) %>% mutate(band = label_90_95)
shade_95_98 <- subset(density_df, x >= p95 & x < p98) %>% mutate(band = label_95_98)
shade_98_100 <- subset(density_df, x >= p98) %>% mutate(band = label_98_100)

shade_df <- bind_rows(shade_90_95, shade_95_98, shade_98_100)
shade_df$band <- factor(shade_df$band, levels = c(label_90_95, label_95_98, label_98_100))

# Max y for truncating CaR lines
half_y <- max(density_df$y) * 0.5

# Final plot
pp <- ggplot(data.frame(loss = loss_200), aes(x = loss)) +
  geom_histogram(aes(y = ..count..), bins = 40, fill = "gray85", color = "white", alpha = 0.6) +
  geom_area(data = shade_df, aes(x = x, y = y, fill = band), alpha = 0.8) +
  geom_line(data = density_df, aes(x = x, y = y), color = "black", size = 0.7) +
  
  # Truncated vertical CaR lines
  annotate("segment", x = p90, xend = p90, y = 0, yend = half_y, color = "lightgreen", linetype = "dashed", size = 1) +
  annotate("segment", x = p95, xend = p95, y = 0, yend = half_y, color = "mediumseagreen", linetype = "dashed", size = 1) +
  annotate("segment", x = p98, xend = p98, y = 0, yend = half_y, color = "darkgreen", linetype = "dashed", size = 1) +
  
  
  annotate("text", x = p90, y = half_y * 1.05, label = "90%", color = "lightgreen", hjust = -0.1, size = 3) +
  annotate("text", x = p95, y = half_y * 1.05, label = "95%", color = "mediumseagreen", hjust = -0.1, size = 3) +
  annotate("text", x = p98, y = half_y * 1.05, label = "98%", color = "darkgreen", hjust = -0.1, size = 3) +
  
  labs(
    x = "% CO₂ Loss",
    y = "Frequency (Number of Simulations)",
    fill = "Carbon at Risk (%)"
  ) +
  scale_fill_manual(
    values = setNames(
      c("lightgreen", "mediumseagreen", "darkgreen"),
      c(label_90_95, label_95_98, label_98_100)
    )) +
  coord_cartesian(xlim = c(0, max(density_df$x)), ylim = c(0, max(density_df$y))) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.title = element_blank(),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

#Print figure pp
pp

ggsave(
  filename = "DACCS_CaR_density_ONSHORE_200yrs.png",
  plot = pp,
  path = "C:/Users/bdg203/Dropbox/Exeter/Offsets/RevisitingRemovals/Nature Paper/DACCS_Figures",
  width = 8,
  height = 6,
  dpi = 300
)



#######################################
########Temp plot for year 1000########
#######################################

temp <- SSCMC(10000)  # or SSCMC(1000) or however many simulations you want

library(ggplot2)
library(dplyr)

# Extract % loss at year 1,000
loss_1000 <- temp[[8]][, 4] / temp[[8]][, 2] * 100
n_sim <- length(loss_1000)

# Histogram bin width
hist_data <- hist(loss_1000, breaks = 30, plot = FALSE)
bin_width <- hist_data$breaks[2] - hist_data$breaks[1]

# Density, scaled to frequency, truncated at x >= 0
dens <- density(loss_1000)
density_df <- data.frame(
  x = dens$x,
  y = dens$y * bin_width * n_sim
) %>% filter(x >= 0)

# Calculate percentiles
p90 <- quantile(loss_1000, 0.90)
p95 <- quantile(loss_1000, 0.95)
p98 <- quantile(loss_1000, 0.98)

# Build custom legend labels with values
label_90_95 <- paste0("90–95%: ≥ ", round(p90, 2), "%")
label_95_98 <- paste0("95–98%: ≥ ", round(p95, 2), "%")
label_98_100 <- paste0("98–100%: ≥ ", round(p98, 2), "%")

# Shaded areas
shade_90_95 <- subset(density_df, x >= p90 & x < p95) %>% mutate(band = label_90_95)
shade_95_98 <- subset(density_df, x >= p95 & x < p98) %>% mutate(band = label_95_98)
shade_98_100 <- subset(density_df, x >= p98) %>% mutate(band = label_98_100)

shade_df <- bind_rows(shade_90_95, shade_95_98, shade_98_100)
shade_df$band <- factor(shade_df$band, levels = c(label_90_95, label_95_98, label_98_100))

# Max y for truncating CaR lines
half_y <- max(density_df$y) * 0.5

# Final plot
pp <- ggplot(data.frame(loss = loss_1000), aes(x = loss)) +
  geom_histogram(aes(y = ..count..), bins = 43, fill = "gray85", color = "white", alpha = 0.6) +
  geom_area(data = shade_df, aes(x = x, y = y, fill = band), alpha = 0.8) +
  geom_line(data = density_df, aes(x = x, y = y), color = "black", size = 0.7) +
  
  # Truncated vertical CaR lines
  annotate("segment", x = p90, xend = p90, y = 0, yend = half_y, color = "lightgreen", linetype = "dashed", size = 1) +
  annotate("segment", x = p95, xend = p95, y = 0, yend = half_y, color = "mediumseagreen", linetype = "dashed", size = 1) +
  annotate("segment", x = p98, xend = p98, y = 0, yend = half_y, color = "darkgreen", linetype = "dashed", size = 1) +
  
  
  annotate("text", x = p90, y = half_y * 1.05, label = "90%", color = "lightgreen", hjust = -0.1, size = 3) +
  annotate("text", x = p95, y = half_y * 1.05, label = "95%", color = "mediumseagreen", hjust = -0.1, size = 3) +
  annotate("text", x = p98, y = half_y * 1.05, label = "98%", color = "darkgreen", hjust = -0.1, size = 3) +
  
  labs(
    x = "% CO₂ Loss",
    y = "Frequency (Number of Simulations)",
    fill = "Carbon at Risk (%)"
  ) +
  scale_fill_manual(
    values = setNames(
      c("lightgreen", "mediumseagreen", "darkgreen"),
      c(label_90_95, label_95_98, label_98_100)
    )) +
  coord_cartesian(xlim = c(0, max(density_df$x)), ylim = c(0, max(density_df$y))) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.title = element_blank(),
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

#Print figure pp
pp

ggsave(
  filename = "DACCS_CaR_density_ONSHORE_1000yrs.png",
  plot = pp,
  path = "C:/Users/bdg203/Dropbox/Exeter/Offsets/RevisitingRemovals/Nature Paper/DACCS_Figures",
  width = 8,
  height = 6,
  dpi = 300
)









