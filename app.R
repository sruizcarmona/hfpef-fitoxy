library(shiny)

#load libraries
library(ggplot2)
library(plyr)
library(gplots)
library(tidyverse)
library(Bolstad2)
library(bvpSolve) 

#############
# NICK HOUSTIS CODE
## FUNCTIONS
#############

#####################
# DmDlSolver (Algorithm 1 in the O2 pathway paper)
# input: pao2, pvo2, q, hb, pA, pmito
# output: dmo2, dlo2
#####################

model <- function(x,y,parms) {
  q <- parms$q
  hb <- parms$hb
  pA <- parms$pA
  pmito <- parms$pmito
  # model represented as a list
  return(list( c(
    (y[3]/(TT*q*odc(y[1],hb,1)))*(pA-y[1]), # y[1] = Lung o2 
    -(y[4]/(TT*q*odc(y[2],hb,1)))*(y[2]-pmito), #y [2] = Muscle o2
    0, # y[3] = dlo2
    0  # y[4] = dmo2
  )))
}

boundCond <- function(i,y,parms) { 
  pvo2 <- parms$pvo2
  pao2 <- parms$pao2 
  # boundary conditions
  if (i==1) return(y[1]-pvo2) # y[1] = pvo2, at t=0 
  if (i==2) return(y[2]-pao2) # y[2] = pao2, at t=0 
  if (i==3) return(y[1]-pao2) # y[1] = pao2, at t=T 
  if (i==4) return(y[2]-pvo2) # y[2] = pvo2, at t=T 
}

DmDlSolver <- function(...,pao2,pvo2,q,hb,pA,pmito,
                       init_dmo2=1, init_dlo2=1,
                       stepsize=0.01, errtol=1e-3, NITER=10000) {
  
  # error checking
  if (pao2>pA) {
    cat("Impossible measurements: pao2 > pA\n")
    return(data.frame(dmo2=NA,dlo2=NA,pmcap=NA))
  }
  if (pmito>pvo2) {
    cat("Impossible parameters: pmito > pvo2; setting pmito to pvo2-0.05\n")
    pmito <- pvo2-0.05
  }
  
  # initial values  
  xguess = seq(0, TT, by = stepsize)
  yguess = matrix(nrow = 4, ncol = length(xguess), data = 0)
  rownames(yguess) <- c("LungO2", "MuscleO2","dlo2","dmo2")
  yguess[1,] <- (pvo2+(pao2-pvo2)*xguess/TT) # initialize the venous blood gas in the lung
  yguess[2,] <- (pao2-(pao2-pvo2)*xguess/TT) # initialize the arterial blood gas in the muscle
  yguess[3,] <- init_dlo2
  yguess[4,] <- init_dmo2
  
  # solver call
  parms <- list(pao2 = pao2,pvo2 = pvo2,q=q,hb=hb,pA=pA,pmito=pmito)
  Sol <- bvptwp(func = model, 
                bound = boundCond, 
                x = seq(0, TT, by = stepsize), 
                ynames = c("L(x)","M(x)","dlo2","dmo2"),
                parms=parms,
                leftbc=2,
                xguess=xguess,
                yguess=yguess,
                verbose=FALSE,
                atol=errtol,
                nmax=10000)
  
  dmo2 <- Sol[1,"dmo2"]
  dlo2 <- Sol[1,"dlo2"]
  pmcap <- sintegral(Sol[,1],Sol[,"M(x)"])$int / TT # mean value of pmcap(x), computed by integration
  dmdlpcap <- data.frame(dmo2=dmo2,dlo2=dlo2,pmcap=pmcap)
  rownames(dmdlpcap) <- NULL
  return(dmdlpcap)
}

#####################
# bgSolver: blood gas solver, with mitochondrial circuit (Algorithm 2 in the O2 pathway paper)
# input: all the parameters that define the physiology = dmo2, dlo2, q, va, hb, p50, vmax
# output: blood gases = pao2, pvo2, and therefore vo2
#####################

model.bg <- function(x,y,parms) {
  q <- parms$q
  hb <- parms$hb
  dmo2 <- parms$dmo2
  dlo2 <- parms$dlo2
  va <- parms$va
  return(list( c(
    (dlo2/(TT*q*odc(y[1],hb,1)))*((PIO2-(q*(odc(y[5],hb,0)-odc(y[4],hb,0)))/(va*K))-y[1]), # y[1] = Lung o2; rather than use an auxiliary variable for pA I have replaced it with its calculated value
    -(dmo2/(TT*q*odc(y[2],hb,1)))*(y[2]-y[3]), # y[2] = Muscle o2
    0, # y[3] = pmito
    0, # y[4] = pvo2 
    0  # y[5] = pao2 
  )))
}

bound.bg <- function(i,y,parms) { 
  q <- parms$q
  hb<- parms$hb
  va <- parms$va
  p50 <- parms$p50
  vmax <- parms$vmax
  if (i==1) return(y[3]-(p50/(-1+(vmax/(q*(odc(y[5],hb,0)-odc(y[4],hb,0))))))) # pmito formula; y[3] = pmito
  if (i==2) return(y[1]-y[4]) # at t=0, y[1] = y[4], the unknown constant for pvo2
  if (i==3) return(y[2]-y[5]) # at t=0, y[2] = y[5], the unknown constant for pao2
  if (i==4) return(y[1]-y[5]) # at t=T, y[1] = pao2 ie y[5]
  if (i==5) return(y[2]-y[4]) # at t=T, y[2] = pvo2 ie y[4]
}

bgSolver <- function(...,va,q,hb,dmo2,dlo2,p50,vmax,
                     pao2_init=120,pvo2_init=15,pmito_init=5,
                     stepsize=0.01,errtol=1e-3,NITER=10000) {
  
  # initial values
  xguess = seq(0,TT,by=stepsize)
  yguess = matrix(nrow = 5, ncol = length(xguess), data = 0)
  rownames(yguess) <- c("L(x)", "M(x)","pmito","pvo2","pao2")
  yguess[1,] <- (pvo2_init+(pao2_init-pvo2_init)*xguess/TT) # initialize the venous blood gas in the lung
  yguess[2,] <- (pao2_init-(pao2_init-pvo2_init)*xguess/TT) # initialize the arterial blood gas in the muscle
  yguess[3,] <- pmito_init
  yguess[4,] <- pvo2_init
  yguess[5,] <- pao2_init
  
  # solver call
  parms <- list(dmo2 = dmo2,dlo2 = dlo2,va=va,q=q,hb=hb,p50=p50,vmax=vmax)
  
  Sol <- bvptwp(func = model.bg, bound=bound.bg, x = seq(0, TT, by = stepsize), ynames = c("L(x)", "M(x)","pmito","pvo2","pao2"),parms=parms,atol=errtol,leftbc=3,xguess=xguess,yguess=yguess,verbose=FALSE,nmax=NITER)
  
  # solution quantities of interest
  pao2 <- Sol[1,3]
  pvo2 <- Sol[1,2] 
  pmito <- Sol[1,4]
  avo2 <- (odc(pao2,hb,0)-odc(pvo2,hb,0))/10
  vo2 <- q*avo2*10
  pA <- PIO2 - vo2/(va*K)
  bg <- data.frame(pao2.alg2=pao2,pvo2.alg2=pvo2,avo2.alg2=avo2,vo2.alg2=vo2,pA.alg2=pA,pmito.alg2=pmito,vmax.alg2=vmax,q.alg2=q)
  rownames(bg) <- NULL
  return(bg)
}


###############
# O2 dissociation curve and its derivative
# Dash-Bassingthwaighte formulation used here, but Kelman is another popular alternative
###############

odc <- function(x,hb,flag) { #wrapper for odcDB
  if (flag==0) return((0.003*x+1.39*hb*odcDB(x,hb,flag))*10) # o2 content: mL O2/ L blood
  if (flag==1) return(odcDB(x,hb,flag)*hb*1.39*10+0.03) # o2 content DERIVATIVE, ie change in o2 content per delta pao2
}


# Simulation of oxyhemoglobin (HbO2) and carbomino hemoglobin (HbCO2)
# dissociation curves and computation of total O2 and CO2 contents in 
# the whole blood (Dash-Bassingthwaighte, ABME 38(4):1683-1701, 2010)

odcDB <- function(pO2,hb,flag) {
  
  # set pCO2, pH, DPG, Temp to default values unless explicitly modeled
  pCO2 <- 40
  pH <- 7.24
  DPG <- 0.00465
  Temp <- 37
  Hct <- hb
  
  # Parameters those are fixed in the model (i.e., water fractions, hemoglobin
  # concentration in RBCs, equilibrium constants, and Hill coefficient)
  Wpl = 0.94             # fractional water space in plasma unitless
  Wrbc = 0.65            # fractional water space in RBCs unitless
  Rrbc = 0.69            # Gibbs-Donnan ratio across RBC membrane unitless
  Hbrbc = 5.18e-3        # hemoglobin concentration in RBCs M
  K2 = 2.95e-5           # CO2 + HbNH2 equilibrium constant unitless
  K2dp = 1.0e-6          # HbNHCOOH dissociation constant M
  K2p = K2/K2dp          # kf2p/kb2p 1/M
  K3 = 2.51e-5           # CO2 + O2HbNH2 equilibrium constant unitless
  K3dp = 1.0e-6          # O2HbNHCOOH dissociation constant M
  K3p = K3/K3dp          # kf3p/kb3p 1/M
  K5dp = 2.63e-8         # HbNH3+ dissociation constant M
  K6dp = 1.91e-8         # O2HbNH3+ dissociation constant M
  nhill = 2.7            # Hill coefficient unitless
  n0 = nhill-1.0         # deviation of Hill coefficient or cooperativity from
  # the stochiometry of O2 for each heme site
  
  # Variables those are privately fixed in the model with the standard 
  # physiological values (i.e., pH0, pCO20, DPG0, Temp0)
  pO20 = 100.0           # standard O2 partial pressure in blood mmHg
  pCO20 = 40.0           # standard CO2 partial pressure in blood mmHg
  pH0 = 7.24             # standard pH in RBCs unitless
  DPG0 = 4.65e-3         # standard 23-DPG concentration in RBCs M
  Temp0 = 37.0           # standard temperature in blood degC
  fact = 1.0e-6/Wpl      # a multiplicative factor M/mmHg
  alphaO20 = fact*1.37	# solubility of O2 in water at 37 C M/mmHg
  alphaCO20 = fact*30.7	# solubility of CO2 in water at 37 C M/mmHg
  O20 = alphaO20*pO20	# standard O2 concentration in RBCs M
  CO20 = alphaCO20*pCO20	# standard CO2 concentration in RBCs M
  Hp0 = 10^(-pH0)        # standard H+ concentration in RBCs M
  pHpl0 = pH0-log10(Rrbc)	# standard pH in plasma unitless
  P500 = 26.8            # standard pO2 at 50# SHbO2 mmHg
  C500 = alphaO20*P500	# standard O2 concentration at 50# SHbO2 M
  
  # Calculation of intermediate variables and the indices n1, n2, n3, and 
  # n4 in the computations of SHbO2 & SHbCO2
  Wbl = (1-Hct)*Wpl + Hct*Wrbc
  pHpl = pH-log10(Rrbc)
  pHpldiff = pHpl-pHpl0
  pHdiff = pH-pH0
  pCO2diff = pCO2-pCO20
  DPGdiff = DPG-DPG0
  Tempdiff = Temp-Temp0
  alphaO2 = fact*(1.37 - 0.0137*Tempdiff + 0.00058*Tempdiff^2)
  alphaCO2 = fact*(30.7 - 0.57*Tempdiff + 0.02*Tempdiff^2)
  pK1 = 6.091 - 0.0434*pHpldiff + 0.0014*Tempdiff*pHpldiff
  K1 = 10^(-pK1)
  O2 = alphaO2*pO2 ##
  CO2 = alphaCO2*pCO2
  Hp = 10^(-pH) ##
  Hppl = 10^(-pHpl)
  
  Term1 = K2p*(1+K2dp/Hp)
  Term2 = K3p*(1+K3dp/Hp)
  Term3 = (1+Hp/K5dp)
  Term4 = (1+Hp/K6dp)
  Term10 = K2p*(1+K2dp/Hp0)
  Term20 = K3p*(1+K3dp/Hp0)
  Term30 = (1+Hp0/K5dp)
  Term40 = (1+Hp0/K6dp)
  Kratio10 = (Term10*CO20+Term30)/(Term20*CO20+Term40)
  Kratio11 = (Term1*CO20+Term3)/(Term2*CO20+Term4)
  Kratio12 = (Term10*alphaCO20*pCO2+Term30)/(Term20*alphaCO20*pCO2+Term40)
  K4dp = Kratio10*O20^n0/C500^nhill
  K4tp = K4dp/O20^n0
  Kratio20 = Kratio10/K4tp
  Kratio21 = Kratio11/K4tp
  Kratio22 = Kratio12/K4tp
  
  P501 = 26.765 - 21.279*pHdiff + 8.872*pHdiff^2
  P502 = 26.80 + 0.0428*pCO2diff + 3.64e-5*pCO2diff^2
  P503 = 26.78 + 795.633533*DPGdiff - 19660.8947*DPGdiff^2
  P504 = 26.75 + 1.4945*Tempdiff + 0.04335*Tempdiff^2 + 0.0007*Tempdiff^3
  C501 = alphaO20*P501
  C502 = alphaO20*P502
  C503 = alphaO20*P503
  C504 = alphaO2*P504
  
  if (abs(pH-pH0) < 1.0e-6) {
    n1 = 1.0       # can be any arbitrary value
  } else {
    n1 = (log10(Kratio21)-nhill*log10(C501))/(pH-pH0)
  }
  if (abs(pCO2-pCO20) < 1.0e-6) {
    n2 = 1.0       # can be any arbitrary value
  } else {
    n2 = (log10(Kratio22)-nhill*log10(C502))/(log10(CO20)-log10(CO2))
  }
  if (abs(DPG-DPG0) < 1.0e-6) {
    n3 = 1.0        # can be any arbitrary value
  } else {
    n3 = (log10(Kratio20)-nhill*log10(C503))/(log10(DPG0)-log10(DPG))
  }
  if (abs(Temp-Temp0) < 1.0e-6) {
    n4 = 1.0       # can be any arbitrary value
  } else {
    n4 = (log10(Kratio20)-nhill*log10(C504))/(log10(Temp0)-log10(Temp))
  }
  Term5 = (Hp0/Hp)^n1*(CO20/CO2)^n2*(DPG0/DPG)^n3*(Temp0/Temp)^n4
  
  # Calculation of Hill coefficients (KHbO2 and KHbCO2) O2 and CO2 saturations
  # of hemoglobin (SHbO2 and SHbCO2) and O2 and CO2 contents in blood. These 
  # are computed as functions of pO2 and pCO2. Also compute the concentrations
  # of all the components of HbO2 and HbCO2.
  K4p = K4dp*(O2/O20)^n0*Term5
  KHbO2 = K4p*(Term2*CO2+Term4)/(Term1*CO2+Term3)
  KHbO2deriv = n0*K4dp*((O2/O20)^(n0-1))*Term5
  KHbCO2 = (Term1+Term2*K4p*O2)/(Term3+Term4*K4p*O2)
  
  ######
  SHbO2 = KHbO2*O2/(1+KHbO2*O2)
  beta = ((KHbO2*(n0+1))/((1+KHbO2*O2)^2))*(alphaO20) # convert back to mmHg from molarity
  ######
  
  SHbCO2 = KHbCO2*CO2/(1+KHbCO2*CO2)
  O2free = Wbl*alphaO2*pO2
  O2bound = 4*Hct*Hbrbc*SHbO2
  O2total = O2free+O2bound
  O2content = 2225.6*O2total
  CO2free = Wbl*alphaCO2*pCO2
  CO2bicarb = ((1-Hct)*Wpl+Hct*Wrbc*Rrbc)*(K1*alphaCO2*pCO2/Hppl)
  CO2bound = 4*Hct*Hbrbc*SHbCO2
  CO2total = CO2free+CO2bicarb+CO2bound
  CO2content = 2225.6*CO2total
  
  if (flag==1) {
    return(beta)
  } else return(SHbO2)
  
  
  
  ## Some extra variables 
  # 	HbNH2 = Hbrbc/((Term1*CO2+Term3)+ K4p*O2*(Term2*CO2+Term4))
  # 	HbNH3p = HbNH2*Hp/K5dp
  # 	O2HbNH2 = K4p*O2*HbNH2
  # 	O2HbNH3p = O2HbNH2*Hp/K6dp
  # 	HbNHCOOH = K2p*CO2*HbNH2
  # 	HbNHCOOm = K2dp*HbNHCOOH/Hp
  # 	O2HbNHCOOH = K3p*CO2*O2HbNH2
  # 	O2HbNHCOOm = K3dp*O2HbNHCOOH/Hp
  # 	SHbO2kin = (O2HbNH2+O2HbNH3p+O2HbNHCOOH+O2HbNHCOOm)/Hbrbc
  # 	SHbCO2kin = (HbNHCOOH+HbNHCOOm+O2HbNHCOOH+O2HbNHCOOm)/Hbrbc
  # 
  # 	O2freepl1 = Wpl*(1-Hct)*alphaO2*pO2
  # 	O2freepl2 = 2225.6*O2freepl1
  # 	O2freerbc1 = Wrbc*Hct*alphaO2*pO2
  # 	O2freerbc2 = 2225.6*O2freerbc1
  # 	O2boundrbc1 = 4*Hct*Hbrbc*SHbO2
  # 	O2boundrbc2 = 2225.6*O2boundrbc1
  # 	CO2freepl1 = Wpl*(1-Hct)*alphaCO2*pCO2
  # 	CO2freepl2 = 2225.6*CO2freepl1
  # 	CO2freerbc1 = Wrbc*Hct*alphaCO2*pCO2
  # 	CO2freerbc2 = 2225.6*CO2freerbc1
  # 	CO2bicarbpl1 = Wpl*(1-Hct)*K1*alphaCO2*pCO2/Hppl
  # 	CO2bicarbpl2 = 2225.6*CO2bicarbpl1
  # 	CO2bicarbrbc1 = Wrbc*Hct*Rrbc*K1*alphaCO2*pCO2/Hppl
  # 	CO2bicarbrbc2 = 2225.6*CO2bicarbrbc1
  # 	CO2boundrbc1 = 4*Hct*Hbrbc*SHbCO2
  # 	CO2boundrbc2 = 2225.6*CO2boundrbc1
  
  #-----------------------------------------------------------------------------
  # The equations for O2 and CO2 saturations of hemoglobin (SHbO2 and SHbCO2)  
  # are derived by considering the various kinetic reactions involving the
  # binding of O2 and CO2 with hemoglobin in RBCs:
  #
  #            kf1p       K1dp
  # 1. CO2+H2O <--> H2CO3 <--> HCO3- + H+  K1=(kf1p/kb1p)*K1dp
  #            kb1p		K1 = 7.43e-7 M K1dp = 5.5e-4 M
  #
  #              kf2p          K2dp
  # 2. CO2+HbNH2 <--> HbNHCOOH <--> HbNHCOO- + H+  K2=(kf2p/kb2p)*K2dp
  #              kb2p		K2 = 2.95e-5 K2dp = 1.0e-6 M
  #
  #                kf3p            K3dp
  # 3. CO2+O2HbNH2 <--> O2HbNHCOOH <--> O2HbNHCOO- + H+ K3=(kf3p/kb3p)*K3dp
  #                kb3p		K3 = 2.51e-5 K3dp = 1.0e-6 M
  #
  #              kf4p          
  # 4. O2+HbNH2 <--> O2HbNH2  K4p=K4dp*func([O2][H+][CO2][DPG]T)
  #              kb4p		K4dp and K4p are to be determined
  #
  #    func = ([O2]/[O2]s)^n0*([H+]s/[H+])^n1*([CO2]s/[CO2])^n2*
  #           ([DPG]s/[DPG])^n3*(Temps/Temp)^n4
  #
  #           K5dp
  # 5. HbNH3+ <--> HbNH2 + H+  K5 = 2.63e-8 M
  #
  #             K6dp
  # 6. O2HbNH3+ <--> O2HbNH2 + H+  K6 = 1.91e-9 M
  #
  # The association and dissociation rate constants of O2 with hemoglobin is
  # assumed to be dependent on [O2] [H+] [CO2] [DPG] and temperature (Temp)
  # such that the equilibrium constant K4p is proportional to ([O2]/[O2]s)^n0
  # ([H+]s/[H+])^n1 ([CO2]s/[CO2])^n2 ([DPG]s/[DPG])^n3 and (Temps/Temp)^n4.
  # The problem is to estimate the values of the proportionality constant K4dp 
  # and the indices n0 n1 n2 n3 and n4 such that SHbO2 is 50# at pO2 = 26.8 
  # mmHg pH = 7.24 pCO2 = 40 mmHg [DPG] = 4.65 mM and Temp = 37 C in RBCs 
  # and the HbO2 dissociation curve shifts appropriately w.r.t. pH and pCO2.
  #----------------------------------------------------------------------------
} 

library(tidyverse)
library(Bolstad2)
library(bvpSolve) 


########################
# CONSTANTS
########################
FIO2 <- 0.21
PIO2 <- (760-47)*FIO2 # partial pressure of inspired o2, Torr; pio2 = 149.7 for fio2=0.21
K <- 1.159 # ml o2 / ml air / Torr
TT = 1 # transit time set arbitrarily to 1 as the calculations are invariant to its value, but useful conceptually

# Mitochondrial constants
P50REF <- 0.24 # mmHg
VRESERVE <- 1.8 # VO2 Knee extension / VO2 cycle, Vo2 normalized to lean mass of exercise muscle; sets a lower bound on vmax (ratios derived from data in Esposito et al 2010)

########################
# Input measurements: pao2 (mmHg), pvo2 (mmHg), paco2 (mmHg), hb (g/dL), vo2 (mL/min), vco2 (mL/min)
# Output O2 pathway parameters: va (L/min), q (L/min), dmo2 (mL/mmHg/min), dlo2 (mL/mmHg/min), vmax (L/min)
# Output O2 tensions: pA (mmHg), pmito (mmHg), average pmcap (mmHg)
########################

# measurements <- data.frame(pao2=97, paco2 = 40, pvo2=21, hb=14, vo2=1600, vco2=1900) # sample values for a single individual
# o2params <- select(data_params,va,q,hb,p50,vmax) %>% mutate(dmo2=dmdlpcap$dmo2,dlo2=dmdlpcap$dlo2)

calc_alg1 <- function(meas){
  data_params <- 
    mutate(meas,
           o2ct.art = mapply(odc,x=pao2,hb=hb,flag=0)/10, # mL O2/ dL blood
           o2ct.ven = mapply(odc,x=pvo2,hb=hb,flag=0)/10, # mL O2/ dL blood
           va=vco2/(K*paco2), # L/min (BTPS), vco2 in mL/min (STPD)
           avo2 = o2ct.art - o2ct.ven, # mL/dL
           q = 0.1*vo2/avo2, # L/min
           o2deliv = q*o2ct.art*10/1000, # L O2/min
           pA = PIO2-vo2/(va*K), # mmHg; vo2 here in mL O2/min (STPD)
           vmax = VRESERVE*vo2, # mL O2/min
           p50 = P50REF, # mmHg
           pmito = p50/((vmax/vo2) - 1)) #mmHg 
  finalparams <- do.call(DmDlSolver,data_params)
  return(cbind(data_params,finalparams))
}

calc_alg2 <- function(meas){
  finalparams <- do.call(bgSolver,meas)
  return(cbind(meas,finalparams))
}

### read preloaded file to make plots
library(readxl)
library(writexl)
precomp_data <- read_excel('data/Combined data_ESC_abstract_060220120_CALCULATED.xlsx')
my_data <- precomp_data
#########
## COR PLOT
create_cor_plot <- function(pData, xVar, yVar,legx,legy,myxlim=c(0,100),myylim=c(0,100),axis=c("title","xlab","ylab")) {
  model <- lm(formula(paste(yVar,"~",xVar)),data=pData)
  model_cor <- cor(pData[,yVar],pData[,xVar])**2
  coeff=coefficients(model)
  eq = paste0("y = ", round(coeff[2],4), "*x + " ,round(coeff[1],2))
  eq_cor=bquote(R^2 ~"=" ~ .(round(model_cor,2)))
  q_plot <- ggplot(pData, aes_string(x=xVar, y=yVar ))
  q_plot + geom_point(aes(color=group)) +
    geom_smooth(method=lm,se=T,color="black") +
    annotate(geom="text", x=legx, y=legy, label=eq, color="black",hjust=0) +
    annotate(geom="text", x=legx, y=legy-(0.05*legy), label=eq_cor, color="black",hjust=0) +
    xlim(myxlim)+
    ylim(myylim) +
    ggtitle(axis[1]) +
    xlab (axis[2]) +
    ylab (axis[3])
}

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
ui <- fluidPage(
  theme = shinythemes::shinytheme("cosmo"),
  titlePanel("HFpEF Project (Houstis Paper Code)"),
  navbarPage("",
######################################################################################################################################################
############### TAB 1
######################################################################################################################################################
    tabPanel("Algorithm 1",
             sidebarLayout(
               sidebarPanel(
                 helpText("Calculate all parameters for a new patient, given the following measurements"),
                 hr(),
                 actionButton("update", "Update Plots"),
                 actionButton("reset", "Reset"),
                 hr(),
                 helpText(h3("Download results:")),
                 downloadButton("dl", "Download"),
                 hr(),
                 helpText(h3("1) Upload values:")),
                 helpText("Please, be sure the headers are same as the precomputed data (download above to check)"),
                 fileInput("ul", "Upload excel file", multiple = FALSE, accept = NULL, width = NULL),
                 helpText(h3("2) Or define manually:")),
                 actionButton("newpatient", "Add Patient"),
                 textInput("newgroup","Group",value = "NEW",width=180),
                 numericInput("vo2","VO2 (ml/min)",value = 1600,width=180),
                 numericInput("vco2","VCO2 (ml/min)",value = 1900,width=180),
                 numericInput("pao2","PaO2 (mmHg)",value = 97,width=180),
                 numericInput("pvo2","PvO2 (mmHg)",value = 21,width=180),
                 numericInput("hb","Hb (g/dL)",value = 14,width=180),
                 # numericInput("q","Q (L/min)",value = 26,width=180),
                 # numericInput("sato2a","SatO2_a (%)",value = 96,width=180),
                 # numericInput("sato2v","SatO2_v (%)",value = 23,width=180),
                 numericInput("paco2","PaCO2 (mmHg)",value = 40,width=180),
                 # numericInput("pvco2","PvCO2 (mmHg)",value = 50,width=180),
                 # numericInput("pha","pH arterial",value = 7.34,width=180),
                 # numericInput("phv","pH venous",value = 7.21,width=180),
                 width=3),
             mainPanel(
                 fluidRow(column(12,DT::dataTableOutput("inDataExcel"),style = "overflow-x: scroll;")),
                 fluidRow(column(12,DT::dataTableOutput("inData"),style = "overflow-x: scroll;")),
                 fluidRow(column(12,DT::dataTableOutput("alldata"),style = "overflow-x: scroll;")),
                 hr(),
                 fluidRow(column(6,plotOutput("plotqvo2")),
                          column(6,plotOutput("plotvavo2"))),
                 fluidRow(column(6,plotOutput("plotdlvo2")),
                          column(6,plotOutput("plotdmvo2")))
               )
             )
    ),
######################################################################################################################################################
############### TAB 2
######################################################################################################################################################
    tabPanel("Algorithm 2",
             sidebarLayout(
               sidebarPanel(
                 helpText("Given all calculated parameters from Algorithm 1, calculate the \"original\" measurements."),
                 hr(),
                 actionButton("update.2", "Correlation Plots"),
                 actionButton("reset.2", "Reset"),
                 hr(),
                 helpText(h3("Download results:")),
                 downloadButton("dl.2", "Download"),
                 hr(),
                 helpText(h3("1) Upload values:")),
                 helpText("Please, be sure the headers are same as the precomputed data (download in Algorithm 1 tab to check)"),
                 fileInput("ul.2", "Upload excel file", multiple = FALSE, accept = NULL, width = NULL),
                 helpText(h3("2) Or define manually:")),
                 actionButton("newpatient.2", "Add Patient"),
                 textInput("newgroup.2","Group",value = "NEW",width=180),
                 numericInput("va.2","VA (L/min)",value = 40.98,width=180),
                 numericInput("q.2","Q (L/min)",value = 12.84,width=180),
                 numericInput("hb.2","Hb (g/dL)",value = 14,width=180),
                 numericInput("vmax.2","Vmax (L/min)",value = 2880,width=180),
                 numericInput("dmo2.2","DM (mL/mmHg/min)",value = 47.35,width=180),
                 numericInput("dlo2.2","DL (mL/mmHg/min)",value = 22.49,width=180),
                 numericInput("p50.2","p50 (mmHg)",value=0.24, width=180),
                 # numericInput("satao2.2","SatO2_a (%)",value = 96,width=180),
                 # numericInput("satcvo2.2","SatO2_v (%)",value = 23,width=180),
                 width=3),
               mainPanel(
                 fluidRow(column(12,DT::dataTableOutput("inDataExcel.2"),style = "overflow-x: scroll;")),
                 hr(),
                 fluidRow(column(12,DT::dataTableOutput("inData.2"),style = "overflow-x: scroll;")),
                 fluidRow(column(12,DT::dataTableOutput("alldata.2"),style = "overflow-x: scroll;")),
                 hr(),
                 fluidRow(column(6,plotOutput("corplot1")),
                         # column(4,plotOutput("corplot3")),
                         column(6,plotOutput("corplot2"))),
                fluidRow(column(6,plotOutput("corplot4")),
                         column(6,plotOutput("corplot5"))),
                 # fluidRow(column(6,plotOutput("corplot1")),
                 #          column(6,plotOutput("corplot2"))),
                 # fluidRow(column(6,plotOutput("corplot3")),
                 #          column(6,plotOutput("corplot4"))),
                 # fluidRow(column(6,plotOutput("corplot5")),
                 #          column(6,plotOutput("corplot6"))),
                 # fluidRow(column(6,plotOutput("corplot7")),
                 #          column(6,plotOutput("corplot8"))),
               )
               )
             ),
######################################################################################################################################################
############### TAB 3
######################################################################################################################################################
    tabPanel("Patient Simulation",
             sidebarLayout(
               sidebarPanel(width = 3,
                            helpText(h3("Step 1:")),
                            helpText("Create a patient, initial data from a random patient."),
                            helpText("Feel free to edit any of the fields!"),
                            actionButton("init_all", "Initiate Patient"),
                            hr(),
                            helpText(h3("Step 2:")),
                            helpText("Run algorithm 1 to calculate all parameters."),
                            actionButton("calcalg1", "Calculate params (Alg 1)"),
                            hr(),
                            helpText(h3("Step 3:")),
                            helpText("Run algorithm 2 to update initial measurements based on Alg 1 results."),
                            helpText("Change fields in previous table to check their effect."),
                            actionButton("calcalg2", "Update measurements (Alg 2)"),
                            hr(),
                            helpText(h4("Repeat steps above to play with different measurements and calculated parameters.")),
               ),
               
               mainPanel(
                 fluidRow(column(12,DT::dataTableOutput("init_meas"),style = "overflow-x: scroll;")),
                 hr(),
                 fluidRow(column(12,DT::dataTableOutput("calc_param"),style = "overflow-x: scroll;")),
                 hr(),
                 fluidRow(column(12,DT::dataTableOutput("final_data"),style = "overflow-x: scroll;")),
                 )
             ))
  )
)

server <- function(input, output,session) {
#######################
## ALGORITHM 1 (TAB 1)
#######################
  output$inData <- DT::renderDataTable( indata())
  indata <- eventReactive(input$newpatient, {
    if(input$newpatient>0){
      newrow <- isolate(c(input$newpatient,
                          input$vo2,
                          input$vco2,
                          input$pao2,
                          input$pvo2,
                          input$hb,
                          input$q,
                          input$sato2a,
                          input$sato2v,
                          input$paco2,
                          input$pvco2,
                          input$pha,
                          input$phv,
                          input$newgroup))
      newtab <- as.data.frame(matrix(data=as.numeric(newrow),ncol=length(newrow),byrow=T))
      newtab[length(newrow)] <- input$newgroup
      colnames(newtab)<-tolower(c("id",
                                  "VO2",
                                  "VCO2",
                                  "PaO2",
                                  "PvO2",
                                  "Hb",
                                  # "Q",
                                  # "satao2",
                                  # "satcvo2",
                                  "PaCO2",
                                  # "PvCO2",
                                  # "pha",
                                  # "phv",
                                  "group"))
      #check if any of the optional variables is there and remove it otherwise
      # if (is.na(newtab$q)) {newtab <- newtab[,names(newtab) != 'q']}
      # if (is.na(newtab$vo2)) {newtab <- newtab[,names(newtab) != 'vo2']}
      # if (is.na(newtab$vco2)) {newtab <- newtab[,names(newtab) != 'vco2']}
      newtab$id <- paste(as.integer(input$newpatient),input$newgroup,sep="_")
      #append new calculations to old data
      tmpres <- calc_alg1(newtab) %>% mutate_if(is.numeric,round,2)
      my_data <<- plyr::rbind.fill(my_data,tmpres)
      #show new patient
      # newtab
      DT::datatable(newtab,options=list(autoWidth=TRUE,dom='t'))
    }
  }, ignoreNULL = FALSE)
  
  output$inDataExcel <- DT::renderDataTable({
    inFile <- input$ul
    if (is.null(inFile)){
      calcdata <<-NULL
      return(NULL)
    }
    inDataExcel <- read_excel(inFile$datapath)
    colnames(inDataExcel) <- tolower(colnames(inDataExcel))
    colnames(inDataExcel)[colnames(inDataExcel) == 'pa'] <- 'pA'
    
    for(i in 1:nrow(inDataExcel)) {
      row <- inDataExcel[i,]
      # do stuff with row
      #clean and remove dmo2 dlo2 and pmcap
      if ("dmo2" %in% names(row)) {row <- row[,names(row) != 'dmo2']}
      if ("dlo2" %in% names(row)) {row <- row[,names(row) != 'dlo2']}
      if ("pmcap" %in% names(row)) {row <- row[,names(row) != 'pmcap']}
      tmpres <- calc_alg1(row) %>% mutate_if(is.numeric,round,2)
      calcdata <<- rbind(calcdata,tmpres)
    }
    calcdata
  })
  
  #print list of all new patients
  output$alldata <- DT::renderDataTable( df())
  df <- eventReactive(input$newpatient, {
    # my_data[seq(38,dim(my_data)[1]),]
    DT::datatable(my_data[seq(38,dim(my_data)[1]),],options=list(autoWidth=TRUE,dom='tlip'))
  })
  
  #reset my_data
  observeEvent(input$reset, {
    my_data <<- precomp_data
    session$reload()
  })
  ################################## PLOTS
  #update plot1
  output$plotqvo2 <- renderPlot({
    plotqvo2()
  })
  plotqvo2 <- eventReactive(input$update, {
    plot_data <- plyr::rbind.fill(my_data,calcdata)
    plot_data$vo2 <- plot_data$vo2/1000
    pTitle <- expression("Correlation between Q and V"["O"[2]])
    pxLab <- expression("V"["O"[2]]*" (L/min)")
    pyLab <- "Q (L/min)"
    create_cor_plot(plot_data,"vo2","q",0,25,c(0,5),c(0,30),c(pTitle,pxLab,pyLab))
  }, ignoreNULL = FALSE)
  #plot 2
  output$plotvavo2 <- renderPlot({
    plotvavo2()
  })
  plotvavo2 <- eventReactive(input$update, {
    plot_data <- plyr::rbind.fill(my_data,calcdata)
    plot_data$vo2 <- plot_data$vo2/1000
    pTitle <- expression("Correlation between V"["A"]*" and V"["O"[2]])
    pxLab <- expression("V"["O"[2]]*" (L/min)")
    pyLab <- expression("V"["A"]* " (L/min)")
    create_cor_plot(plot_data,"vo2","va",0,130,c(0,5),c(0,150),c(pTitle,pxLab,pyLab))
  },ignoreNULL=F)
  
  #plot 3
  output$plotdlvo2 <- renderPlot({
    plotdlvo2()
  })
  plotdlvo2 <- eventReactive(input$update, {
    plot_data <- plyr::rbind.fill(my_data,calcdata)
    plot_data$vo2 <- plot_data$vo2/1000
    pTitle <- expression("Correlation between D"["L"]*" and V"["O"[2]])
    pxLab <- expression("V"["O"[2]]*" (L/min)")
    pyLab <- expression("D"["L"]* " (mL/min" %.% "mmHg)")
    create_cor_plot(plot_data[c(-7,-2,-3),],"vo2","dlo2",0,40,c(0,5),c(0,50),c(pTitle,pxLab,pyLab))
  },ignoreNULL=F)
  
  #plot 4
  output$plotdmvo2 <- renderPlot({
    plotdmvo2()
  })
  plotdmvo2 <- eventReactive(input$update, {
    plot_data <- plyr::rbind.fill(my_data,calcdata)
    plot_data$vo2 <- plot_data$vo2/1000
    pTitle <- expression("Correlation between D"["M"]*" and V"["O"[2]])
    pxLab <- expression("V"["O"[2]]*" (L/min)")
    pyLab <- expression("D"["M"]* " (mL/min" %.% "mmHg)")
    create_cor_plot(plot_data,"vo2","dmo2",0,100,c(0,5),c(0,120),c(pTitle,pxLab,pyLab))
  },ignoreNULL=F)
  
  output$dl <- downloadHandler(
    filename = function() { "outputfile.xlsx"},
    content = function(file) {write_xlsx(plyr::rbind.fill(my_data,calcdata), path = file)}
  )
  ##########################################################################################################################################
  ## ALGORITHM 2 (TAB 2)
  ##########################################################################################################################################
  resdata <- NULL
  output$inData.2 <- DT::renderDataTable( indata2())
  indata2 <- eventReactive(input$newpatient.2, {
    if(input$newpatient.2>0){
      newrow <- isolate(c(input$newpatient.2, input$va.2, input$q.2, input$hb.2,input$vmax.2,
                          input$dmo2.2, input$dlo2.2,input$p50.2,input$satao2.2,input$satcvo2.2,
                          input$newgroup.2))
      newtab <- as.data.frame(matrix(data=as.numeric(newrow),ncol=length(newrow),byrow=T))
      newtab[length(newrow)] <- input$newgroup.2
      colnames(newtab)<-tolower(c("id",
                                  "VA",
                                  "q",
                                  "hb",
                                  "vmax",
                                  "dmo2",
                                  "dlo2",
                                  "p50",
                                  # "satao2",
                                  # "satcvo2",
                                  "group"))
      newtab$id <- paste(as.integer(input$newpatient.2),input$newgroup.2,sep="_")
      # o2params <- select(newtab,va,q,hb,p50,vmax,dmo2,dlo2)
      #show new patient
      tmpres <- calc_alg2(newtab) %>% mutate_if(is.numeric,round,2)
      resdata <<- plyr::rbind.fill(resdata,tmpres)
      DT::datatable(newtab,options=list(autoWidth=TRUE,dom='t'))
    }
  }, ignoreNULL = FALSE)
  
  output$inDataExcel.2 <- DT::renderDataTable({
    inFile.2 <- input$ul.2
    if (is.null(inFile.2)){
      calcdata.2 <<-NULL
      return(NULL)
    }
    inDataExcel.2 <- read_excel(inFile.2$datapath)
    colnames(inDataExcel.2) <- tolower(colnames(inDataExcel.2))
    for(i in 1:nrow(inDataExcel.2)) {
      row <- inDataExcel.2[i,]
      # o2params <- select(row,va,q,hb,p50,vmax,dmo2,dlo2)
      tmpres <- calc_alg2(row) %>% mutate_if(is.numeric,round,2)
      calcdata.2 <<- rbind(calcdata.2,tmpres)
    }
    calcdata.2
  })
  
  #print list of all new patients
  output$alldata.2 <- DT::renderDataTable( df2())
  df2 <- eventReactive(input$newpatient.2, {
    # my_data[seq(38,dim(my_data)[1]),]
    # resdata
    # DT::datatable(my_data[seq(38,dim(my_data)[1]),],options=list(autoWidth=TRUE,dom='tlip'))
    DT::datatable(resdata,options=list(autoWidth=TRUE,dom='tlip'))
  })
  
  #reset my_data
  observeEvent(input$reset.2, {
    resdata <<- NULL
    session$reload()
  })
  
  output$dl.2 <- downloadHandler(
    filename = function() { "outputfile.xlsx"},
    content = function(file) {write_xlsx(plyr::rbind.fill(resdata,calcdata.2), path = file)}
  )
  
  ################# PLOTS
  ###### PLOT 1
  output$corplot1 <- renderPlot({
    cplot1()
  })
  cplot1 <- eventReactive(input$update.2, {
    plot_data <- calcdata.2
    plot_data$vo2 <- plot_data$vo2/1000
    plot_data$va <- plot_data$va/1000
    vartoplot <- names(plot_data)[str_detect(names(plot_data),"alg2")]
    pyLab <- vartoplot[1]
    pxLab <- str_remove(vartoplot[1],".alg2")
    pTitle <- paste0(pxLab," - Experimental vs Algorithm 2 correlation")
    lim <- 200
    create_cor_plot(plot_data,pxLab,pyLab,0,0.8*lim,c(0,lim),c(0,lim),c(pTitle,pxLab,pyLab))
  })
  
  ###### PLOT 2
  output$corplot2 <- renderPlot({
    cplot2()
  })
  cplot2 <- eventReactive(input$update.2, {
    plot_data <- calcdata.2
    plot_data$vo2 <- plot_data$vo2/1000
    plot_data$va <- plot_data$va/1000
    vartoplot <- names(plot_data)[str_detect(names(plot_data),"alg2")]
    pyLab <- vartoplot[2]
    pxLab <- str_remove(vartoplot[2],".alg2")
    pTitle <- paste0(pxLab," - Experimental vs Algorithm 2 correlation")
    lim <- 40
    create_cor_plot(plot_data,pxLab,pyLab,0,0.8*lim,c(0,lim),c(0,lim),c(pTitle,pxLab,pyLab))
  })
  
  # ###### PLOT 3
  # output$corplot3 <- renderPlot({
  #   cplot3()
  # })
  # cplot3 <- eventReactive(input$update.2, {
  #   plot_data <- calcdata.2
  #   plot_data$vo2 <- plot_data$vo2/1000
  #   plot_data$va <- plot_data$va/1000
  #   vartoplot <- names(plot_data)[str_detect(names(plot_data),"alg2")]
  #   pyLab <- vartoplot[3]
  #   pxLab <- str_remove(vartoplot[3],".alg2")
  #   pTitle <- paste0(pxLab," - Experimental vs Algorithm 2 correlation")
  #   lim <- 25
  #   create_cor_plot(plot_data,pxLab,pyLab,0,0.8*lim,c(0,lim),c(0,lim),c(pTitle,pxLab,pyLab))
  # })
  
  ###### PLOT 4
  output$corplot4 <- renderPlot({
    cplot4()
  })
  cplot4 <- eventReactive(input$update.2, {
    plot_data <- calcdata.2
    plot_data$vo2 <- plot_data$vo2/1000
    plot_data$vo2.alg2 <- plot_data$vo2.alg2/1000
    vartoplot <- names(plot_data)[str_detect(names(plot_data),"alg2")]
    pyLab <- vartoplot[4]
    pxLab <- str_remove(vartoplot[4],".alg2")
    pTitle <- paste0(pxLab," - Experimental vs Algorithm 2 correlation")
    lim <- 5
    create_cor_plot(plot_data,pxLab,pyLab,0,0.8*lim,c(0,lim),c(0,lim),c(pTitle,pxLab,pyLab))
  })
  
  ###### PLOT 5
  output$corplot5 <- renderPlot({
    cplot5()
  })
  cplot5 <- eventReactive(input$update.2, {
    plot_data <- calcdata.2
    plot_data$vo2 <- plot_data$vo2/1000
    plot_data$va <- plot_data$va/1000
    vartoplot <- names(plot_data)[str_detect(names(plot_data),"alg2")]
    pyLab <- vartoplot[5]
    pxLab <- str_remove(vartoplot[5],".alg2")
    pTitle <- paste0(pxLab," - Experimental vs Algorithm 2 correlation")
    lim <- 150
    create_cor_plot(plot_data,pxLab,pyLab,0,0.8*lim,c(0,lim),c(0,lim),c(pTitle,pxLab,pyLab))
  })
  
  # ###### PLOT 6
  # output$corplot6 <- renderPlot({
  #   cplot6()
  # })
  # cplot6 <- eventReactive(input$update.2, {
  #   plot_data <- calcdata.2
  #   plot_data$vo2 <- plot_data$vo2/1000
  #   plot_data$va <- plot_data$va/1000
  #   vartoplot <- names(plot_data)[str_detect(names(plot_data),"alg2")]
  #   pTitle <- "Correlation Plot"
  #   pyLab <- vartoplot[6]
  #   pxLab <- str_remove(vartoplot[6],".alg2")
  #   lim <- 0.5
  #   create_cor_plot(plot_data,pxLab,pyLab,0,0.8*lim,c(0,lim),c(0,lim),c(pTitle,pxLab,pyLab))
  # })
  
  # ###### PLOT 7
  # output$corplot7 <- renderPlot({
  #   cplot7()
  # })
  # cplot7 <- eventReactive(input$update.2, {
  #   plot_data <- calcdata.2
  #   plot_data$vo2 <- plot_data$vo2/1000
  #   plot_data$va <- plot_data$va/1000
  #   vartoplot <- names(plot_data)[str_detect(names(plot_data),"alg2")]
  #   pTitle <- "Correlation Plot"
  #   pyLab <- vartoplot[7]
  #   pxLab <- str_remove(vartoplot[7],".alg2")
  #   lim <- 8
  #   create_cor_plot(plot_data,pxLab,pyLab,0,0.8*lim,c(0,lim),c(0,lim),c(pTitle,pxLab,pyLab))
  # })
  
  # ###### PLOT 8
  # output$corplot8 <- renderPlot({
  #   cplot8()
  # })
  # cplot8 <- eventReactive(input$update.2, {
  #   plot_data <- calcdata.2
  #   plot_data$vo2 <- plot_data$vo2/1000
  #   plot_data$va <- plot_data$va/1000
  #   vartoplot <- names(plot_data)[str_detect(names(plot_data),"alg2")]
  #   pTitle <- "Correlation Plot"
  #   pyLab <- vartoplot[8]
  #   pxLab <- str_remove(vartoplot[8],".alg2")
  #   lim <- 30
  #   create_cor_plot(plot_data,pxLab,pyLab,0,0.8*lim,c(0,lim),c(0,lim),c(pTitle,pxLab,pyLab))
  # })
  
  
  
  
  
  ##########################################################################################################################################
  ## Patient Simulation (TAB 3)
  ##########################################################################################################################################
  # table 1 - INPUT PATIENT
  output$init_meas <- DT::renderDataTable( df1tab3())
  df1tab3 <- eventReactive(input$init_all, {
    # patient <<- my_data[1,3:14]
    patient <<- data.frame(pao2=97, paco2 = 40, pvo2=21, hb=14, vo2=1600, vco2=1900) # sample values for a single individual
    DT::datatable(patient,editable=T,rownames=F,options=list(dom='t'))
  })
  
  # table 2 - CALC PARAMS
  output$calc_param <- DT::renderDataTable( df2tab3())
  df2tab3 <- eventReactive(input$calcalg1, {
    # if (is.na(patient$q)) {patient <- patient[,names(patient) != 'q']}
    # if (is.na(patient$vo2)) {patient <- patient[,names(patient) != 'vo2']}
    # if (is.na(patient$vco2)) {patient <- patient[,names(patient) != 'vco2']}
    patres <<- calc_alg1(patient) %>% mutate_if(is.numeric,round,2)
    patrestab <<- select(patres,q,avo2,hb,va,o2deliv,pA,vmax,p50,pmito,dmo2,dlo2,pmcap)
    # patres <<- patres[,12:23]
    DT::datatable(patrestab,editable=T,rownames=F,options=list(dom='t'))
  },ignoreNULL = TRUE)
  
  #table 3 - RECALC INPUT
  output$final_data <- DT::renderDataTable( df3tab3())
  df3tab3 <- eventReactive(input$calcalg2, {
    alg2tmp <- select(patrestab,va,q,hb,vmax,p50,dmo2,dlo2)
    alg2res <- calc_alg2(alg2tmp) %>% mutate_if(is.numeric,round,2)
    alg2res <- select(alg2res,pao2.alg2,pvo2.alg2,avo2.alg2,vo2.alg2,pA.alg2,pmito.alg2,vmax.alg2,q.alg2)
    DT::datatable(alg2res,rownames=F,options=list(dom='t'))
  },ignoreNULL = TRUE)
  
  #table 4 - PLUS/MINUS changed things!!
  
  
  #observeevents to take into account updated values
  observeEvent(input$init_meas_cell_edit,{
    patient[1,input$init_meas_cell_edit$col+1] <<- as.numeric(input$init_meas_cell_edit$value)
  })
  observeEvent(input$calc_param_cell_edit,{
    patrestab[1,input$calc_param_cell_edit$col+1] <<- as.numeric(input$calc_param_cell_edit$value)
  })
  
}
shinyApp(ui = ui, server = server)
