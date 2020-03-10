library(shiny)

#load libraries
# library(bvpSolve)
library(ggplot2)
library(plyr)
library(gplots)
library(tidyverse)
library(Bolstad2)
library(bvpSolve) 

#load input data
my_data <-  c("4632","3655","35","80","21","14.1","25.99","96","23.3","3776","3061","42","97","26","13.3","16.85","98.3","41.7","3661","3048","32","111","26","14.4","21.29","98.7","39.2","5383","4390","34","91","22","15.4","20.62","98.1","29.9","2619","2093","37","89","33","15.9","11.73","97.8","66.6","3617","3039","39","64","27","16.5","17.28","93.1","43.4","4287","5239","35","125","51","16.3","24.46","98.4","35","3479","4631","39","99","27","14.5","26.89","97.8","41","3968","4284","38","104","30","14.1","24.31","98.1","38.8","4745","5857","37","84","24","15.4","28.65","97","28.7","3364","3956","35","86","25","15.4","21.22","97.2","34.7","2091","2767","39","82","26","14.6","12.77","96.6","37.6","3544","4797","39","85","25","15.3","29.24","97.4","40.1","3733","4480","36","105","38","15.7","24.71","98.1","55.8","1338","1456","29","53","21","18.2","8.44","91.2","40","1677","1395","30","67","25","15.3","10.46","92.6","33.1","800","755","28","90","16","13.3","5.48","98","17.3","1128","1130","34","61","26","17.1","9.24","89.2","37.6","1370","1365","23","43","21","15.8","13.55","86.8","45.3","816","859","23","51","19","15.1","5.36","91","35.5","1332","1129","21","60","18","15.7","10.8","91.8","22.2")
my_data <- as.data.frame(matrix(data=as.numeric(my_data),ncol=9,byrow=T))
rownames(my_data) <-c(paste("C",seq(14),sep=""),paste("CTEPH",seq(7),sep=""))
colnames(my_data) <- c("VO2","VCO2","PaCO2","PaO2","PvO2","Hb","Q_exp","satO2_a","satO2_v")

#############
# NICK HOUSTIS CODE
## FUNCTIONS
#############

#####################
# DmDlSolver
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
    (y[3]/(TT*q*o2ct(y[1],hb,1)))*(pA-y[1]), # y[1] = Lung o2 
    -(y[4]/(TT*q*o2ct(y[2],hb,1)))*(y[2]-pmito), #y [2] = Muscle o2
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
  
  return(data.frame(dmo2=dmo2,dlo2=dlo2,pmcap=pmcap))
}


###############
# O2 dissociation curve and its derivative
# Dash-Bassingthwaighte formulation used here, but Kelman is another popular alternative
###############

o2ct <- function(x,hb,flag) { #wrapper for odcDB
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
} 

#####
# PATHWAY.R
# CALCULAIONS
# library(tidyverse)
# library(Bolstad2)
# library(bvpSolve) 


########################
# CONSTANTS
########################
FIO2 <- 0.21
PIO2 <- (760-47)*FIO2 # partial pressure of inspired o2, Torr; pio2 = 149.7 for fio2=0.21
K <- 1.159 # ml o2 / ml air / Torr
TT = 1 # transit time set arbitrarily to 1 as the calculations are invariant to its value, but useful conceptually

# Mitochondrial constants
P50REF <- 0.24 # mmHg
VRESERVE <- 1.8 # VO2 Knee extension / VO2 cycle; sets a lower bound on vmax (derived from Esposito et al 2010)

########################
# Input measurements: pao2 (mmHg), pvo2 (mmHg), paco2 (mmHg), hb (g/dL), vo2 (mL/min), vco2 (mL/min)
# Output O2 pathway parameters: va (L/min), q (L/min), dmo2 (mL/mmHg/min), dlo2 (mL/mmHg/min), vmax (L/min)
# Output O2 tensions: pA (mmHg), pmito (mmHg), average pmcap (mmHg)
########################
# SRC FUNCTIONS TO CALCULATE ALL
calc_params <- function(meas) {
  dp <- meas
  dp <- mutate(dp,
               o2ct.art = 0.0032*pao2+1.4*hb*satao2/100, #SRC
               o2ct.ven = 0.0032*pvo2+1.4*hb*satcvo2/100, #SRC
               avo2 = o2ct.art - o2ct.ven) # mL/dL
  if (!has_name(dp,"q")) {
    if (!has_name(dp,"vo2")) {
      cat("Impossible to calculate without Q or VO2\n")
      return()
    }
    dp <- mutate(dp,
                 q = 0.1*vo2/avo2)
  }
  if (!has_name(dp,"vo2")) {
    dp <- mutate(dp,
                 vo2 = q*avo2/0.1)
  }
  if (!has_name(dp,"vco2")){
      #check missing variables
      if (is.null(dp$pha) || is.null(dp$phv) || is.null(dp$paco2) || is.null(dp$pvco2) || is.null(dp$satao2) || is.null(dp$satcvo2)) {
        cat("Impossible to calculate VCO2, missing values! Please check your input\n")
        return()
      }
      #define vars
      plasmatemp <- 37
      pha <- meas$pha
      phv <- meas$phv
      paco2 <- meas$paco2
      pvco2 <- meas$pvco2
      satao2 <- meas$satao2
      satcvo2 <- meas$satcvo2
      hb <- meas$hb
      q <- meas$q
      #calculate art and ven co2 sol
      co2.s <- 0.0307+(0.00057*(37-plasmatemp))+(0.00002*(37-plasmatemp)^2)
      #calculate apparent pk, pkprime, art and ven
      co2.pkp.art <- 6.086+(0.042*(7.4-pha))+((38-plasmatemp)*(0.00472+0.00139*(7.4-pha)))
      co2.pkp.ven <- 6.086+(0.042*(7.4-phv))+((38-plasmatemp)*(0.00472+0.00139*(7.4-phv)))
      #plasma co2 content, art and ven
      co2.plasma.art <- 2.226*co2.s*paco2*(1+10^(pha-co2.pkp.art))
      co2.plasma.ven <- 2.226*co2.s*pvco2*(1+10^(pha-co2.pkp.ven))
      #blood co2 content, art and ven
      co2ct.art <- co2.plasma.art*(1-(0.0289*hb)/((3.352-0.456*satao2)*(8.142-pha)))
      co2ct.ven <- co2.plasma.ven*(1-(0.0289*hb)/((3.352-0.456*satcvo2)*(8.142-phv)))
      dp <- mutate(dp,
                    co2ct.art=co2ct.art,
                    co2ct.ven=co2ct.ven,
                    vco2=10*q*(co2ct.ven-co2ct.art))
  }
  dp <- mutate(dp,
               # o2ct.art = 0.0032*pao2+1.4*hb*satao2/100, #SRC
               # o2ct.ven = 0.0032*pvo2+1.4*hb*satcvo2/100, #SRC
               # o2ct.art = mapply(o2ct,x=pao2,hb=hb,flag=0)/10, # mL O2/ dL blood
               # o2ct.ven = mapply(o2ct,x=pvo2,hb=hb,flag=0)/10, # mL O2/ dL blood
               # avo2 = o2ct.art - o2ct.ven, # mL/dL
               # q = 0.1*vo2/avo2,q, # L/min
               va=vco2/(K*paco2), # L/min (BTPS), vco2 in mL/min (STPD)
               o2deliv = q*o2ct.art*10/1000, # L O2/min
               pA = PIO2-vo2/(va*K), # vo2 in mL/min (STPD)
               vmax = VRESERVE*vo2, # mL O2/min
               p50 = P50REF, # mmHg
               pmito = p50/((vmax/vo2) - 1)) #mmHg
  return (dp)
}
calc_all <- function(meas){
  dataparams <- calc_params(meas)
  finalparams <- do.call(DmDlSolver,dataparams)
  return(cbind(dataparams,finalparams))
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
  titlePanel("HFpEF Project"),
  navbarPage("",
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
                 textInput("newgroup","Group",value = "NEW",width=120),
                 numericInput("vo2","VO2 (ml/min)",value = 4500,width=120),
                 numericInput("vco2","VCO2 (ml/min)",value = 3600,width=120),
                 numericInput("pao2","PaO2 (mmHg)",value = 80,width=120),
                 numericInput("pvo2","PvO2 (mmHg)",value = 20,width=120),
                 numericInput("hb","Hb (g/dL)",value = 14,width=120),
                 numericInput("q","Q (L/min)",value = 26,width=120),
                 numericInput("sato2a","SatO2_a (%)",value = 96,width=120),
                 numericInput("sato2v","SatO2_v (%)",value = 23,width=120),
                 numericInput("paco2","PaCO2 (mmHg)",value = 35,width=120),
                 numericInput("pvco2","PvCO2 (mmHg)",value = 50,width=120),
                 numericInput("pha","pH arterial",value = 7.34,width=120),
                 numericInput("phv","pH venous",value = 7.21,width=120),
                 width=3),
             mainPanel(
                 fluidRow(column(12,tableOutput("inDataExcel"))),
                 fluidRow(column(12,tableOutput("inData"))),
                 fluidRow(column(12,tableOutput("alldata"))),
                 hr(),
                 fluidRow(column(6,plotOutput("plotqvo2")),
                          column(6,plotOutput("plotvavo2"))),
                 fluidRow(column(6,plotOutput("plotdlvo2")),
                          column(6,plotOutput("plotdmvo2")))
               )
             )
    ),
    tabPanel("Algorithm 2",
             sidebarLayout(
               sidebarPanel(width = 3,
                            actionButton("update", "Update Plots"),
               ),
               
               mainPanel(
                 # plotOutput("plotqvo22"),
               )
             )),
    tabPanel("Patient Simulation",
             sidebarLayout(
               sidebarPanel(width = 3,
                            actionButton("update", "Update Plots"),
               ),
               
               mainPanel(
                 # plotOutput("plotqvo22"),
               )
             ))
  )
)
server <- function(input, output,session) {
  output$inData <- renderTable( indata())
  indata <- eventReactive(input$newpatient, {
    if(input$newpatient>0){
      newrow <- isolate(c(input$newpatient, input$vo2,input$vco2,input$pao2,input$pvo2,
                          input$hb,input$q,input$sato2a,input$sato2v,input$paco2,input$pvco2,input$pha,input$phv,
                          input$newgroup))
      newtab <- as.data.frame(matrix(data=as.numeric(newrow),ncol=14,byrow=T))
      newtab[14] <- input$newgroup
      colnames(newtab)<-tolower(c("id","VO2","VCO2","PaO2","PvO2","Hb","Q","satao2","satcvo2","PaCO2","PvCO2","pha","phv","group"))
      #check if any of the optional variables is there and remove it otherwise
      if (is.na(newtab$q)) {newtab <- newtab[,names(newtab) != 'q']}
      if (is.na(newtab$vo2)) {newtab <- newtab[,names(newtab) != 'vo2']}
      if (is.na(newtab$vco2)) {newtab <- newtab[,names(newtab) != 'vco2']}
      newtab$id <- paste(as.integer(input$newpatient),input$newgroup,sep="_")
      #append new calculations to old data
      my_data <<- plyr::rbind.fill(my_data,calc_all(newtab))
      #show new patient
      newtab
    }
  }, ignoreNULL = FALSE)
  
  output$inDataExcel <- renderTable({
    inFile <- input$ul
    if (is.null(inFile)){
      calcdata <<-NULL
      return(NULL)
    }
    inDataExcel <- read_excel(inFile$datapath)
    colnames(inDataExcel) <- tolower(colnames(inDataExcel))
    for(i in 1:nrow(inDataExcel)) {
      row <- inDataExcel[i,]
      # do stuff with row
      calcdata <<- rbind(calcdata,calc_all(row))
    }
    calcdata
  })
  
  #print list of all new patients
  output$alldata <- renderTable( df())
  df <- eventReactive(input$newpatient, {
    my_data[seq(38,dim(my_data)[1]),]
  })
  
  #reset my_data
  observeEvent(input$reset, {
    my_data <<- precomp_data
    session$reload()
  })
  
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
  
}
shinyApp(ui = ui, server = server)
