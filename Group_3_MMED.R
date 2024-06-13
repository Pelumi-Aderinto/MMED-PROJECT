rm(list=ls(all=T))
library(deSolve) ## differential equation solver library

## CONTROL PARAMETERS
SAVEPLOTS <- FALSE
if(SAVEPLOTS){
  try(system('mkdir Figures'))
  if(file.exists("Figures")){
    plotPath <- 'Figures'
  }else{
    warning("Figures directory does not exist; plotPath will not be defined, and figures will be plotted but not saved.")
  }
}
seiarModel <- function(t,y,params){
  with(c(as.list(y),params), {
    N <- S + E + Isurv + Idead + A + R ## Total population size
    nu<-mu*N                           ## Births (zero in current model, as are deaths)
    dS <- nu - beta*S*(Isurv+Idead)/N - mu*S ## Susceptible
    dE <- beta*S*(Isurv+Idead)/N - mu*E - sigma*E ## Exposed (incubating)
    dIsurv <- symp*(1-cfr)*sigma*E - mu*Isurv - gamma*Isurv ## symptomatic surviving infections
    dIdead <- symp*cfr*sigma*E - mu*Idead - gamma*Idead  ## symptomatic dying infections
    dA <- (1-symp)*sigma*E - mu*A - gamma*A ## asymptomatic individuals
    dR <- gamma*(Isurv + A) - mu*R          ## recovered & immune
    dcumExp <- beta*S*(Isurv+Idead)/N  ## cumulative exposed
    dcumInc <- symp*sigma*E  ## cumulative incidence of symptomatic infections
    dcumMort <- gamma*Idead ## cumulative mortality
    return(list(c(dS=dS,dE=dE,dIsurv=dIsurv,dIdead=dIdead,dA=dA,dR=dR, dcumExp=dcumExp, dcumInc=dcumInc, dcumMort=dcumMort)))
  })
}


R0 <- 2
N0 <- 4*10^6 ## Liberia's population
N0K <- N0/10^3 ## population size in thousands
## Initialize with one symptomatic infectious individual who will die.
init <- c(S = N0-1,
          E = 0,
          Isurv = 0,
          Idead = 1,
          A = 0,
          R = 0,
          cumExp = 0, cumInc = 0, cumMort = 0)
times<-seq(0,500,1) ## Simulate for 500 days

param.vals<-c( ## Other parameters
  beta= NA, ## Calculated based on R0 and other parameters, see below.
  N0=N0,
  mu= 0, ## Assume no birth/death for now, though it doesn't affect this toy model much. For 50 yr life expect mu=.02/365.25
  sigma=1/9.1, ## progression rate = 1/incubation or 1/latent period (assumed to be the
  ## same for Ebola). Lancet estimat 9.1 days; CDC estimate 6 days.
  symp = NA, ## symptomatic proportion
  gamma=1/6,  ## 1/infectious period. CDC estimate 6 days
  cfr = .7) ## case fatality rate. Lancet for Liberia = 72.3%

Ro.calc<-function(params) { ## Analytical solution for R0
  with(as.list(params),
       beta*(sigma/(mu+sigma)) * symp * (cfr/(mu+gamma) + (1-cfr)/(mu+gamma))
  )}
beta.calc<-function(Ro,params) { ## Solve above function for beta
  with(as.list(params),
       Ro/((sigma/(mu+sigma)) * symp * (cfr/(mu+gamma) + (1-cfr)/(mu+gamma)))
  )}

runSEIAR <- function(sympProp, paramVals = param.vals, basicReproNum = R0, browse=F){
  if(browse) browser()
  paramVals['symp'] <- sympProp
  paramVals['beta'] <- beta.calc(basicReproNum,paramVals)
  print(paste("Calculated beta value for ",sympProp*100,"% symptomatic: ",round(paramVals['beta'],3),".",sep=""))
  tc <- data.frame(lsoda(init, times, seiarModel, paramVals))  ## Run ODE model
  tc$N <-  rowSums(tc[,c('S','E','Isurv','Idead','A','R')])    ## Calculate total population size
  tc[,-1] <- tc[,-1] / 10^3                                    ## Show numbers (other than time) in thousands
  tc$Reff <- R0*(tc$S/tc$N)                                    ## Calculate R_effective
  return(tc)
}

sympVals <- c(1,.8,0.5,.3)
tcSymp <- runSEIAR(sympVals[1])     ## 100% symptomatic
tcSymp_cumInc <- tail(tcSymp$cumInc, n=1)
tcAsymp <- runSEIAR(sympVals[2])   ## 80% symptomatic
tcAsymp_cumInc <- tail(tcAsymp$cumInc, n=1)
tcAsymp1 <- runSEIAR(sympVals[3])  ## 50% symptomatic
tcAsymp1_cumInc <- tail(tcAsymp1$cumInc, n=1)
tcAsymp2 <- runSEIAR(sympVals[4])  ## 30% symptomatic
tcAsymp2_cumInc <- tail(tcAsymp2$cumInc, n=1)
## Compare calculated beta values. Note beta is bigger to make up for lower symptomatic proportion.

tail(tcSymp$N + tcSymp$cumMort) ## Check (population size + cumulative mortality) is constant
tail(tcAsymp$N + tcAsymp$cumMort) ## Check (population size + cumulative mortality) is constant
tail(tcAsymp1$N + tcAsymp1$cumMort)
tail(tcAsymp2$N + tcAsymp2$cumMort)
tshow <- c(1:5,nrow(tcSymp))
show <- c('cumInc','cumExp','cumMort')
tcAsymp[tshow,show]/tcSymp[tshow,show]/tcAsymp1[tshow,show]/tcAsymp2[tshow,show]


## Set calendar time. Meltzer et al. (2014) MMWR estimates 3915 EVD cases in Liberia, Aug 28, 2014 after
## applying a correction factor for unreported EVD cases. We set the day in our model closest to this value
## to be August 28.
aug28 <- tcSymp$time[which.min(abs(tcSymp$cumInc*10^3 - 3915))]
days <- as.Date('2014-08-28') + (tcSymp$time - aug28)
tcSymp$days <- tcAsymp$days <- tcAsymp1$days <- tcAsymp2$days <- days ## add calendar days to both modeled time series



crit <- function(vaccPropNeeded = 0.5, ## required to reduce Reff < 1 for R0 = 2
                 cumIncSymp, ## cumulative incidence of symptomatic cases when vaccination starts
                 totalPop = N0K,  ## total population
                 cfr = 0.6,  ## case fatality rate
                 symp = 0.5) { ## symptomatic proportion
  # Calculate the number of asymptomatic immune individuals
  numAsympImmune <- cumIncSymp * (1 - symp) / symp
  # Calculate the number of symptomatic survivors (immune individuals)
  numSympSurvImmune <- (1 - cfr) * cumIncSymp
  # Calculate the total immune individuals
  totalImmune <- numAsympImmune + numSympSurvImmune
  # Calculate the proportion immune
  propImmune <- totalImmune / totalPop
  # Calculate the required vaccination coverage
  return(pmax(vaccPropNeeded - propImmune, 0) / (1 - propImmune))
}

# Calculate the vaccination threshold using dynamic cumIncSymp
vacc_threshold_Symp_1 <- crit(cumIncSymp = tcSymp_cumInc, totalPop = N0K, symp = 1)
vacc_threshold_Symp_2 <- crit(cumIncSymp = tcAsymp_cumInc, totalPop = N0K, symp = 0.8)
vacc_threshold_Symp_3 <- crit(cumIncSymp = tcAsymp1_cumInc, totalPop = N0K, symp = 0.5)
vacc_threshold_Symp_4 <- crit(cumIncSymp = tcAsymp2_cumInc, totalPop = N0K, symp = 0.3)

# Display the vaccination thresholds
print(vacc_threshold_Symp_1)
print(vacc_threshold_Symp_2)
print(vacc_threshold_Symp_3)
print(vacc_threshold_Symp_4)

####################################################################################################
## Figures
####################################################################################################

## Figure in Lancet letter
sel <- days > as.Date('2014-09-01') & days < as.Date('2015-01-10')  ## show Sep 2014 - Feb 2015
if(SAVEPLOTS) png(file.path(plotPath, 'rel cumInc 2 panel.png'), w = 4, 5, units='in', res = 300)
par('ps' = 11, mar = c(4.5, 5, .5, 1), lwd = 2, mgp = c(3, 1, 0), mfrow = c(2, 1))
## Comparing cumulative EVD cases with and without accounting for asymptomatic proportion.
mains <- c("(A) Cumulative # of Cases", '(B) Vaccination Coverage Needed for Elimination')
mains <- rep('', 2)
## Percent difference in projected cumulative incidence of symptomatic EVD
## cases between symptomatic an asymptomatic models (Panel A)
ylab <- 'Cumulative Cases\n(Thousands)'
plot(days[sel], tcSymp$cumInc[sel], type = 'n', xlab = '', ylab = ylab, las = 1, main = '', xaxt = 'n')
mtext(mains[1], side = 3, line = 1)
lines(days[sel], tcSymp$cumInc[sel], lty = 1, col = "red")
lines(days[sel], tcAsymp$cumInc[sel], lty = 2, col = "orange")
lines(days[sel], tcAsymp1$cumInc[sel], lty = 3, col = "blue")
lines(days[sel], tcAsymp2$cumInc[sel], lty = 4, col = "green")

mth <- seq.Date(as.Date('2014-01-01'), as.Date('2015-02-01'), by = 'month')
axis.Date(side = 1, at = mth, format = '%b %e', las = 2)
legend('topleft', legend = paste0(sympVals * 100, '% Symptomatic'), lty = 1:4, col = c("red", "orange", "blue", "green"), bty = 'n', cex = 1, bg = 'white')

## Total vaccination coverage needed (Panel B)
cRange <- seq(0, 1, by = 0.2)
plot(cRange * N0K, crit(cumIncSymp = cRange*N0K, symp = 1), type = "l", lty = 1, col = "red", ylim = c(0, .5),
     xlab = "Cumulative Cases\n(Thousands)", ylab = "Target Vaccination \nCoverage")
lines(cRange * N0K, crit(cumIncSymp = cRange*N0K, symp = 0.8), lty = 2, col = "orange")
lines(cRange * N0K, crit(cumIncSymp = cRange*N0K, symp = 0.5), lty = 3, col = "blue")
lines(cRange * N0K, crit(cumIncSymp = cRange*N0K, symp = 0.3), lty = 4, col = "green")
mtext(mains[2], side = 3, line = 1, at = 0.09 * N0K)
legend('topright', legend = paste0(sympVals * 100, '% Symptomatic'), lty = 1:4, col = c("red", "orange", "blue", "green"), bty = 'n', cex = 1, bg = 'white')
if(SAVEPLOTS) graphics.off()