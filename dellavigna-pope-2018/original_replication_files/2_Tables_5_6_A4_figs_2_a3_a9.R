###############################
# This script generates the results for most of Table 5 
# Columns A1,A3,B1,B1,B2,B4,B5,B7
# Also creates the forecasters individual implied parameters (for Appendix figure 9)
# Table 6: Implied forecast probabilities, and Panel B
# Appendix Table 4: B2, B4, B6, B8
#
# Written by: Avner Shlain
###############################
# set working directory here
setwd("C:/Users/Avner/Dropbox/forecastExp/Data Analysis/Analysis/behavioral/")
rm(list=ls())

library(rootSolve)
library(dplyr)
library(stargazer)
library(foreign)
library(readstata13)
library(stringr)
library(grid)
#######
# Bootstrap Sample the results
#######
set.seed(3212017)

# Read raw data, sort by treatment
raw.data <- read.dta13(file = "dtafiles/MTurkCleanedData.dta") %>% arrange(treatment) %>%
  select(buttonpresses, treatment, treatmentname)

# Keep track of treatments
bs_legend <- data.frame(unique(raw.data[,2:3])) %>% arrange(treatment)

#Write files for later robustness calculations of the low pay treatment (read by 3_Table5_panelA_Col4_LowPaySimulations.m)
lowpay <- filter(raw.data,treatmentname == "Very Low Pay") %>% select(buttonpresses)
write.table(as.matrix(lowpay),file = "dtafiles/buttonpressesLowPay.csv",row.names = FALSE,col.names = FALSE)

nopay <- filter(raw.data,treatmentname == "No Payment") %>% select(buttonpresses)
write.table(as.matrix(nopay),file = "dtafiles/buttonpressesNoPay.csv",row.names = FALSE,col.names = FALSE)

# 1000 bootstrap samples, initalize:
bs_sample_full <- matrix(data = 0,nrow=1000,ncol=18)
colnames(bs_sample_full) <- bs_legend$treatmentname

# Sampling - 550 workers with replacements
for (i in 1:nrow(bs_sample_full)) {
  temp_means <- raw.data %>%
    group_by(treatment) %>%
    sample_n(size = 550,replace = TRUE) %>%
    summarise(means = mean(buttonpresses))
  bs_sample_full[i,] <- temp_means$means
}

write.csv(bs_sample_full,file = "dtafiles/bootstrap1000A.csv", row.names = FALSE)

sample_means <- raw.data %>%
  group_by(treatment) %>%
  summarise(means = round(mean(buttonpresses)))
sample_means$treatment <- c("p1","p10","p0","p4","ge","p01","ch01","ch10","t2",
                      "t4","g40","l40","g80","pr1","pr50","ci","fb","try")
mean_effort <- data.frame(t(sample_means$means))
names(mean_effort) <- sample_means$treatment


bs_sample <- data.frame(bs_sample_full)
names(bs_sample) <- c("p1","p10","p0","p4","ge","p01","ch01","ch10","t2",
                      "t4","g40","l40","g80","pr1","pr50","ci","fb","try")


#########
## Read forecasts
#########

forecasts <- read.dta(file = "input/ExpertForecast.dta")[,c(1,3:21)] %>%
  filter(completed == "Yes") %>% select(-completed)
names(forecasts) <- c("id","p0","p1","p10","p4","p01","ch01","ch10","t2","t4",
                      "g40","l40","g80","pr1","pr50","ci","fb","try","ge")
forecasts <- forecasts[,names(bs_sample)]
                      
#################
## FUNCTIONS
#################

## 3 EQUATIONS WITH 3 UNKOWNS:
focs.eqs <- function(theta,params,spec,...) {
  log.s <- theta[1]
  log.k <- theta[2]
  log.g <- theta[3]
  e1 <- params[1,1]
  p1 <- params[1,2]
  e2 <- params[2,1]
  p2 <- params[2,2]
  e3 <- params[3,1]
  p3 <- params[3,2]

  if (spec == "exp") {
  equations <- c(Eq.p0 = exp(log.g)*(e1) - log(exp(log.s)+p1) + log.k,
                 Eq.p1 = exp(log.g)*(e2) - log(exp(log.s)+p2) + log.k,
                 Eq.p10 = exp(log.g)*(e3) - log(exp(log.s)+p3) + log.k)
  }
  if (spec == "power") {
  equations <- c(Eq.p0 = exp(log.g)*log(e1) - log(exp(log.s)+p1) + log.k,
                 Eq.p1 = exp(log.g)*log(e2) - log(exp(log.s)+p2) + log.k,
                 Eq.p10 = exp(log.g)*log(e3) - log(exp(log.s)+p3) + log.k)
  }
  return(equations)
}

# SOLVE THE EQUATIONS SET NUMERICALLY
Min.Dist <- function(E=c(1521,2029,2175),P=c(0,0.01,0.10),Spec,...) {
  Parameters <- matrix(c(E,P),ncol=2)
  
  #Guesstimate with s=0
  if (Spec == "power") {
    log_k <- (log(P[3]) - log(P[2])*log(E[3])/log(E[2]))/(1 - log(E[3])/log(E[2]))
    log_gamma <- log((log(P[2]) - log_k)/log(E[2]))
    log_s <- exp(log_gamma)*log(E[1]) + log_k
  }

  if (Spec == "exp") {
    log_k <- (log(P[3]) - log(P[2])*(E[3])/(E[2]))/(1 - (E[3])/(E[2]))
    log_gamma <- log((log(P[2]) - log_k)/(E[2]))
    log_s <- exp(log_gamma)*(E[1]) + log_k
  }
  
  #Solve with the above as initial guess
  rootz <- multiroot(focs.eqs,
                     start = c(log_s,log_k,log_gamma),parms=Parameters,spec=Spec,
                     rtol = 1e-16,atol=1e-16,ctol=1e-16)
  
  roots <- rootz$root
  roots <- exp(roots)
  return(c(s = roots[1], k = roots[2], gamma = roots[3],precision = rootz$estim.precis))
}

# OUT OF SAMPLE PREDICTIONS
oos <- function(s,k,g, spec) {
  e <- c(NA,NA)
  # predict 4 cents per 100 and 1 cent per 1000
  
  if (spec == "power") { 
    e[1] <- ((s + 0.04)/k)^(1/g)
    e[2] <- ((s + 0.001)/k)^(1/g)
  }
  if (spec == "exp") { 
    e[1] <- log((s + 0.04)/k)*(1/g)
    e[2] <- log((s + 0.001)/k)*(1/g)
  }
  return(c(e4 = e[1], e01 = e[2]))
}


# PARAMETERS CALCULATION
Beh.Params <- function(E,s,k,g,spec,solve.again=FALSE) {
  #Should we solve for the behavioral params again
  if (solve.again == TRUE) {
    temp <- Min.Dist(as.numeric(E[c(3,1,2)]),Spec=spec)
    s <- temp[1]
    k <- temp[2]
    g <- temp[3]
  }
  if (spec == "exp") {
    Eg <- exp(E*g)
    # alpha - altruism; a - warm glow
    ech01g <- Eg[7] #E$ch01
    ech10g <- Eg[8] #E$ch10
    Alpha <- 100/9*k*(ech10g - ech01g)
    A <- 100*k*ech01g - 100*s - Alpha
    
    # DeltaS - gift exchange shift
    S.ge <- k*Eg[5] - s
    
    # time params - delta and beta
    Delta <- sqrt((k*Eg[10] - s)/(k*Eg[9] - s))
    Beta <- 100*(k*Eg[9] - s)/(Delta^2)
    
    # Probability weighting
    Pi1 <- (k*Eg[14] - s)/(0.01)
    
  }
  if (spec == "power") {

  # alpha - altruism; a - warm glow
  ech01 <- E[7] #E$ch01
  ech10 <- E[8] #E$ch10
  Alpha <- 100/9*k*(ech10^g - ech01^g)
  A <- 100*k*ech01^g - 100*s - Alpha

  # DeltaS - gift exchange shift
  S.ge <- k*E[5]^g - s

  # time params - delta and beta
  Delta <- sqrt((k*E[10]^g - s)/(k*E[9]^g - s))
  Beta <- 100*(k*E[9]^g - s)/(Delta^2)

  # Probability weighting
  Pi1 <- (k*E[14]^g - s)/(0.01)
  }
  
  # Risk aversion - lambda
  if (E[13]-E[11] > 10) {
    Lambda <- 2*(E[12] - E[11])/(E[13] - E[11]) + 1
    if (Lambda < 0) {
      Lambda <- NaN
    }
  } else {
    Lambda <- NaN
  }
  return(c(alpha=Alpha, a=A, s.ge = S.ge, beta = Beta, delta = Delta, 
           pi1 = Pi1, lambda = Lambda))
}

#########################################################
####      TABLE 5         ###############################
#########################################################

# Panel A
T5A <- data.frame(rows = c("gamma","gamma se","k","k se","s","s se","Implied 4","Implied low-pay"),
                  col1= NA, col2 = NA, col3 = NA, col4 = NA)

spec = "power"
results <- apply(X = bs_sample[,c("p0","p1","p10")],MARGIN = 1,Min.Dist, Spec="power")
results <- data.frame(t(results))
names(results) <- c("s","k","gamma","precision")

results_clean <- results %>%
  filter(precision < 1) %>%
  filter(k > 0) %>%
  filter(s > 0)


T5A$col1[c(5,3,1)] <- Min.Dist(Spec="power")[1:3]
T5A$col1[c(6,4,2)] <- apply(results_clean[,1:3],2,sd)
T5A$col1[7:8] <- oos(T5A$col1[5],T5A$col1[3],T5A$col1[1],"power")

spec = "exp"
results <- apply(X = bs_sample[,c("p0","p1","p10")],MARGIN = 1,Min.Dist, Spec="exp")
results <- data.frame(t(results))
names(results) <- c("s","k","gamma","precision")

results_clean <- results %>%
  filter(precision < 1) %>%
  filter(k > 0) %>%
  filter(s > 0)



T5A$col3[c(5,3,1)] <- Min.Dist(Spec="exp")[1:3]
T5A$col3[c(6,4,2)] <- apply(results_clean[,1:3],2,sd)
T5A$col3[7:8] <- oos(T5A$col3[5],T5A$col3[3],T5A$col3[1],"exp")


# Read NLS results - column 2
nls_power <- read.csv("output/NLS_results_Table5_POWER.csv",as.is=TRUE)
T5A$col2[c(1,3,5)] <- as.numeric(str_sub(nls_power[c(3,6,9),2],2)) #read gamma,k,s
T5A$col2[3] <- T5A$col2[3]/as.numeric(str_sub(nls_power[32,2],2))  #rescale k
T5A$col2[5] <- T5A$col2[5]/as.numeric(str_sub(nls_power[33,2],2))  #rescale s
T5A$col2[c(2,4,6)] <- str_sub(nls_power[c(4,7,10),2],2)            #read se.s NOT RESCALED
T5A$col2[7:8] <-  (as.numeric(str_sub(nls_power[c(42,43),2],2)))

# Read NLS results - column 4
nls_expon <- read.csv("output/NLS_results_Table5_EXPON.csv",as.is=TRUE)
T5A$col4[c(1,3,5)] <- as.numeric(str_sub(nls_expon[c(3,6,9),2],2)) #read gamma,k,s
T5A$col4[3] <- T5A$col4[3]/as.numeric(str_sub(nls_expon[32,2],2))  #rescale k
T5A$col4[5] <- T5A$col4[5]/as.numeric(str_sub(nls_expon[33,2],2))  #rescale s
T5A$col4[c(2,4,6)] <- str_sub(nls_expon[c(4,7,10),2],2)            #read se.s NOT RESCALED
T5A$col4[7:8] <-  (as.numeric(str_sub(nls_expon[c(42,43),2],2)))


write.table(x = "Table 5 Panel A",file = "output/Table5.csv",sep = ",",row.names = FALSE, col.names=FALSE, append=FALSE)
write.table(x = T5A,file = "output/Table5.csv",sep = ",",row.names = FALSE , col.names=FALSE, append=TRUE)
write.table(x = c("","Table 5 Panel B"),file = "output/Table5.csv",sep = ",",row.names = FALSE , col.names=FALSE, append=TRUE)
T5A


#T5 Panel B
T5B <- data.frame(rows = c("alpha","alpha ci","a","a ci","gift exchange","gift exchange ci","beta","beta ci","delta","delta ci"),
                  col1=NA, col2 = NA, col3=NA, col4 = NA, col5 = NA, col6 = NA, col7 = NA)
# Point estimates - column 1
T5B$col1[c(1,3,5,7,9)] <- Beh.Params(E=mean_effort, s=T5A$col1[5], 
                                     k=T5A$col1[3], g=T5A$col1[1],  
                                     spec="power")[1:5]

# confidence intervals - column 1, based on bootstrap
results <- apply(X = bs_sample,MARGIN = 1,Beh.Params,s=T5A$col1[5], 
                 k=T5A$col1[3], g=T5A$col1[1], spec="power", solve.again = TRUE)
results <- data.frame(t(results))[1:5]
names(results) <- c("alpha","a","ge","beta","delta")
CIs.p <- signif(apply(results,2,quantile,p=c(0.025,0.975),na.rm=TRUE),2)
T5B$col1[c(2,4,6,8,10)] <- paste("(",apply(CIs.p,2,paste,collapse=","),")",sep="")

# column 2 - forecasts 25p, median, 75p
results <- apply(X = forecasts,MARGIN = 1,Beh.Params,s=T5A$col1[5], 
                 k=T5A$col1[3], g=T5A$col1[1], spec="power")
results <- data.frame(t(results))[1:5]
results <- results[complete.cases(results),]
names(results) <- c("alpha","a","ge","beta","delta")
T5B$col2[c(1,3,5,7,9)] <- signif(apply(results,2,quantile,p=0.5),2)
IQr <- signif(apply(results,2,quantile,p=c(0.25,0.75)),2)
T5B$col2[c(2,4,6,8,10)] <- paste("(",apply(IQr,2,paste,collapse=","),")",sep="")

# column 3 - nls estimates
T5B$col3[c(1,3,5,7,9)] <- as.numeric(str_sub(nls_power[c(21,24,12,15,18),3],2)) #read alpha,a,s.ge,beta,delta
ses <- as.numeric(str_sub(nls_power[c(22,25,13,16,19),3],3,-2))
T5B$col3[c(2,4,6,8,10)] <- paste(rep("(",5),signif(T5B$col3[c(1,3,5,7,9)] - 1.96*ses,3),
      ",",signif(T5B$col3[c(1,3,5,7,9)] + 1.96*ses,3),rep(")",5),sep="")

# Point estimates - column 4, based on Min-Dist estimates
T5B$col4[c(1,3,5,7,9)] <- Beh.Params(E=mean_effort, s=T5A$col3[5], 
                                     k=T5A$col3[3], g=T5A$col3[1],  
                                     spec="exp")[1:5]

# confidence intervals - column 4, based on bootstrap with Min-Dist estimates
results <- apply(X = bs_sample,MARGIN = 1,Beh.Params,s=T5A$col3[5], 
                 k=T5A$col3[3], g=T5A$col3[1], spec="exp", solve.again = TRUE)
results <- data.frame(t(results))[1:5]
names(results) <- c("alpha","a","ge","beta","delta")
CIs.p <- signif(apply(results,2,quantile,p=c(0.025,0.975)),2)
T5B$col4[c(2,4,6,8,10)] <- paste("(",apply(CIs.p,2,paste,collapse=","),")",sep="")

# column 5 - forecasts 25p, median, 75p, based on Min-Dist estimates, EXP
results <- apply(X = forecasts,MARGIN = 1,Beh.Params,s=T5A$col3[5], 
                 k=T5A$col3[3], g=T5A$col3[1], spec="exp")
results <- data.frame(t(results))[1:5]
results <- results[complete.cases(results),]
names(results) <- c("alpha","a","ge","beta","delta")
T5B$col5[c(1,3,5,7,9)] <- signif(apply(results,2,quantile,p=0.5),2)
IQr <- signif(apply(results,2,quantile,p=c(0.25,0.75)),2)
T5B$col5[c(2,4,6,8,10)] <- paste("(",apply(IQr,2,paste,collapse=","),")",sep="")

# column 6 - nls estimates
T5B$col6[c(1,3,5,7,9)] <- as.numeric(str_sub(nls_expon[c(21,24,12,15,18),3],2)) #read alpha,a,s.ge,beta,delta
ses <- as.numeric(str_sub(nls_expon[c(22,25,13,16,19),3],3,-2))
T5B$col6[c(2,4,6,8,10)] <- paste(rep("(",5),signif(T5B$col6[c(1,3,5,7,9)] - 1.96*ses,3),
                                 ",",signif(T5B$col6[c(1,3,5,7,9)] + 1.96*ses,3),rep(")",5),sep="")

# column 7 - forecasts 25p, median, 75p, based on NLS estimates, EXP
NLS.s <- as.numeric(str_sub(nls_expon[9,3],2))/as.numeric(str_sub(nls_expon[33,2],2))
NLS.k <- as.numeric(str_sub(nls_expon[6,3],2))/as.numeric(str_sub(nls_expon[32,2],2))
NLS.g <- as.numeric(str_sub(nls_expon[3,3],2))

results <- apply(X = forecasts,MARGIN = 1,Beh.Params,s=NLS.s, 
                 k=NLS.k, g=NLS.g, spec="exp")
results <- data.frame(t(results))[1:5]
results <- results[complete.cases(results),]
names(results) <- c("alpha","a","ge","beta","delta")
T5B$col7[c(1,3,5,7,9)] <- signif(apply(results,2,quantile,p=0.5),2)
IQr <- signif(apply(results,2,quantile,p=c(0.25,0.75)),2)
T5B$col7[c(2,4,6,8,10)] <- paste("(",apply(IQr,2,paste,collapse=","),")",sep="")

T5B <- apply(T5B,2,as.character)
write.table(x = T5B,file = "output/Table5.csv",sep = ",", col.names=FALSE, append=TRUE,row.names = FALSE)

#########################################################
####      TABLE 6         ###############################
#########################################################
T6A <- data.frame(rows = c("25th","median","75th"),
                  col4= NA, col5 = NA, col6 = NA)

# Column 4 - take NLS estimates, get implied weights
nls_expon <- read.csv("output/estimation_results_T6_EXPON.csv",as.is=TRUE)
NLS.s <- as.numeric(str_sub(nls_expon[9,2],2))/as.numeric(str_sub(nls_expon[24,2],2))
NLS.k <- as.numeric(str_sub(nls_expon[6,2],2))/as.numeric(str_sub(nls_expon[23,2],2))
NLS.g <- as.numeric(str_sub(nls_expon[3,2],2))

results <- apply(X = forecasts,MARGIN = 1,Beh.Params,s=NLS.s, 
                 k=NLS.k, g=NLS.g, spec="exp")
results <- data.frame(t(results))[6]
curv1 <- results
results <- results[complete.cases(results),]
names(results) <- c("Pi(1 percent)")
T6A$col4 <- signif(quantile(results,p=c(0.25,0.5,0.75)),2)#,"%",sep="")

# Column 5 - take NLS estimates, get implied weights
NLS.s <- as.numeric(str_sub(nls_expon[9,3],2))/as.numeric(str_sub(nls_expon[24,3],2))
NLS.k <- as.numeric(str_sub(nls_expon[6,3],2))/as.numeric(str_sub(nls_expon[23,3],2))
NLS.g <- as.numeric(str_sub(nls_expon[3,3],2))

results <- apply(X = forecasts,MARGIN = 1,Beh.Params,s=NLS.s, 
                 k=NLS.k, g=NLS.g, spec="exp")
results <- data.frame(t(results))[6]
curv88 <- results
results <- results[complete.cases(results),]
names(results) <- c("Pi(1 percent)")
T6A$col5 <- signif(quantile(results,p=c(0.25,0.5,0.75)),2)#,"%",sep="")


# Column 6 - take NLS estimates, get implied weights
NLS.s <- as.numeric(str_sub(nls_expon[9,4],2))/as.numeric(str_sub(nls_expon[24,4],2))
NLS.k <- as.numeric(str_sub(nls_expon[6,4],2))/as.numeric(str_sub(nls_expon[23,4],2))
NLS.g <- as.numeric(str_sub(nls_expon[3,4],2))

results <- apply(X = forecasts,MARGIN = 1,Beh.Params,s=NLS.s, 
                 k=NLS.k, g=NLS.g, spec="exp")
results <- data.frame(t(results))[6]
curv44 <- results
results <- results[complete.cases(results),]
names(results) <- c("Pi(1 percent)")
T6A$col6 <- signif(quantile(results,p=c(0.25,0.5,0.75)),2)#,"%",sep="")


write.table(x = "Table 6 Panel A - Pi results",file = "output/Table6.csv",sep = ",",row.names = FALSE, col.names=FALSE, append=FALSE)
write.table(x = T6A,file = "output/Table6.csv",sep = ",",row.names = FALSE , col.names=FALSE, append=TRUE)
write.table(x = c("","Table 6 Panel B - Lambda estimation and predictions"),
  file = "output/Table6.csv",sep = ",",row.names = FALSE , col.names=FALSE, append=TRUE)
# PANEL B
T6B <- data.frame(rows = c("lambda","ci"),col1=NA,col2=NA)

# T6B- column 1
T6B$col1[1] <- signif(Beh.Params(E = sample_means$means, s=0.1, k=1, g=1, spec="power")[7],3)
results <- apply(X = bs_sample,MARGIN = 1,Beh.Params,s=NLS.s, 
      k=NLS.k, g=NLS.g, spec="power",solve.again=TRUE)
results <- data.frame(t(results))[7] #lambda 
results <- results[complete.cases(results),]
names(results) <- c("Lambda")
IQr <- signif(quantile(results,p=c(0.025,0.975)),2)
T6B$col1[2] <- paste("(",paste(IQr,collapse=","),")",sep="")

# T6B- column 2
results <- apply(X = forecasts,MARGIN = 1,Beh.Params,s=NLS.s, 
                 k=NLS.k, g=NLS.g, spec="power")
results <- data.frame(t(results))[7] #lambda 
results <- results[complete.cases(results),]
names(results) <- c("Lambda")
T6B$col2[1] <- signif(quantile(results,p=0.5),3)
IQr <- signif(quantile(results,p=c(0.25,0.75)),2)
T6B$col2[2] <- paste("(",paste(IQr,collapse=","),")",sep="")

write.table(x = T6B,file = "output/Table6.csv",sep = ",",row.names = FALSE , col.names=FALSE, append=TRUE)


#########################################################
####    APPENDIX TABLE 4  PANEL B       #################
#########################################################
TA4B <- data.frame(rows = c("alpha","alpha ci","a","a ci","gift exchange","gift exchange ci","beta","beta ci","delta","delta ci"),
                  col1=NA, col2 = NA, col3=NA, col4 = NA, col5 = NA, col6 = NA, col7 = NA, col8 = NA)

nls_expon <- read.csv("output/estimation_results_AT4_.csv",as.is=TRUE)

# TA4B - column 1 - NLS estimates
TA4B$col1[c(1,3,5,7,9)] <- as.numeric(str_sub(nls_expon[c(18,21,9,12,15),2],2)) #read alpha,a,s.ge,beta,delta
ses <- as.numeric(str_sub(nls_expon[c(19,22,10,13,16),2],3,-2))
TA4B$col1[c(2,4,6,8,10)] <- paste(rep("(",5),signif(TA4B$col1[c(1,3,5,7,9)] - 1.96*ses,3),
                                 ",",signif(TA4B$col1[c(1,3,5,7,9)] + 1.96*ses,3),rep(")",5),sep="")


# TA4B - column 2
NLS.s <- as.numeric(str_sub(nls_expon[6,2],2))/as.numeric(str_sub(nls_expon[33,2],2))
NLS.k <- as.numeric(str_sub(nls_expon[3,2],2))/as.numeric(str_sub(nls_expon[32,2],2))
NLS.g <- 0.01

results <- apply(X = forecasts,MARGIN = 1,Beh.Params,s=NLS.s, 
                 k=NLS.k, g=NLS.g, spec="exp")
results <- data.frame(t(results))[1:5]
results <- results[complete.cases(results),]
names(results) <- c("alpha","a","ge","beta","delta")
TA4B$col2[c(1,3,5,7,9)] <- signif(apply(results,2,quantile,p=0.5),2)
IQr <- signif(apply(results,2,quantile,p=c(0.25,0.75)),2)
TA4B$col2[c(2,4,6,8,10)] <- paste("(",apply(IQr,2,paste,collapse=","),")",sep="")

# TA4B - column 3 - NLS estimates
TA4B$col3[c(1,3,5,7,9)] <- as.numeric(str_sub(nls_expon[c(18,21,9,12,15),3],2)) #read alpha,a,s.ge,beta,delta
ses <- as.numeric(str_sub(nls_expon[c(19,22,10,13,16),3],3,-2))
TA4B$col3[c(2,4,6,8,10)] <- paste(rep("(",5),signif(TA4B$col3[c(1,3,5,7,9)] - 1.96*ses,3),
                                  ",",signif(TA4B$col3[c(1,3,5,7,9)] + 1.96*ses,3),rep(")",5),sep="")

# TA4B - column 4
NLS.s <- as.numeric(str_sub(nls_expon[6,3],2))/as.numeric(str_sub(nls_expon[33,3],2))
NLS.k <- as.numeric(str_sub(nls_expon[3,3],2))/as.numeric(str_sub(nls_expon[32,3],2))
NLS.g <- 0.02

results <- apply(X = forecasts,MARGIN = 1,Beh.Params,s=NLS.s, 
                 k=NLS.k, g=NLS.g, spec="exp")
results <- data.frame(t(results))[1:5]
results <- results[complete.cases(results),]
names(results) <- c("alpha","a","ge","beta","delta")
TA4B$col4[c(1,3,5,7,9)] <- signif(apply(results,2,quantile,p=0.5),2)
IQr <- signif(apply(results,2,quantile,p=c(0.25,0.75)),2)
TA4B$col4[c(2,4,6,8,10)] <- paste("(",apply(IQr,2,paste,collapse=","),")",sep="")

# TA4B - column 5 - NLS estimates
TA4B$col5[c(1,3,5,7,9)] <- as.numeric(str_sub(nls_expon[c(18,21,9,12,15),4],2)) #read alpha,a,s.ge,beta,delta
ses <- as.numeric(str_sub(nls_expon[c(19,22,10,13,16),4],3,-2))
TA4B$col5[c(2,4,6,8,10)] <- paste(rep("(",5),signif(TA4B$col5[c(1,3,5,7,9)] - 1.96*ses,3),
                                  ",",signif(TA4B$col5[c(1,3,5,7,9)] + 1.96*ses,3),rep(")",5),sep="")

# TA4B - column 6
NLS.s <- as.numeric(str_sub(nls_expon[6,4],2))/as.numeric(str_sub(nls_expon[33,4],2))
NLS.k <- as.numeric(str_sub(nls_expon[3,4],2))/as.numeric(str_sub(nls_expon[32,4],2))
NLS.g <- as.numeric(str_sub(nls_expon[24,4],2))

results <- apply(X = forecasts,MARGIN = 1,Beh.Params,s=NLS.s, 
                 k=NLS.k, g=NLS.g, spec="exp")
results <- data.frame(t(results))[1:5]
results <- results[complete.cases(results),]
names(results) <- c("alpha","a","ge","beta","delta")
TA4B$col6[c(1,3,5,7,9)] <- signif(apply(results,2,quantile,p=0.5),2)
IQr <- signif(apply(results,2,quantile,p=c(0.25,0.75)),2)
TA4B$col6[c(2,4,6,8,10)] <- paste("(",apply(IQr,2,paste,collapse=","),")",sep="")

# TA4B - column 7 - NLS estimates
TA4B$col7[c(1,3,5,7,9)] <- as.numeric(str_sub(nls_expon[c(18,21,9,12,15),5],2)) #read alpha,a,s.ge,beta,delta
ses <- as.numeric(str_sub(nls_expon[c(19,22,10,13,16),5],3,-2))
TA4B$col7[c(2,4,6,8,10)] <- paste(rep("(",5),signif(TA4B$col7[c(1,3,5,7,9)] - 1.96*ses,3),
                                  ",",signif(TA4B$col7[c(1,3,5,7,9)] + 1.96*ses,3),rep(")",5),sep="")

# TA4B - column 8
NLS.s <- as.numeric(str_sub(nls_expon[6,5],2))/as.numeric(str_sub(nls_expon[33,5],2))
NLS.k <- as.numeric(str_sub(nls_expon[3,5],2))/as.numeric(str_sub(nls_expon[32,5],2))
NLS.g <- as.numeric(str_sub(nls_expon[24,5],2))

results <- apply(X = forecasts,MARGIN = 1,Beh.Params,s=NLS.s, 
                 k=NLS.k, g=NLS.g, spec="exp")
results <- data.frame(t(results))[1:5]
results <- results[complete.cases(results),]
names(results) <- c("alpha","a","ge","beta","delta")
TA4B$col8[c(1,3,5,7,9)] <- signif(apply(results,2,quantile,p=0.5),2)
IQr <- signif(apply(results,2,quantile,p=c(0.25,0.75)),2)
TA4B$col8[c(2,4,6,8,10)] <- paste("(",apply(IQr,2,paste,collapse=","),")",sep="")

write.table(x = TA4B,file = "output/TableA4B.csv",sep = ",",row.names = FALSE , col.names=FALSE, append=FALSE)


########################
# Reproduce figures 1, 2a-b and Appendix figure 3a-b

k <- T5A$col1[3]
gamma <- T5A$col1[1]

fig_power <- expand.grid(button =seq(1500,2200,by=1),
                        MC.MB=c("MC","MB: Low Piecerate","MB: High Piecerate"),
                        stringsAsFactors = FALSE,KEEP.OUT.ATTRS = F) %>%
  mutate(value=(k*button^gamma)*(MC.MB=="MC") +
                0.04*(MC.MB == "MB: Low Piecerate") + 
           0.1*(MC.MB == "MB: High Piecerate"))

fig1 <- ggplot(fig_power, aes(x=button, y=value, group = MC.MB, linetype=MC.MB)) +  
  geom_line(size=1) +
  theme_classic() + 
  scale_y_continuous(limits=c(0,0.12),breaks=seq(0,0.12,by=0.02), expand=c(0,0)) +
  scale_x_continuous(limits=c(1500,2200),breaks=seq(1500,2200,by=100), expand=c(0,0)) +
  xlab("Button Presses e") + ylab("MB and MC") +
  scale_linetype_manual("",values=c("MB: Low Piecerate" = "dotdash",
                                  "MB: High Piecerate" = "dotted",
                                 "MC" = "solid")) +
  theme(legend.position=c(0, .55), legend.justification=c(0,0.5),
        legend.title=element_blank(),
        legend.background = element_rect(linetype="solid",color="black")) +
  ggtitle("Illustrating the Model: Marginal Benefits and Cost Curves")
  
fig1
ggsave("output/fig1.eps",fig1,width=7, height=5)
ggsave("output/fig1.pdf",fig1,width=7, height=5)


Tfig_power <- expand.grid(button =seq(1500,2200,by=1),
                         MC.MB=c("MC","MB: s","MB: s+0.01","MB: s+0.04","MB: s+0.1"),
                         stringsAsFactors = FALSE,KEEP.OUT.ATTRS = F) %>%
  mutate(value=(k*button^gamma)*(MC.MB=="MC") +
           0.01*(MC.MB == "MB: s+0.01") + 
           0.04*(MC.MB == "MB: s+0.04") + 0.1*(MC.MB == "MB: s+0.1"))

data_points <- data.frame(button = sample_means$means[1:4], 
                          MC.MB = "Observed Effort",
                          value = c(0.01,0.1,0,0.04))

library(grid)
fig2a <- ggplot(fig_power, aes(x=button, y=value, group = MC.MB, linetype=MC.MB)) +  
  geom_line(size=1) +
  geom_point(data = data_points, aes(x=button, y=value), shape=17, color="red",size=3) +
  theme_classic() + 
  scale_y_continuous(limits=c(0,0.12),breaks=seq(0,0.12,by=0.02), expand=c(0,0)) +
  scale_x_continuous(limits=c(1500,2200),breaks=seq(1500,2200,by=100), expand=c(0,0)) +
  xlab("Button Presses e") + ylab("MB and MC") +
  scale_linetype_manual("",values=c("MB: s" = "solid",
                                    "MB: s+0.01" = "longdash",
                                    "MB: s+0.04" = "dotdash",
                                    "MB: s+0.1" = "dotted",
                                    "MC" = "solid",
                                    "Observed Effort" = "1F")) +
  theme(legend.position=c(0, 0.6), legend.justification=c(0,0.6),
        legend.title=element_blank(),
        legend.background = element_rect(linetype="solid",color="black")) +
  ggtitle("Simulated Marginal Benefits and Costs with Power Cost Function") +
  annotate("text",x=1550,y=0.02,hjust=0,label="Data 0c/100") +
  annotate("text",x=1900,y=0.02,hjust=0,label="Data 1c/100") +
  annotate("text",x=2080,y=0.02,hjust=0,label="Data 4c/100") +
  annotate("text",x=1950,y=0.055,hjust=0,label="Prediction 4c/100") +
  annotate("text",x=2020,y=0.09,hjust=0,label="Data 10c/100") +
  annotate("segment", x = 1580, y = 0.015, xend = 1525, yend = 0.002, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 1980, y = 0.015, xend = 2025, yend = 0.012, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 2132, y = 0.025, xend = 2132, yend = 0.038, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 2080, y = 0.05, xend = 2115, yend = 0.042, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 2120, y = 0.095, xend = 2170, yend = 0.0995, arrow = arrow(length = unit(0.2, "cm")))

           
fig2a
ggsave("output/fig2a.eps",fig2a,width=7, height=5)
ggsave("output/fig2a.pdf",fig2a,width=7, height=5)


fig_power.b <- expand.grid(button =seq(1500,2030,by=1),
                           MC.MB=c("MC","MB: s","MB: s+0.001","MB: s+0.01"),
                           stringsAsFactors = FALSE,KEEP.OUT.ATTRS = F) %>%
  mutate(value=(k*button^gamma)*(MC.MB=="MC") +
           0.01*(MC.MB == "MB: s+0.01") + 
           0.001*(MC.MB == "MB: s+0.001"))

data_points.b <- data.frame(button = sample_means$means[c(1,3,6)], 
                            MC.MB = "Observed Effort",
                            value = c(0.01,0,0.001))

fig2b <- ggplot(fig_power.b, aes(x=button, y=value, group = MC.MB, linetype=MC.MB)) +  
  geom_line(size=1) +
  geom_point(data = data_points.b, aes(x=button, y=value), shape=17, color="red",size=3) +
  theme_classic() + 
  scale_y_continuous(limits=c(0,0.0105),breaks=seq(0,0.01,by=0.001), expand=c(0,0)) +
  scale_x_continuous(limits=c(1500,2050),breaks=seq(1500,2050,by=50), expand=c(0,0)) +
  xlab("Button Presses e") + ylab("MB and MC") +
  scale_linetype_manual("",values=c("MB: s" = "solid",
                                    "MB: s+0.001" = "dotted",
                                    "MB: s+0.01" = "longdash",
                                    "MC" = "solid",
                                    "Observed Effort" = "1F")) +
  theme(legend.position=c(0, 0.6), legend.justification=c(0,0.6),
        legend.title=element_blank(),
        legend.background = element_rect(linetype="solid",color="black")) +
  ggtitle("Simulated Marginal Benefits and Costs with Power Cost Function") +
  annotate("text",x=1550,y=0.0025,hjust=0,label="Data 0c/100") +
  annotate("text",x=1900,y=0.009,hjust=0,label="Data 1c/100") +
  annotate("text",x=1780,y=0.0025,hjust=0,label="Data 1c/1000") +
  annotate("text",x=1845,y=0.004,hjust=0,label="Prediction 1c/1000") +
  annotate("text",x=1630,y=0.0016,hjust=0,label="Crowding-out Region", color="navyblue") +
  annotate("segment", x = 1550, y = 0.002, xend = 1525, yend = 0.0002, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 1980, y = 0.0095, xend = 2025, yend = 0.0098, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 1850, y = 0.002, xend = 1880, yend = 0.0012, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 1893, y = 0.0035, xend = 1893, yend = 0.0012, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 1501, y = 0.0013, xend = 1883, yend = 0.0013, arrow = arrow(ends="both",angle=90, length = unit(0.2, "cm")),color="navyblue")
fig2b
ggsave("output/fig2b.eps",fig2b,width=7, height=5)
ggsave("output/fig2b.pdf",fig2b,width=7, height=5)


k <- T5A$col3[3]
gamma <- T5A$col3[1]

fig_expon <- expand.grid(button =seq(1500,2200,by=1),
                         MC.MB=c("MC","MB: s","MB: s+0.01","MB: s+0.04","MB: s+0.1"),
                         stringsAsFactors = FALSE,KEEP.OUT.ATTRS = F) %>%
  mutate(value=(k*exp(gamma*button))*(MC.MB=="MC") +
           0.01*(MC.MB == "MB: s+0.01") + 
           0.04*(MC.MB == "MB: s+0.04") + 0.1*(MC.MB == "MB: s+0.1"))

data_points <- data.frame(button = sample_means$means[1:4], 
                          MC.MB = "Observed Effort",
                          value = c(0.01,0.1,0,0.04))


fig3a <- ggplot(fig_expon, aes(x=button, y=value, group = MC.MB, linetype=MC.MB)) +  
  geom_line(size=1) +
  geom_point(data = data_points, aes(x=button, y=value), shape=17, color="red",size=3) +
  theme_classic() + 
  scale_y_continuous(limits=c(0,0.12),breaks=seq(0,0.12,by=0.02), expand=c(0,0)) +
  scale_x_continuous(limits=c(1500,2200),breaks=seq(1500,2200,by=100), expand=c(0,0)) +
  xlab("Button Presses e") + ylab("MB and MC") +
  scale_linetype_manual("",values=c("MB: s" = "solid",
                                    "MB: s+0.01" = "longdash",
                                    "MB: s+0.04" = "dotdash",
                                    "MB: s+0.1" = "dotted",
                                    "MC" = "solid",
                                    "Observed Effort" = "1F")) +
  theme(legend.position=c(0, 0.6), legend.justification=c(0,0.6),
        legend.title=element_blank(),
        legend.background = element_rect(linetype="solid",color="black")) +
  ggtitle("Simulated Marginal Benefits and Costs with Exponential Cost Function") +
  annotate("text",x=1550,y=0.02,hjust=0,label="Data 0c/100") +
  annotate("text",x=1900,y=0.02,hjust=0,label="Data 1c/100") +
  annotate("text",x=2080,y=0.02,hjust=0,label="Data 4c/100") +
  annotate("text",x=1950,y=0.055,hjust=0,label="Prediction 4c/100") +
  annotate("text",x=2020,y=0.09,hjust=0,label="Data 10c/100") +
  annotate("segment", x = 1580, y = 0.015, xend = 1525, yend = 0.002, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 1980, y = 0.015, xend = 2025, yend = 0.012, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 2132, y = 0.025, xend = 2132, yend = 0.038, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 2080, y = 0.05, xend = 2115, yend = 0.042, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 2120, y = 0.095, xend = 2170, yend = 0.0995, arrow = arrow(length = unit(0.2, "cm")))


fig3a
ggsave("output/figOA3a.eps",fig3a,width=7, height=5)
ggsave("output/figOA3a.pdf",fig3a,width=7, height=5)

fig_expon.b <- expand.grid(button =seq(1500,2030,by=1),
                           MC.MB=c("MC","MB: s","MB: s+0.001","MB: s+0.01"),
                           stringsAsFactors = FALSE,KEEP.OUT.ATTRS = F) %>%
  mutate(value=(k*exp(gamma*button))*(MC.MB=="MC") +
           0.01*(MC.MB == "MB: s+0.01") + 
           0.001*(MC.MB == "MB: s+0.001"))

data_points.b <- data.frame(button = sample_means$means[c(1,3,6)], 
                            MC.MB = "Observed Effort",
                            value = c(0.01,0,0.001))

fig3b <- ggplot(fig_expon.b, aes(x=button, y=value, group = MC.MB, linetype=MC.MB)) +  
  geom_line(size=1) +
  geom_point(data = data_points.b, aes(x=button, y=value), shape=17, color="red",size=3) +
  theme_classic() + 
  scale_y_continuous(limits=c(0,0.0105),breaks=seq(0,0.01,by=0.001), expand=c(0,0)) +
  scale_x_continuous(limits=c(1500,2050),breaks=seq(1500,2050,by=50), expand=c(0,0)) +
  xlab("Button Presses e") + ylab("MB and MC") +
  scale_linetype_manual("",values=c("MB: s" = "solid",
                                    "MB: s+0.001" = "dotted",
                                    "MB: s+0.01" = "longdash",
                                    "MC" = "solid",
                                    "Observed Effort" = "1F")) +
  theme(legend.position=c(0, 0.6), legend.justification=c(0,0.6),
        legend.title=element_blank(),
        legend.background = element_rect(linetype="solid",color="black")) +
  ggtitle("Simulated Marginal Benefits and Costs with Exponential Cost Function") +
  annotate("text",x=1550,y=0.0025,hjust=0,label="Data 0c/100") +
  annotate("text",x=1900,y=0.009,hjust=0,label="Data 1c/100") +
  annotate("text",x=1780,y=0.0025,hjust=0,label="Data 1c/1000") +
  annotate("text",x=1845,y=0.004,hjust=0,label="Prediction 1c/1000") +
  annotate("text",x=1630,y=0.0016,hjust=0,label="Crowding-out Region", color="navyblue") +
  annotate("segment", x = 1550, y = 0.002, xend = 1525, yend = 0.0002, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 1980, y = 0.0095, xend = 2025, yend = 0.0098, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 1850, y = 0.002, xend = 1880, yend = 0.0012, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 1893, y = 0.0035, xend = 1884, yend = 0.0012, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = 1501, y = 0.0013, xend = 1883, yend = 0.0013, arrow = arrow(ends="both",angle=90, length = unit(0.2, "cm")),color="navyblue")
fig3b
ggsave("output/figOA3b.eps",fig3b,width=7, height=5)
ggsave("output/figOA3b.pdf",fig3b,width=7, height=5)


#######################
### Generate data for histograms of individual parameters (appendix figures 9a-f)
### Later to be read by 3_ExpertParameters_Histograms.do
#######################

results <- apply(X = forecasts,MARGIN = 1,Beh.Params,s=T5A$col1[5], 
                 k=T5A$col1[3], g=T5A$col1[1], spec="power")
results <- data.frame(t(results))[c(1,2,4,5,7)]
names(results) <- c("expertAltruism","expertWarmGlow","expertBeta","expertDelta",
                    "expertLambda")
results$curv1 <- curv1[,1]
results$curv88 <- curv88[,1]
results$curv44 <- curv44[,1]

write.csv(x=results, file="dtafiles/HistogramData.csv",row.names = FALSE,na = ".")
