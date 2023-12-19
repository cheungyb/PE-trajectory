### Estimation4ParPKPD.R
### R version: R 4.2.1
### R code for estimation of a 4-parameter PK/PD model for analysis of recurrent events 
### after a single dose or multiple doses of a treatment/exposure 
### Read the comments below and Example_data_2023.csv for the data format expected

rm(list=ls())

### specify the Workspace
setwd("C:/nBox/Stat/Nonlinear/PKPD/Manuscript/resubmit/Github")

### load functions used in the analysis
source("Functions4ParPKPD.R")

### Data format: one person has multiple rows with the same "id".
### Each record ends with "event" = 1 if the outcome event occurs; otherwise "event" = 0 (censored).
### "treat" is the binary tretment/exposure variable whose effect is to be represented by the 4-parameter function.
### "doe" is the start date of follow-up, in yyyy-mm-dd format.
### "startT" and "endT" are time in months since doe; they define the start and end of each spell of follow-up/at-risk time.
### The data file must include the variables "id", "treat", "startT", "endT" and "event". 
### Covariates can be time-varying (e.g. current age)
### d1, d2, ... are variable names for timing of doses; time is months since doe
### doseX = NULL means a single dose at time "0"; equivalently, specify doseX = c("d1") where d1=0
### Missing doses are represented by a blank cell in the csv file

### Read example data from a csv file
df_E <- read.csv("Example_data_2023.csv",header=TRUE)
df_E[1:10,]

### the nonlinear 4-parameter PKPD function: 
### g(t)=gamma*log(C50)-log(C50^gamma+(Ka/(Ka-1)*(exp(-t)-exp(-Ka*t)))^gamma)+delta*(1-exp(-Ka*t))
### with gamma, C50>0, Ka>0 and Ka!=1 for representing the treatment effect at time t since time 0 (first dose)
### if there is multiple doses, then G(t)=g(t-d1)+g(t-d2)+g(t-d3)+...
### if someone hadn't received the scheduled second dose, i.e. d2=NA, then g(t-d2)=0 for t>d2

### specify the date of doses of intervention
### E.g., in the example file, "d1","d2","d3","d4"
# single dose: doseX = NULL or doseX = c("d1")
# if four doses are specified
doseX = c("d1","d2","d3","d4")

### specify the covariates included in the regression
### time-contant covariate: "gender", "village"
### time-varying covariate: "agegr"
### data has been splitted at the time point where the time-varying covariate changed
X = c("gender", "village", "agegr")

#Convert "village" to a factor, as "village" is a categorical covariate 
df_E$village = factor(df_E$village)
levels(df_E$village) 
#here 4 levels, "village1" is the reference group 

#Finally, there are 5 estimators of coefficients for covariates
Covariates = c("Gender","Village2","Village3","Village4","Agegr")

### Those having incomplete covariates information are excluded in the analysis
df_E <- df_E[complete.cases(df_E[,X]),]

### After the estimation, we can specify the choice of SE
### Naive SE or Cluster Robust SE 

#---------------------------- The Model ----------------------------------#

### specify the recurrent event data
dic.pd <- dicp.dose(data=df_E, doseX=doseX)

### optionally, specify initial values in optimization
q <- length(Covariates)
par0 <- c(rep(0.1,3),rep(0,q+1))

### use the "optim" command to minimize the negative log-likelihood function
fit.pd <- optim(par0, likpd.dose, fX=X, data=df_E
                ,diclst=dic.pd, fun=gfun, method="BFGS"
                ,hessian=TRUE, control=list(maxit=1000, trace=TRUE))
fit.pd$par

### compute loglikelihood value 
lik <- -likpd.dose(fit.pd$par, fX=X, data=df_E, diclst=dic.pd, fun=gfun)
cat("log likelihood= ", lik, "\n")

### show AIC;BIC 
AIC <- -2*lik+2*length(fit.pd$par)
BIC <- -2*lik+log(nrow(df_E[df_E$event==1,]))*length(fit.pd$par)
cat("AIC= ", AIC, "\n")
cat("BIC= ", BIC, "\n")

### save estimation results
### constraints: gamma, ka and c50 are non-negative values
gamma <- exp(fit.pd$par[1])
ka <- exp(fit.pd$par[2])
c50 <- exp(fit.pd$par[3])
delta <- fit.pd$par[4]
H <- fit.pd$hessian

### summary results
coef.pd <- c(gamma,ka,c50,delta,fit.pd$par[-c(1:4)])
#Naive SE
#SEs <- sqrt(diag(solve(H)))


#----------------- Calculate the clustered standard errors (SE) ---------------------#
### Cluster robust SE with allowing for clustering of observations/events within person 

clusterID <- unique(df_E$id)
error <- Robust.error(fit.pd$par, Hess=H, clusterID=clusterID, data=df_E, fX=X, doseX=doseX, fun=gfun)

## Robust SE allowing for multiple events within cluster/person
SEs <- sqrt(diag(error$ClusterH))

### Delta method are used for computing SEs of (gamma,ka,c50) 
se.pd <- numeric(length(fit.pd$par))
se.pd[1:3] <- c(gamma,ka,c50)* SEs[1:3]
se.pd[-c(1:3)] <- SEs[-c(1:3)]

z.pd <- coef.pd/se.pd
tab.pd <- cbind(coef.pd, se.pd, z.pd)

### labels the outputs
rownames(tab.pd)= c("gamma","Ka","C50","delta",Covariates)

### show results
print(round(tab.pd,4))

### ------- Peak information ----------- ###
Infor <- peak.information(fit.pd$par[1:4], tau=24)
cat("Time in month to peak PE = ", Infor[1], "\n")
cat("Peak PE = ", Infor[2], "\n")
cat("Time in month to PE/2 after reaching peak PE  = ", Infor[3], "\n")
cat("AUC = ", Infor[4], "\n")

#The end!

