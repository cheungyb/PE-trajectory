### R version: R 4.2.1
rm(list=ls())

### specify the Workspace
setwd("C:/nBox/Stat/Nonlinear/PKPD/Manuscript/resubmit/Github")

### load functions used in the simulation study
source("Simu_Program4ParPKPD.R")
source("Simu_Functions4ParPKPD.R")

# Defines parameters values (those varying)

library(snowfall)

par1 <- c(1000, 9, 3)
par2 <- c(1000, 3, 3)
par3 <- c(1400, 9, 3)
par4 <- c(1400, 3, 3)
param=rbind(par1,par2,par3,par4)

g = 1:4
nproc <- 4

for (dose in 1:3) {

    sfInit(parallel = T, cpus = nproc)
    sfExportAll()
    
	run.simul <- function(valeur) {
      
      for (k in valeur) {

        set.seed(20230727)

        write.table(t(c(try(simu_gfun4parPKPD(n=param[k,1],ka=param[k,2],dose=dose), TRUE))), 
          paste("./work/simu_gfun4parPKPD_hete_n=", param[k,1], "_dose=", dose,"_ka=", param[k,2], ".q", sep = ""), 
		append = F, quote = F, col.names = F, row.names = F)
        
        for (m in 1:499) {         
        write.table(t(c(try(simu_gfun4parPKPD(n=param[k,1],ka=param[k,2],dose=dose), TRUE))), 
          paste("./work/simu_gfun4parPKPD_hete_n=", param[k,1], "_dose=", dose,"_ka=", param[k,2], ".q", sep = ""), 
		append = T, quote = F, col.names = F, row.names = F)                  
        }#end of m
  	  }#end of k
	}#end of run.simul 
    
    Big.simul <- sfClusterApplyLB(g, run.simul)
    
    sfStop()
}#end of for-loop
###end of snowfall

##------------- show results -----------------------##

for (dose in 1:3) {
  for (k in g) {
	results <- read.table( paste("./work/simu_gfun4parPKPD_hete_n=", param[k,1], "_dose=", dose,  
				"_ka=", param[k,2], ".q", sep = ""), header = F, fill = TRUE)
	print(dim(results))

	names(results)[1:4] = c("Size","N0","N1","Dose")
	names(results)[5:8] = c("true_lngamma","true_lnka","true_delta","true_lnC50")
	names(results)[9:12] = c("ln_gamma","ln_ka","delta","ln_C50")
	names(results)[13:16] = c("seN_lngamma","seN_lnka","seN_delta","seN_lnC50")
	names(results)[17:22] = c("se_lngamma","se_lnka","se_delta","se_lnC50","AIC","BIC")
	names(results)[23:26] = c("T.star","peakPE","T.half","AUC")
	# relative bias in %
	Bias <- round(colMeans(results[,9:12]/results[,5:8]-1)*100,1)
	# ASE/ESD
	R.se <- round(colMeans(results[,17:20])/apply(results[,9:12], 2, sd),2)
	# CP in %
	cover <- ifelse(abs(results[,9:12]-results[,5:8])/results[,17:20]<= qnorm(1-0.05/2),1,0)
	CP <- round(colMeans(cover)*100,1)
	# RMSE
	RMSE <- round(sqrt(colMeans((results[,9:12]-results[,5:8])^2)),3)
	print(cbind(Bias, R.se, CP, RMSE))	
	
	# here time in days
	results[,c(23,25)] <- results[,c(23,25)]*30
	print(colMeans(results[,23:26]))
	cat("\n")
	
	# written as a csv file
	write.csv(results,file = paste("./work/simu_gfun4parPKPD_hete_csv_n=", param[k,1], "_dose=", dose,
				"_ka=", param[k,2], ".csv", sep = ""), row.names=F)
	}
}#end!

