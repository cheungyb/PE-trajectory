### R version: R 4.2.1
rm(list=ls())

### specify the Workspace
setwd("C:/nBox/Stat/Nonlinear/PKPD/Manuscript/resubmit/Github")

### load functions used in the simulation study
source("Simu_Program5ParPKPD.R")
source("Simu_Functions4ParPKPD.R")

# Defines parameters values (those varying)

library(snowfall)

par1 <- c(1000, 3, 2)   #here ka=3, k=2
par2 <- c(1000, 9, 0.5) #here ka=9; k=0.5
par3 <- c(1400, 3, 2)  
par4 <- c(1400, 9, 0.5)
param=rbind(par1,par2,par3,par4)

g = 1:4
nproc <- 4

for (dose in 1:2) {

    sfInit(parallel = T, cpus = nproc)
    sfExportAll()
    
  run.simul <- function(valeur) {
      
      for (k in valeur) {

        set.seed(20230727)

        write.table(t(c(try(simu_gfun5parPKPD(n=param[k,1],ka=param[k,2],k=param[k,3],dose=dose), TRUE))), 
          paste("./work/simu_gfun5parPKPD_hete_n=", param[k,1], "_dose=", dose,"_ka=", param[k,2], ".q", sep = ""), 
		append = F, quote = F, col.names = F, row.names = F)
        
        for (m in 1:499) {         
        write.table(t(c(try(simu_gfun5parPKPD(n=param[k,1],ka=param[k,2],k=param[k,3],dose=dose), TRUE))), 
          paste("./work/simu_gfun5parPKPD_hete_n=", param[k,1], "_dose=", dose,"_ka=", param[k,2], ".q", sep = ""), 
		append = T, quote = F, col.names = F, row.names = F)                  
        } #end of m
  	 }#end of k
	}#end of run.simul 
    
   Big.simul <- sfClusterApplyLB(g, run.simul)
    
   sfStop()

}#end of for-loop
###end of snowfall

##------------- show results -----------------------##

for (dose in 1:2) {
  for (k in g) {
	results <- read.table( paste("./work/simu_gfun5parPKPD_hete_n=", param[k,1], "_dose=", dose,  
				"_ka=", param[k,2], ".q", sep = ""), header = F, fill = TRUE)
	print(dim(results))
	names(results)[1:4] = c("Size","N0","N1","Dose")
	names(results)[5:9] = c("true_gamma","true_ka","true_delta","true_C50","true_k")
	names(results)[10:15] = c("ln_gamma","ln_ka","delta","ln_C50","AIC","BIC")
	names(results)[16:19] = c("T.star","peakPE","T.half","AUC")
	
	#here time in days
	results[,c(16,18)] <- results[,c(16,18)]*30
	print(colMeans(results[,16:19]))

	# written as a csv file
	write.csv(results,file = paste("./work/simu_gfun5parPKPD_hete_csv_n=", param[k,1], "_dose=", dose,
				"_ka=", param[k,2], ".csv", sep = ""), row.names=F)
	}
}#end!

