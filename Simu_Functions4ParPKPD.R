### Written functions for our simulation study only
## There is no covariates adjustment
## Cluster robust SE is calculated
#-----------------------------------------------------------#
#### Extended A-G model ####
#4-par PK/PD model: a=gamma>0,b=ka>0,d=C50>0, c=delta;
gfunE <- function(t, a,b,c,d) {
	f = b/(b-1)*(exp(-t)-exp(-b*t))*(t>0)
	g = a*log(d)-log(d^a+f^a)+c*(1-exp(-b*t))
	ifelse(is.na(t),0,g*(t>0))
}

#-----------------------------------------------------------#

# parametric method: specify the risk set at each distinct failure time
dicp.dose <- function(data=dataE, doseX=NULL){

## the first dataset must include variables named "id","event","startT","endT"
  id <- data$id
  event <- data$event
  startT <- data$startT
  endT <- data$endT
  sort.t <- sort(unique(endT[event==1]))

## the second dataset must include variable "treat"

  treat <- data$treat

  if(length(doseX)>=1) {
	Dmatrix=as.matrix(data[,doseX])
	} else {
		Dmatrix=matrix(0,nrow(data))
	}

  dic.list <- vector("list", length(sort.t))
  
 for (i in 1:length(sort.t)){
    f.sub <- which( (endT==sort.t[i]) & (event==1) ) 
    Rset <- which( (startT<sort.t[i]) & (endT>=sort.t[i]) )  
    
    dis.matrix <- apply(Dmatrix, 2, function(w) (sort.t[i]-w)*((sort.t[i]-w)>0)*(treat==1) )
  
    dic.list[[i]] <- list(sort.t=sort.t[i]
                          ,f.sub=f.sub
                          ,Rset=Rset
                          ,dis.matrix=dis.matrix )
}
  return(dic.list)
}#end of function

#------------------------------------------------------------------------#

# parametric method: log likelihood function 
likpd.dose <- function(parvec, data=dataE, diclst, fun=fun){
  a.iter <- exp(parvec[1])
  b.iter <- exp(parvec[2])
  c.iter <- parvec[3]
  d.iter <- exp(parvec[4])
  
  lik <- 0
  
  for (i in 1:length(diclst)){

    Rset <- diclst[[i]]$Rset
    f.sub <- diclst[[i]]$f.sub

    gt <- apply(diclst[[i]]$dis.matrix, 2,fun, a=a.iter, b=b.iter, c=c.iter, d=d.iter)
    Gt <- rowSums(gt)
       
    A <- sum(Gt[f.sub])
    B <- sum(exp(Gt[Rset]))
    lik <- lik+A-length(f.sub)*log(B)
  }
  return(-lik)  
}#end of function

#-----------------------------------------------------#

dicp.dose.robust <- function(data=dataE, doseX=NULL){

## the first dataset must include variables named "id","event","startT","endT"
  id <- data$id
  event <- data$event
  startT <- data$startT
  endT <- data$endT

  sort.t <- endT

## the second dataset must include variable "treat"
## the second dataset must include variable "doe" if season=TRUE

  treat <- data$treat

  if(length(doseX)>=1) {
	Dmatrix=as.matrix(data[,doseX])
	} else {
		Dmatrix=matrix(0,nrow(data))
	}

  dic.robust.list <- vector("list", length(sort.t))
  
  for (i in 1:length(sort.t)){
    f.sub <- i
    Rset <- which( (startT<sort.t[i]) & (endT>=sort.t[i]) )  
    Rofset <- which( (sort.t>startT[i]) & (sort.t <= endT[i]) &(event==1) )  

    dis.matrix <- apply(Dmatrix, 2, function(w) (sort.t[i]-w)*((sort.t[i]-w)>0)*(treat==1) )
  
    dic.robust.list[[i]] <- list(sort.t=sort.t[i]
                          ,f.sub=f.sub
                          ,Rset=Rset
				  ,Rofset=Rofset
                          ,dis.matrix=dis.matrix )
    }
  return(dic.robust.list)
}#end of function


#-----------------------------------------------------#

# gradient of log likelihood 
gr.pd.robust <- function(parvec, data=dataE, diclst, fun=fun){
  a.iter <- exp(parvec[1])
  b.iter <- exp(parvec[2])
  c.iter <- parvec[3]
  d.iter <- exp(parvec[4])

  U <- matrix(0, nrow=length(diclst), ncol=length(parvec))
  Gt <- matrix(0, nrow=length(diclst), ncol=nrow(data))
  dGda <- matrix(0, nrow=length(diclst), ncol=nrow(data))
  dGdb <- matrix(0, nrow=length(diclst), ncol=nrow(data))
  dGdc <- matrix(0, nrow=length(diclst), ncol=nrow(data))
  dGdd <- matrix(0, nrow=length(diclst), ncol=nrow(data))

  S0 <- rep(0,length(diclst))
  S1a <- rep(0,length(diclst))
  S1b <- rep(0,length(diclst))
  S1c <- rep(0,length(diclst))
  S1d <- rep(0,length(diclst))

  for (i in 1:length(diclst)){
    f.sub <- diclst[[i]]$f.sub
    Rset <- diclst[[i]]$Rset
 
    gt <- apply(diclst[[i]]$dis.matrix, 2, fun, a=a.iter, b=b.iter, c=c.iter, d=d.iter)
    Gt[i,] <- rowSums(gt)

    gam <- a.iter
    ka <- b.iter
    c50 <- d.iter
    delta <- c.iter

    # non-linear part
    dGda[i,] <- rowSums( apply(diclst[[i]]$dis.matrix, 2,
		function(w) ifelse(is.na(w)|w<=0, 0,  gam*(w>0)*(ka/(ka-1)*(exp(-w)-exp(-ka*w))*(w>0))^gam*
				(log(c50)-log(ka/(ka-1)*(exp(-w)-exp(-ka*w))*(w>0)))/(c50^gam+(ka/(ka-1)*(exp(-w)-exp(-ka*w))*(w>0))^gam) )))  
 
    dGdb[i,] <- rowSums( apply(diclst[[i]]$dis.matrix, 2,
		function(w) ifelse(is.na(w)|w<=0, 0, ka*(w>0)*( delta*w*exp(-ka*w)*(w>0)-gam*(ka/(ka-1)*(exp(-w)-exp(-ka*w))*(w>0))^(gam-1)/
			(c50^gam+(ka/(ka-1)*(exp(-w)-exp(-ka*w))*(w>0))^gam)*(w>0)*(((ka^2*w-ka*w+1)*exp(-ka*w)*(w>0)-exp(-w))/(ka-1)^2) )  )))

    dGdc[i,] <- rowSums( apply(diclst[[i]]$dis.matrix, 2,
		function(w) ifelse(is.na(w)|w<=0, 0, (w>0)*( 1-exp(-ka*w) ))))

    dGdd[i,] <- rowSums( apply(diclst[[i]]$dis.matrix, 2,
		function(w) ifelse(is.na(w)|w<=0, 0, c50*(w>0)*( gam/c50-gam*c50^(gam-1)/(c50^gam+(ka/(ka-1)*(exp(-w)-exp(-ka*w))*(w>0))^gam) ))))


    S0[i] = sum(exp(Gt[i,Rset]))
    S1a[i] = sum(exp(Gt[i,Rset])*dGda[i,Rset])
    S1b[i] = sum(exp(Gt[i,Rset])*dGdb[i,Rset])
    S1c[i] = sum(exp(Gt[i,Rset])*dGdc[i,Rset])
    S1d[i] = sum(exp(Gt[i,Rset])*dGdd[i,Rset])

    U[i, 1] <- dGda[i,f.sub]-S1a[i]/S0[i]
    U[i, 2] <- dGdb[i,f.sub]-S1b[i]/S0[i]
    U[i, 3] <- dGdc[i,f.sub]-S1c[i]/S0[i]
    U[i, 4] <- dGdd[i,f.sub]-S1d[i]/S0[i]
    }

 V <- matrix(0, nrow=length(diclst), ncol=length(parvec))

 for (i in 1:length(diclst)){

    f.sub <- diclst[[i]]$f.sub
    Rofset <- diclst[[i]]$Rofset

    Vk <- rep(0,4)

    wt = exp(Gt[Rofset,f.sub])/S0[Rofset]

    Vk[1] <- sum((dGda[Rofset,f.sub]-S1a[Rofset]/S0[Rofset])*wt)
    Vk[2] <- sum((dGdb[Rofset,f.sub]-S1b[Rofset]/S0[Rofset])*wt)
    Vk[3] <- sum((dGdc[Rofset,f.sub]-S1c[Rofset]/S0[Rofset])*wt)
    Vk[4] <- sum((dGdd[Rofset,f.sub]-S1d[Rofset]/S0[Rofset])*wt)

    V[i, 1] <- U[i,1]-Vk[1]
    V[i, 2] <- U[i,2]-Vk[2]
    V[i, 3] <- U[i,3]-Vk[3]
    V[i, 4] <- U[i,4]-Vk[4]
  }
 V
}#end of function


#-----------------------------------------------------------#
Robust.error <- function(coefficients, Hess, data=dataE, doseX=NULL,fun=fun) {

  dic.pd.plus <- dicp.dose.robust(data=data,doseX=doseX)
  U <- gr.pd.robust(coefficients,data=data,diclst=dic.pd.plus,fun=fun)

  #gradient at estimators
  score <- colSums(U[data$event==1,])

  dataI <- data[!duplicated(data$id),]

  ## Cluster varaince matrix

  clusterM <- matrix(0, nrow=length(coefficients), ncol=length(coefficients))
  for(i in 1:nrow(dataI)) {
    pos <- which(data$id==dataI$id[i]&data$event==1 )
    if(length(pos)==1) 
	  clusterM = clusterM + (U[pos,])%*%t(U[pos,])
    else
	  clusterM = clusterM +(colSums(U[pos,]))%*%t(colSums(U[pos,]))
   }

  clusterH <- solve(Hess)%*%clusterM%*%solve(Hess)

  return(list(Score=score,ClusterH=clusterH))
}#end of function

#-----------------------------------------------------------#

peak.information <- function(parvec,doseN) {

  tau = 3*doseN+2
  a <- exp(parvec[1])
  b <- exp(parvec[2])
  c <- parvec[3]
  d <- exp(parvec[4])

  tp = -log(b)/(1-b)


 findroot <- function(t) {
	f = b/(b-1)*(exp(-t)-exp(-b*t))
	-a*f^(a-1)*b/(b-1)*(b*exp(-b*t)-exp(-t))/(d^a+f^a)+c*b*exp(-b*t)
	}
	
 #estimated gfun.Est(t):
 gfun.Est <- function(t)  gfunE(t, a=a, b=b, c=c, d=d)
 halfroot <- function(t) 1+exp(gfun.Est(tstar))-2*exp(gfun.Est(t))
 integrand <- function(t) 1-exp(gfun.Est(t))

  if(findroot(tp/4)*findroot(tp*4)<0) {
  	tstar = uniroot(findroot, c(tp/4,tp*4))$root
   } else if(findroot(tp/2)*findroot(tp*2)<0) {
	tstar = uniroot(findroot, c(tp/2,tp*2))$root
   } else if(findroot(0.0001)*findroot(tp*2)<0) {
	tstar = uniroot(findroot, c(0.0001,tp*2))$root
   } else tstar=0

  #halfroot(tau) shou be nagative!
  if( halfroot(tstar)*halfroot(tau)<0 ) {
     thalf = uniroot(halfroot, c(tstar,tau))$root
    } else thalf = tau

  peakPE = 1-exp(gfun.Est(tstar))
  area = integrate(integrand, lower = 0, upper = 3)$value

  as.numeric(c(tstar,peakPE,thalf,area))
}#end of function



