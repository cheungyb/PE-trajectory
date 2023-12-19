### Written functions
#-----------------------------------------------------------#
#### Extended A-G model with no-linear function g(t) with c50>0, ka>0 but ka!=1, and gamma>0 ####
#### g(t): PK/PD model
gfun <- function(t, par) {
	gam = par[1]
	ka = par[2]
	c50 = par[3]
	delta = par[4]
	ct = ka/(ka-1)*(exp(-t)-exp(-ka*t))*(t>0)
	g = gam*log(c50)-log(c50^gam+ct^gam)+delta*(1-exp(-ka*t))
	ifelse(is.na(t),0, g*(t>0))
	}

#-----------------------------------------------------------#
# parametric method: specify the risk set at each distinct failure time
dicp.dose <- function(data=data, doseX=NULL){

## the first dataset must include variables named "id","event","startT","endT"
  id <- data$id
  event <- data$event
  startT <- data$startT
  endT <- data$endT
  sort.t <- sort(unique(endT[event==1]))

## the second dataset must include variable "treat"

  treat <- data$treat

  if(length(doseX)>1) {
	Dmatrix = as.matrix(data[,doseX])
	} else {
		Dmatrix = matrix(0,nrow(data))
	}

  dic.list <- vector("list", length(sort.t))
  for (i in 1:length(sort.t)){
    f.sub <- which( (endT==sort.t[i]) & (event==1) ) 
    Rset <- which( (startT<sort.t[i]) & (endT>=sort.t[i]) )  
    
    dis.matrix <- apply(Dmatrix, 2, function(w) (sort.t[i]-w)*((sort.t[i]-w)>0)*(treat==1) )
  
    dic.list[[i]] <- list(sort.t=sort.t[i],f.sub=f.sub
                          ,Rset=Rset,dis.matrix=dis.matrix )
	}
  return(dic.list)
}#end of function

#------------------------------------------------------------------------#
# parametric method: log likelihood function 
likpd.dose <- function(parvec, fX=NULL, data=dataE, diclst, fun=gfun){
  par.iter <- exp(parvec)
  par.iter[4] <- parvec[4]
  beta.iter <- parvec[-c(1:4)]
  
  if(length(fX)>1) {
  	X0 <- data[, fX]
  	for (j in 1:length(fX)){
  		if(is.factor(X0[,j])) {
			dX <- data.frame(X0[,j])
  			repX <- model.matrix(~ X0[,j], dX)[,-1]
			X0 <- cbind(X0[,-j], repX)
		}
	}
  }

  X0 <- as.matrix(X0)
  lik <- 0

  for (i in 1:length(diclst)){

    Rset <- diclst[[i]]$Rset
    f.sub <- diclst[[i]]$f.sub

    EX <-  X0 %*% beta.iter

    gt = apply(diclst[[i]]$dis.matrix, 2, fun, par=par.iter)
    Gt <- rowSums(gt)
       
    A <- sum(EX[f.sub]+Gt[f.sub])
    B <- sum(exp(EX[Rset]+Gt[Rset]))
    lik <- lik+A-length(f.sub)*log(B)
  }
  return(-lik)  
}#end of function

#-----------------------------------------------------#

# parametric method: specify the risk set at each observed time
dicp.dose.robust <- function(data=dataE, doseX=NULL){

## the first dataset must include variables named "id","event","startT","endT"
  id <- data$id
  event <- data$event
  startT <- data$startT
  endT <- data$endT

  tt <- endT

## the second dataset must include variable "treat"

  treat <- data$treat

  if(length(doseX)>0) {
	Dmatrix=as.matrix(data[,doseX])
	} else {
		Dmatrix=matrix(0,nrow(data))
	}

  dic.robust.list <- vector("list", length(tt))
  
  for (i in 1:length(tt)){
    f.sub <- i
    Rset <- which( (startT<tt[i]) & (endT>=tt[i]) )  
    Rofset <- which( (tt>startT[i]) & (tt <= endT[i]) &(event==1) )  
    
    dis.matrix <- apply(Dmatrix, 2, function(w) (tt[i]-w)*((tt[i]-w)>0)*(treat==1) )
  
    dic.robust.list[[i]] <- list(f.sub=f.sub
                          ,Rset=Rset
				  ,Rofset=Rofset
                          ,dis.matrix=dis.matrix )
	}
  return(dic.robust.list)
}#end of function

#-----------------------------------------------------#
# gradient of log likelihood 
gr.pd.robust <- function(parvec, fX=NULL, data, diclst, fun=gfun){
  par.iter <- exp(parvec)
  par.iter[4] <- parvec[4]

  gam <- par.iter[1]
  ka <- par.iter[2]
  c50 <- par.iter[3]
  delta <- par.iter[4]
  beta.iter <- parvec[-c(1:4)]

  if(length(fX)>1) {
  	X0 <- data[, fX]
  	for (j in 1:length(fX)){
  		if(is.factor(X0[,j])) {
			dX <- data.frame(X0[,j])
  			repX <- model.matrix(~ X0[,j], dX)[,-1]
			X0 <- cbind(X0[,-j], repX)
		}
	}
  }

  X0 <- as.matrix(X0)
  U <- matrix(0, nrow=length(diclst), ncol=length(parvec))
  Gt <- matrix(0, nrow=length(diclst), ncol=nrow(data))

  part1 <- matrix(0, nrow=length(diclst), ncol=length(beta.iter))
  part2 <- matrix(0, nrow=length(diclst), ncol=length(beta.iter))

  dGd.gam <- matrix(0, nrow=length(diclst), ncol=nrow(data))
  dGd.ka <- matrix(0, nrow=length(diclst), ncol=nrow(data))
  dGd.c50 <- matrix(0, nrow=length(diclst), ncol=nrow(data))
  dGd.delta <- matrix(0, nrow=length(diclst), ncol=nrow(data))

  S0 <- rep(0,length(diclst))
  S1.gam <- rep(0,length(diclst))
  S1.ka <- rep(0,length(diclst))
  S1.c50 <- rep(0,length(diclst))
  S1.delta <- rep(0,length(diclst))

  for (i in 1:length(diclst)){
    f.sub <- diclst[[i]]$f.sub
    Rset <- diclst[[i]]$Rset

    X <-  X0

    gt = apply(diclst[[i]]$dis.matrix, 2, fun, par=par.iter)
    Gt[i,] <- rowSums(gt)
    Xt <- as.vector(as.matrix(X[Rset,])%*% beta.iter)

    # linear part: vector
    S0[i] <- sum(exp(Xt+as.vector(Gt[i,Rset])))
    part1[i,] <- X[f.sub, ]

    if(length(Rset)==1 | length(beta.iter)==1)
   	 part2[i,] <- sum(exp(Xt+as.vector(Gt[i,Rset]))*X[Rset, ])
    else
   	 part2[i,] <- colSums(exp(Xt+as.vector(Gt[i,Rset]))*X[Rset, ])

    if(ncol(U)>4)    U[i, 5:ncol(U)] <- part1[i,]-part2[i,]/S0[i]

    # gradiant of non-linear part corresponding to fun=gfun:
    dGd.gam[i,] <- rowSums( apply(diclst[[i]]$dis.matrix, 2,
		function(w) ifelse(is.na(w)|w<=0, 0,  gam*((ka/(ka-1)*(exp(-w)-exp(-ka*w)))^gam*
			(log(c50)-log(ka/(ka-1)*(exp(-w)-exp(-ka*w))))/(c50^gam+(ka/(ka-1)*(exp(-w)-exp(-ka*w)))^gam)) )))  
 
    dGd.ka[i,] <- rowSums( apply(diclst[[i]]$dis.matrix, 2,
		function(w) ifelse(is.na(w)|w<=0, 0, ka*((delta*w*exp(-ka*w)-gam*(ka/(ka-1)*(exp(-w)-exp(-ka*w))*(w>0))^(gam-1)/
			(c50^gam+(ka/(ka-1)*(exp(-w)-exp(-ka*w))*(w>0))^gam)*(((ka^2*w-ka*w+1)*exp(-ka*w)-exp(-w))/(ka-1)^2)) )  )))

    dGd.c50[i,] <- rowSums( apply(diclst[[i]]$dis.matrix, 2,
		function(w) ifelse(is.na(w)|w<=0, 0, c50*((gam/c50-gam*c50^(gam-1)/(c50^gam+(ka/(ka-1)*(exp(-w)-exp(-ka*w))*(w>0))^gam)) ))))

    dGd.delta[i,] <- rowSums( apply(diclst[[i]]$dis.matrix, 2,
		function(w) ifelse(is.na(w)|w<=0, 0, (1-exp(-ka*w)) )))


    S1.gam[i] <- sum(exp(Xt+as.vector(Gt[i,Rset]))*dGd.gam[i,Rset])
    S1.ka[i] <- sum(exp(Xt+as.vector(Gt[i,Rset]))*dGd.ka[i,Rset])
    S1.c50[i] <- sum(exp(Xt+as.vector(Gt[i,Rset]))*dGd.c50[i,Rset])
    S1.delta[i] <- sum(exp(Xt+as.vector(Gt[i,Rset]))*dGd.delta[i,Rset])

    U[i, 1] <- dGd.gam[i,f.sub]-S1.gam[i]/S0[i]
    U[i, 2] <- dGd.ka[i,f.sub]-S1.ka[i]/S0[i]
    U[i, 3] <- dGd.c50[i,f.sub]-S1.c50[i]/S0[i]
    U[i, 4] <- dGd.delta[i,f.sub]-S1.delta[i]/S0[i]
  }

 V <- matrix(0, nrow=length(diclst), ncol=length(parvec))

 for (i in 1:length(diclst)){

    f.sub <- diclst[[i]]$f.sub
    Rofset <- diclst[[i]]$Rofset

    X <-  X0

    Xt <- as.vector(as.matrix(part1[Rofset,])%*%beta.iter)

    Vk <- rep(0,4)

    wt = exp(Xt+as.vector(Gt[Rofset,f.sub]))/S0[Rofset]

    if(ncol(V)>4 )    {
		V[i, 5:ncol(V)] <- U[i, 5:ncol(U)]-sum((part1[Rofset,]-part2[Rofset,]/S0[Rofset])*wt)
	}

    Vk[1] <- sum((dGd.gam[Rofset,f.sub]-S1.gam[Rofset]/S0[Rofset])*wt)
    Vk[2] <- sum((dGd.ka[Rofset,f.sub]-S1.ka[Rofset]/S0[Rofset])*wt)
    Vk[3] <- sum((dGd.c50[Rofset,f.sub]-S1.c50[Rofset]/S0[Rofset])*wt)
    Vk[4] <- sum((dGd.delta[Rofset,f.sub]-S1.delta[Rofset]/S0[Rofset])*wt)

    V[i, 1] <- U[i,1]-Vk[1]
    V[i, 2] <- U[i,2]-Vk[2]
    V[i, 3] <- U[i,3]-Vk[3]
    V[i, 4] <- U[i,4]-Vk[4]
  }
 V
}#end of function


#-----------------------------------------------------------#
Robust.error <- function(coefficients, Hess, clusterID=NULL, data=dataE, fX=NULL, doseX=NULL, fun=gfun) {

  dic.pd.plus <- dicp.dose.robust(data=data,doseX=doseX)

  U <- gr.pd.robust(coefficients,fX=fX,data=data,diclst=dic.pd.plus,fun=fun)

  #gradient at estimators
  score <- colSums(U[data$event==1,])

  if(is.null(clusterID)) 
	clusterID = unique(data$id)
  else
	clusterID = unique(clusterID)

  ## Cluster varaince matrix
  clusterM <- matrix(0, nrow=length(coefficients), ncol=length(coefficients))
  for(i in 1:length(clusterID)) {
    pos <- which(data$id==clusterID[i] & data$event==1)
    if(length(pos)==1) 
	  clusterM =clusterM + (U[pos,])%*%t(U[pos,])
    else
	  clusterM =clusterM +(colSums(U[pos,]))%*%t(colSums(U[pos,]))
   }

  clusterH <- solve(Hess)%*%clusterM%*%solve(Hess)

  return(list(Score=score, ClusterH=clusterH))
}#end of function

#-----------------------------------------------------------#

peak.information <- function(parvec,tau) {

  parvec = as.numeric(parvec)
  par = exp(parvec)
  par[4] = parvec[4]

  gam = par[1]
  ka = par[2]
  c50 = par[3]
  delta = par[4]

  tp = -log(ka)/(1-ka)

findroot <- function(t) {
	ct = ka/(ka-1)*(exp(-t)-exp(-ka*t))
	-gam*ct^(gam-1)*ka/(ka-1)*(ka*exp(-ka*t)-exp(-t))/(c50^gam+ct^gam)+delta*ka*exp(-ka*t)
	}
	
#estimated gfun.Est(t):
gfun.Est <- function(t)  gfun(t, par)
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


