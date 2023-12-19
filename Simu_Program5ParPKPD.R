#gfunE: 4-par PKPD model as the analysis model

#### Extended A-G model
## simulation function 

simu_gfun5parPKPD <- function(n=1000, dose=1, ka=9, gamma=3, k=2, c50=0.4, delta=0.1, sig=20){

z <- rbinom(n,1,0.5) 
#Tau in month
Tau = 2*dose+3

# funcitons for generation of recurrent event times only

gt <- function(t)  {
	C=(ka/(ka-k))*(exp(-k*t)-exp(-ka*t))*(t>0) 
	g=gamma*log(c50)-log(c50^gamma+C^gamma)+delta*(1-exp(-ka*t))
	ifelse(is.na(t),0,g*(t>0))
	}

d1 <- rep(0,n)
d2 <- runif(n,1,2)
d3 <- runif(n,3,4)

#single dose; multi-dose (up to three doses)
Gt <- function(t, dose) {
	(dose==1)*gt(t-d1)+(gt(t-d1)+gt(t-d2))*(dose==2)+(gt(t-d1)+gt(t-d2)+gt(t-d3))*(dose==3)

}

Ctime <- c(rep(1, 0.8*n),runif(0.2*n,0.8,1))*Tau

w <- rgamma(n,shape=sig,scale=1/sig)

lambda0 = 0.2
lambda.ft <- function(t) lambda0*exp(z*Gt(t,dose=dose))*w
lambda = 0.28*w

# choose a value of lambda such that lambda.ft <= lambda for any value of t

data <- NULL
for( i in 1:n) {

	T = vector()
	Tstar = 0
	j = 0

	while (Tstar < Ctime[i]) {
		R <- rexp(1,lambda[i])
		Tstar = Tstar+R
		v <- runif(1)
	if ((lambda.ft(rep(Tstar,n))[i] >= v*lambda[i]) & Tstar < Ctime[i]){
		j = j+1
		T[j] = Tstar
		}
	}

id <- rep(i,j+1)
Z <- rep(z[i],j+1)
D1 <- rep(d1[i],j+1)
D2 <- rep(d2[i],j+1)
D3 <- rep(d3[i],j+1)

start <- c(0,T)
end <- c(T,Ctime[i])
status <- c(rep(1,j),0)
add <- data.frame(cbind(id=id,start=start,end=end,status=status,Z=Z,D1=D1,D2=D2,D3=D3))
data <- rbind(data, add)

}

N1 = nrow(data[data$Z==1&data$status==1,])
N0 = nrow(data[data$Z==0&data$status==1,])

df_E = data
names(df_E)[c(2:5)] <- c("startT", "endT","event","treat")
names(df_E)[6:8] <- c("d1","d2","d3")

 
####====================== model fit ======================####

doseX=NULL 
if(dose==2)  doseX=c("d1","d2")
if(dose==3)  doseX=c("d1","d2","d3")

dic.pd <- dicp.dose(data=df_E,doseX=doseX)

### specify initial values in optimization
ln_a0 <- 0.1
ln_b <- 0.1
c0 <- 0
ln_d <- 0.1
par0 <- c(ln_a0, ln_b, c0, ln_d)

## there is no covariates adjustment in the program
### use the "optim" command to minimize the negative log-likelihood function 
### no SE is required: hessian=FALSE

fit.pd <- optim(par0, likpd.dose, data=df_E, diclst=dic.pd
                ,fun=gfunE, method="BFGS", hessian=FALSE
		    ,control=list(maxit=1000, trace=FALSE))

est = fit.pd$par

### compute loglikelihood value 
lik = -likpd.dose(fit.pd$par,data=df_E,diclst=dic.pd,fun=gfunE)

### show AIC, BIC
aic = -2*lik+2*length(fit.pd$par)
bic = -2*lik+log(nrow(df_E[df_E$event==1,]))*length(fit.pd$par)

## peak information
peakV <- peak.information(fit.pd$par,doseN=dose)

## true values
trueV <- c(dose,gamma,ka,delta,c50,k)

#return all values
c(n,N0,N1,trueV,est,aic,bic,peakV)

}# end of simu function