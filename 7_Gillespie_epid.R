# SIR -----------------------------------------------------

miGillespie<-function(tmax,beta,gama,x0,t0=0){
  v<-rbind(c(-1,1,0),c(0,-1,1))
  xm<-c(t=t0,x0)
  x<-x0
  t<-t0
  while(t<tmax){
    r1<-runif(1)
    r2<-runif(1)
    a<-c(beta*x['S']*x['I']/sum(x0),gama*x['I'])
    if(any(a>0)){
      tau<-log(1/r1)/sum(a)
      s<-min(which(cumsum(a)>=r2*sum(a)))
      t<-t+tau
      x<-x+v[s,]
      xm<-rbind(xm,c(t,x))    
    }else{
      t<-tmax+1
    }
  }
  return(xm)
}

## Ejemplo

beta<-2
gama<-0.2

N<-100
S0<-N-2
I0<-2
R0<-0
x0<-c(S=S0,I=I0,R=R0)

tmax<-30

set.seed(52); xm1<-miGillespie(tmax,beta,gama,x0)
set.seed(2); xm2<-miGillespie(tmax,beta,gama,x0)
set.seed(10); xm3<-miGillespie(tmax,beta,gama,x0)

matplot(xm1[,1],xm1[,-1],t="l",xlab="Tiempo",
        ylab="Individuos")
matplot(xm2[,1],xm2[,-1],t="l",add=TRUE)
matplot(xm3[,1],xm3[,-1],t="l",add=TRUE)
legend("topright", c("S","I","R"),col=1:3,lty=1,lwd=2)

# SIS --------------------------------------------------------

miGillespie_SIS<-function(tmax,beta,gama,x0,t0=0){
  v<-rbind(c(-1,1),c(1,-1))
  xm<-c(t=t0,x0)
  x<-x0
  t<-t0
  while(t<tmax){
    r1<-runif(1)
    r2<-runif(1)
    a<-c(beta*x['S']*x['I']/sum(x0),gama*x['I'])
    if(any(a>0)){
      tau<-log(1/r1)/sum(a)
      s<-min(which(cumsum(a)>=r2*sum(a)))
      t<-t+tau
      x<-x+v[s,]
      xm<-rbind(xm,c(t,x))    
    }else{
      t<-tmax+1
    }
  }
  return(xm)
}

## Ejemplo
beta<-2
gama<-0.5

N<-200
S0<-N-2
I0<-2
x0<-c(S=S0,I=I0)

tmax<-50

set.seed(52) 
xm1<-miGillespie_SIS(tmax,beta,gama,x0)
xm2<-miGillespie_SIS(tmax,beta,gama,x0)
xm3<-miGillespie_SIS(tmax,beta,gama,x0)

matplot(xm1[,1],xm1[,-1],t="l",xlab="Tiempo",
        ylab="Individuos")
matplot(xm2[,1],xm2[,-1],t="l",add=TRUE)
matplot(xm3[,1],xm3[,-1],t="l",add=TRUE)
legend("topright", c("S","I"),col=1:2,lty=1,lwd=2)


#punto de equilibrio del sistema determinista
#S/N=gama/beta
equil_s<-gama/beta*N
abline(h=equil_s)
equil_i<-N-equil_s
abline(h=equil_i,col="red")

pi1<-gama/(gama+beta*(equil_s/200))
pi2<-1-pi1
pi1
pi2

# Modelo Determinista y Estocástico --------------------------------------------------------

# Load deSolve library
library(deSolve)

#### Solver del modelo SIR determinista
SIR<-function(theta,t=tiempos){
  SIRmod <- function(t, x, theta) { 
    with(as.list(c(theta, x)), 
         { 
           dS <- -bet*S*I/(S+I+R) 
           dI <- bet*S*I/(S+I+R)-gam*I
           dR <- gam*I
           res <- c(dS, dI, dR) 
           list(res) }
    ) }
  ## Solver
  out <- lsoda(x0, t, SIRmod, theta)
  out[which(out[,3]<0),3]<-0
  return(out)
}

beta<-2
gama<-0.2
theta=c(bet=beta,gam=gama)

N<-500
S0<-N-2
I0<-2
R0<-0
x0<-c(S=S0,I=I0,R=R0)

tmax<-30
t<-seq(0,tmax,by=1)

res<-SIR(theta,t)

#Plot epidemic curves
par(mar=c(5,5,5,1))
matplot(res[,1],res[,-1],t="l",lwd=3,xlab="Time",ylab="Population")
legend("topright",c("S","I","R"),col=1:3,lty=1:3,lwd=3)

#simulaciones estocásticas

set.seed(52)
xm1<-miGillespie(tmax,beta,gama,x0)
xm2<-miGillespie(tmax,beta,gama,x0)
xm3<-miGillespie(tmax,beta,gama,x0)
xm4<-miGillespie(tmax,beta,gama,x0)
xm5<-miGillespie(tmax,beta,gama,x0)

matplot(xm1[,1],xm1[,-1],t="l",add=TRUE)
matplot(xm2[,1],xm2[,-1],t="l",add=TRUE)
matplot(xm3[,1],xm3[,-1],t="l",add=TRUE)
matplot(xm4[,1],xm4[,-1],t="l",add=TRUE)
matplot(xm5[,1],xm5[,-1],t="l",add=TRUE)


# Aumentamos el tamaño de la población----------------------------------

N<-5000
S0<-N-2
I0<-2
R0<-0
x0<-c(S=S0,I=I0,R=R0)

res<-SIR(theta,t)

#Plot epidemic curves
par(mar=c(5,5,5,1))
matplot(res[,1],res[,-1],t="l",lwd=3,xlab="Time",ylab="Population")
legend("topright",c("S","I","R"),col=1:3,lty=1:3,lwd=3)

#simulaciones estocásticas
set.seed(52)
xm1<-miGillespie(tmax,beta,gama,x0)
xm2<-miGillespie(tmax,beta,gama,x0)
xm3<-miGillespie(tmax,beta,gama,x0)
xm4<-miGillespie(tmax,beta,gama,x0)
xm5<-miGillespie(tmax,beta,gama,x0)

matplot(xm1[,1],xm1[,-1],t="l",add=TRUE)
matplot(xm2[,1],xm2[,-1],t="l",add=TRUE)
matplot(xm3[,1],xm3[,-1],t="l",add=TRUE)
matplot(xm4[,1],xm4[,-1],t="l",add=TRUE)
matplot(xm5[,1],xm5[,-1],t="l",add=TRUE)

# ------------------------------------------------------------------------

N<-10000
S0<-N-2
I0<-2
R0<-0
x0<-c(S=S0,I=I0,R=R0)

res<-SIR(theta,t)

#Plot epidemic curves
par(mar=c(5,5,5,1))
matplot(res[,1],res[,-1],t="l",lwd=3,xlab="Time",ylab="Population")
legend("topright",c("S","I","R"),col=1:3,lty=1:3,lwd=3)

#simulaciones estocásticas
set.seed(52)
xm1<-miGillespie(tmax,beta,gama,x0)
xm2<-miGillespie(tmax,beta,gama,x0)
xm3<-miGillespie(tmax,beta,gama,x0)

matplot(xm1[,1],xm1[,-1],t="l",add=TRUE)
matplot(xm2[,1],xm2[,-1],t="l",add=TRUE)
matplot(xm3[,1],xm3[,-1],t="l",add=TRUE)


