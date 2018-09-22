library(rootSolve)
library(manipulate)
library(deSolve)
library(pracma)
library(FME)
library(gdata)

#First we evaluate each species in it's isolated environment.
data=read.csv(choose.files())
tR=data$t..days
Rhiz=data$Rhizopertha
x=Rhiz[1]

data=read.csv(choose.files())
tO=data$t..days.
Oryz=data$Oryzaephilus
y=Oryz[1]

#We we will use non-linear regression methods to develop a model for
#each species.To more easily compare the data with the model, 
#the time intervals given in the data are used.

#First, for the Rhizopertha population
R=nls(Rhiz~(Kr*x)/(x+(Kr-x)*exp(-Rr*tR)), start=list(Rr=.05, Kr=335))
yr=predict(R, list(S=tR))

#We can then plot the data vs time for the isolated rhizopertha
#population along the model.
plot(tR,Rhiz,ylab="Rhizopertha Population",xlab="time (days)",main="Isolated Rhizopertha Population")
points(tR, yr, col="Red")
legend("topleft",c("Data","Rhiz Model"),pch=19,col=c("black","red"))

#Then, for the Oryzaephilus population.
O=nls(Oryz~(Ko*y)/(y+(Ko-y)*exp(-Ro*tO)), start=list(Ro=.05, Ko=480))
yo=predict(O, list(S=tO))

#We can then plot the model and data for the Oryzaephilus population.
plot(tO,Oryz,main="Isolated Oryzaephilus Population",ylab="Oryzaephilus Population",xlab="time (days)",ylim=c(0,500))
points(tO, yo, col="blue")
legend("topleft",c("Data","Oryz Model"),pch=19,col=c("black","blue"))

#We can now define the intrinsic growth rate and carrying capacity
#the Rhizopertha population,
CR=coef(R)
Rr=CR[1]
Kr=CR[2]

#and the Oryzaephilus population
CO=coef(O)
Ro=CO[1]
Ko=CO[2]

#################
#Now we can develop a model for the competitive scenario based 
#off of these values.
data=read.csv(choose.files())
x1=data$Rhizo
y1=data$Oryz
tro=data$t..days.

#We first define our competitive equations.
M = function(tro,y,parms){ 
  with(as.list(c(y,parms)),{   
    list(c(
      0.0574*X*(1-X/331)-ar*X*Y,
      0.0646*Y*(1-Y/445)-ao*X*Y))
  })
}

#Then, modfit is used to best estimate the competitive coefficients
#for each equation.
CompModel <- function(params0){ 
  with(as.list(c(params0)),{
    parms=c(ar=ar0,ao=ao0)
    yini=c(X=X0,Y=Y0);
    out = ode(yini,times=tro,func=M,parms)
    return((out[,2]-x1)^2+(out[,3]-y1)^2 )     
  })
}
params0=c(ar0=.000001,ao0=.000001, X0=2,Y0=2);
fit = modFit(CompModel,params0,lower=c(0),upper=c(5))
fit$par

#Obtain model prediction
times = seq(0,300,length=1000)
ar0 = fit$par[[1]]; ao0=fit$par[[2]];  X_fit=fit$par[[3]]; Y_fit=fit$par[[4]];
parms=c(ar=ar0,ao=ao0); yini=c(X=X_fit,Y=Y_fit);
Model = ode(yini,times=times,func=M,parms)

#Plot the competitive data.
plot(tro,x1,ylim=c(0,500),xlim=c(0,300),col="red",main="Inter-specific Competition Data",xlab="time(days)",ylab="population")
points(tro,y1,col="blue")

#Plot the competitive model.
lines(Model[,1],Model[,2],col="red",xlim=c(0,300),ylim=c(0,500),type='l',main="Inter-specific Competition Model",xlab="time(days)",ylab="Population")
lines(Model[,1],Model[,3],col="blue")
legend("topleft",pch=1,c("Rhizopertha Data","Oryzaephilus Data"),col=c("red","blue"),bty="n")
legend("left",lty=1,c("Rhiz Model", "Oryz Model"),col=c("red","blue"),bty="n")
#Then we will define the competition coefficients.
ar=fit$par[1]
ao=fit$par[2]

#Now we can perform stability analysis using null-clines.
#First, we will start by solving for the equillibrium value of
#interest.
f=function (X,Y) Rr*X*(1-X/Kr)-ar*X*Y
g=function (X,Y) Ro*Y*(1-Y/Ko)-ao*X*Y

N = function(t,y,parms){ 
  with(as.list(c(y,parms)),{   
    list(c(
      f(X,Y),
      g(X,Y)))
  })
}

#The null-clines are then plotted.
xn=seq(0,600,length=1000)
x_null = Rr*(Kr-xn)/(ar*Kr);
y_null = (Ko*Ro-ao*Ko*xn)/Ro;
plot(xn,x_null,main="Null-cline Plot for Logistic Inter-specific Competition",col="red",type="l",xlab="Rhizopertha Population",ylab="Oryzaephilus population",ylim = c(-200,800),xlim = c(0,500))  #Plot V-nullicline
lines(xn,y_null,col="blue")
abline(h=0,col="blue")
abline(v=0,col="red")
legend("topright",c("Rhizo","Oryz"),pch=19,col=c("red","blue"))

trajectory = function(x0,y0)   # Write a function that takes I.C. as 
{                               # input and plots trajectory with arrows
  times = seq(0,800,length=1000)
  init = c(X = x0, Y = y0)
  out = ode(init,times,N,NULL)
  lines (out[,2],out[,3],type="l")
  arrows(out[10,2],out[10,3],out[11,2],out[11,3],
         length=0.1,lwd=2)
}

plot(xn,x_null,main="Null-cline Plot for Logistic Inter-specific Competition",col="red",type="l",xlab="Rhizopertha Population",ylab="Oryzaephilus population",ylim = c(378,500),xlim = c(180,286))  #Plot V-nullicline
lines(xn,y_null,col="blue")
legend("topright",c("Rhizo","Oryz"),pch=19,col=c("red","blue"))
trajectory (180,380)
trajectory (180,500)
trajectory (200,500)
trajectory (215,380)
trajectory (220,500)
trajectory (240,380)
trajectory (240,500)
trajectory (275,380)
trajectory (280,500)


plot(xn,x_null,main="Null-cline Plot for Logistic Inter-specific Competition",col="red",type="l",xlab="Rhizopertha Population",ylab="Oryzaephilus population",ylim = c(-200,200),xlim = c(281,381))  #Plot V-nullicline
abline(h=0,col="blue")
legend("topright",c("Rhizo","Oryz"),pch=19,col=c("red","blue"))
trajectory (280,10)
trajectory (380,10)
trajectory (280,-10)
trajectory (380,-10)


#To evaluate the nature of the equillibrium, we can also solve
#for eigenvalues.
EquilFn = function (z) c(f(z[1],z[2]),g(z[1],z[2])) # Define set of equations to find roots
equil_out=fsolve(EquilFn,c(200,300)) # Find roots of equations (i.e., equilibrium points), include initial guess
equil=c(X=equil_out$x[1],Y=equil_out$x[2])

EquilFn2 = function (z) c(f(z[1],z[2]),g(z[1],z[2]))
equil_out2=fsolve(EquilFn2,c(330,0)) # Find roots of equations (i.e., equilibrium points), include initial guess
equil2=c(X=equil_out2$x[1],Y=equil_out2$x[2])

Jacob = jacobian.full(equil,N)  #Calculate Jacobian at equilibrium point
ei = eigen(Jacob) #Calculate eigenvalues and eigenvectors
ei

Jacob2 = jacobian.full(equil2,N)  #Calculate Jacobian at equilibrium point
ei2 = eigen(Jacob2) #Calculate eigenvalues and eigenvectors
ei2

equil3=c(X=0,Y=0)
Jacob3 = jacobian.full(equil3,N)  #Calculate Jacobian at equilibrium point
ei3 = eigen(Jacob3) #Calculate eigenvalues and eigenvectors
ei3
