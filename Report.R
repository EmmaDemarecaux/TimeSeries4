## This empties the work space ####
rm(list=ls())

## change directory

setwd("~/Travail/DTU/Time_series/Ass4")

#setwd("~/DTU courses/Time series/A4")
source("functions.R")

library(RColorBrewer)
library(xtable)
library(marima)
library(matlib)

Data <- read.table(file="A4_annual.txt", sep="\t", header=TRUE)
n<-length(Data$year)
n_pred<-2050-Data$year[n]

Time <- data.frame('t'=c(1:n))
Data <- cbind(Time,Data)
Year <- 1850:2050

#---------------------------------------------------------- Question 1

blues <- rev(brewer.pal(11, "RdBu"))

par(mfrow=c(2,1))
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))

plot(Data$year, Data$nh,xlab='Year', ylab='',
     main="Temperature anomalies - Northern hemisphere",
     type = 'o',pch=1,lty=5,bty='l',cex=0.5,
     col=blues[2], cex.axis = 0.7)

plot(Data$year, Data$sh,xlab='Year', ylab='',
     main="Temperature anomalies - Southern hemisphere",
     type = 'o',pch=1,lty=5,bty='l',cex=0.5,
     col=blues[9], cex.axis = 0.7)

par(mfrow=c(1,1))
par(mgp=c(2, 0.7,0), mar=c(3,3,2,1))
plot(Data$year,Data$nh,type = 'o',pch=1,lty=1,bty='l',cex=0.7,
     col=blues[3], cex.axis = 1,main='Temparature anomalies',ylab='',xlab='Year')
lines(Data$year,Data$sh,type = 'o',cex=0.8,col=blues[9])

legend('topleft',legend = c('Northern hemisphere','Southern hemisphere'),
  col = c(blues[3],blues[9]),lty = c(5,5), pch = c(1,1),cex=0.7,
  bty = n)

#------------------------------------------------------------------- Question 3
#Manual kalmann filter: 
## Simulating the process:
TS0 <- -0.4
TN0 <- -0.3
A=matrix(nrow = 2, ncol = 2,byrow=T,
           c(1,0,
             0,1))
C=matrix(nrow = 2, ncol = 2, byrow = T,
          c(1,0,
            0,1))
Sigma_1 <- matrix(c(0.01,0,0,0.01),nrow=2)
Sigma_2 <- matrix(c(0.01,0,0,0.01),nrow=2)
V0=Sigma_1

X <- matrix(nrow=2,ncol=n)
X[,1] <- c(TS0,TN0)
Xhat0<- X[,1]

Y<-cbind(Data$sh,Data$nh)

kal <- kalman(Y=Y, A=A, C=C, Sigma_1=Sigma_1, Sigma_2=Sigma_2, Xhat0=Xhat0, V0=V0, n_pred=n_pred)

#----------------------sh
Reconstructions(kal,0)

#---------------------- nh
Reconstructions(kal,1)

Loglike(kal)
#sh
Pred<-kal$pred
Pred<-data.frame('x'=Pred[,1],'y'=Pred[,2])
Pred_Var<-kal$Sigma_xx_pred

sd_pred<-sqrt(Pred_Var[1,1,])
t95_minus_pred<-qnorm(0.025)*sd_pred
t95_plus_pred<-qnorm(0.975)*sd_pred

Pred_plus<-Pred$x+t95_plus_pred
Pred_minus<-Pred$x+t95_minus_pred

table <- matrix(cbind(Pred$x[c(171,181,191,201)],Pred_minus[c(171,181,191,201)],
                      Pred_plus[c(171,181,191,201)]),nrow=3,byrow=TRUE)
rownames(table) <- c("Predictions","Lower bounderies","Upper bounderies")
colnames(table) <- c("2020","2030","2040","2050")
table <- as.table(table)
table

#nh
sd_pred<-sqrt(Pred_Var[2,2,])
t95_minus_pred<-qnorm(0.025)*sd_pred
t95_plus_pred<-qnorm(0.975)*sd_pred

Pred_plus<-Pred$y+t95_plus_pred
Pred_minus<-Pred$y+t95_minus_pred

table <- matrix(cbind(Pred$y[c(171,181,191,201)],Pred_minus[c(171,181,191,201)],
                      Pred_plus[c(171,181,191,201)]),nrow=3,byrow=TRUE)
rownames(table) <- c("Predictions","Lower bounderies","Upper bounderies")
colnames(table) <- c("2020","2030","2040","2050")
table <- as.table(table)
table

#----------------------------------------------------------- Question 4 done !!!!
xpar = -0.4
ypar = -0.3

stdpar = log(0.01)
Sigma1_par = log(0.01)
Sigma2_par = log(0.01)

PAR <- c(xpar,ypar,stdpar,Sigma1_par,Sigma2_par)

Kmopt <- optim(PAR,my.obj.4,method = "L-BFGS-B")

Xhat0 = Kmopt$par[1:2]
V0=diag(exp(Kmopt$par[3]),2)
Sigma_1 = diag(exp(Kmopt$par[4]),2)
Sigma_2 = diag(exp(Kmopt$par[5]),2)

kal <- kalman(Y=Y, A=A, C=C, Sigma_1=Sigma_1, Sigma_2=Sigma_2, 
              Xhat0=Xhat0, V0=V0, n_pred=n_pred)
Loglike(kal)
Reconstructions(kal,0)
Reconstructions(kal,1)

#sh
Pred<-kal$pred
Pred<-data.frame('x'=Pred[,1],'y'=Pred[,2])
Pred_Var<-kal$Sigma_xx_pred

sd_pred<-sqrt(Pred_Var[1,1,])
t95_minus_pred<-qnorm(0.025)*sd_pred
t95_plus_pred<-qnorm(0.975)*sd_pred

Pred_plus<-Pred$x+t95_plus_pred
Pred_minus<-Pred$x+t95_minus_pred

table <- matrix(cbind(Pred$x[c(171,181,191,201)],Pred_minus[c(171,181,191,201)],
                      Pred_plus[c(171,181,191,201)]),nrow=3,byrow=TRUE)
rownames(table) <- c("Predictions","Lower bounderies","Upper bounderies")
colnames(table) <- c("2020","2030","2040","2050")
table <- as.table(table)
table

#nh
sd_pred<-sqrt(Pred_Var[2,2,])
t95_minus_pred<-qnorm(0.025)*sd_pred
t95_plus_pred<-qnorm(0.975)*sd_pred

Pred_plus<-Pred$y+t95_plus_pred
Pred_minus<-Pred$y+t95_minus_pred

table <- matrix(cbind(Pred$y[c(171,181,191,201)],Pred_minus[c(171,181,191,201)],
                      Pred_plus[c(171,181,191,201)]),nrow=3,byrow=TRUE)
rownames(table) <- c("Predictions","Lower bounderies","Upper bounderies")
colnames(table) <- c("2020","2030","2040","2050")
table <- as.table(table)
table

#----------------------------------------------------------- Question 5 done !!!!
#This time we had a rho qui doit appartenir ? [-1,1]
xpar = -0.3
ypar = -0.4

stdpar = log(0.01)
Sigma1_par = log(0.01)
Sigma2_par = log(0.01)
rho = 0.2
PAR <- c(xpar,ypar,stdpar,Sigma1_par,Sigma2_par,rho)

Kmopt <- optim(PAR,my.obj.5,method = "L-BFGS-B")

Xhat0 = Kmopt$par[1:2]
V0=diag(exp(Kmopt$par[3]),2)
rho=Kmopt$par[6]
Sigma_1 = matrix(c(exp(Kmopt$par[4])*exp(Kmopt$par[4]),exp(Kmopt$par[4])*exp(Kmopt$par[4])*rho,
                   exp(Kmopt$par[4])*exp(Kmopt$par[4])*rho,exp(Kmopt$par[4])*exp(Kmopt$par[4])),nrow=2)
Sigma_2 = diag(exp(Kmopt$par[5]),2)


kal <- kalman(Y=Y, A=A, C=C, Sigma_1=Sigma_1, Sigma_2=Sigma_2, 
              Xhat0=Xhat0, V0=V0, n_pred=n_pred)
Loglike(kal)
#----------------------sh
Reconstructions(kal,0)
#---------------------- nh
Reconstructions(kal,1)

#sh
Pred<-kal$pred
Pred<-data.frame('x'=Pred[,1],'y'=Pred[,2])
Pred_Var<-kal$Sigma_xx_pred

sd_pred<-sqrt(Pred_Var[1,1,])
t95_minus_pred<-qnorm(0.025)*sd_pred
t95_plus_pred<-qnorm(0.975)*sd_pred

Pred_plus<-Pred$x+t95_plus_pred
Pred_minus<-Pred$x+t95_minus_pred

table <- matrix(cbind(Pred$x[c(171,181,191,201)],Pred_minus[c(171,181,191,201)],
                      Pred_plus[c(171,181,191,201)]),nrow=3,byrow=TRUE)
rownames(table) <- c("Predictions","Lower bounderies","Upper bounderies")
colnames(table) <- c("2020","2030","2040","2050")
table <- as.table(table)
table

#nh
sd_pred<-sqrt(Pred_Var[2,2,])
t95_minus_pred<-qnorm(0.025)*sd_pred
t95_plus_pred<-qnorm(0.975)*sd_pred

Pred_plus<-Pred$y+t95_plus_pred
Pred_minus<-Pred$y+t95_minus_pred

table <- matrix(cbind(Pred$y[c(171,181,191,201)],Pred_minus[c(171,181,191,201)],
                      Pred_plus[c(171,181,191,201)]),nrow=3,byrow=TRUE)
rownames(table) <- c("Predictions","Lower bounderies","Upper bounderies")
colnames(table) <- c("2020","2030","2040","2050")
table <- as.table(table)
table
#---------------------------------------------------Question 6

TS0 <- -0.3
TN0 <- -0.16
TU0 <- 0
#A <- matrix(c(1,0,0,1),nrow=2)
A=matrix(nrow = 3, ncol = 3,byrow=T,
          c(1,0,1,
            0,1,1,
            0,0,1))

#C <- matrix(c(1,0,0,1),nrow=2)
C=matrix(nrow = 2, ncol = 3, byrow = T,
          c(1,0,0,
            0,1,0))

Sigma_1 <- matrix(c(0.01,0,0,0,0.01,0,0,0,0.0001),nrow=3)
Sigma_2 <- matrix(c(0.01,0,0,0.01),nrow=2)

V0<-matrix(c(0.0001,0,0,0,0.0001,0,0,0,0.0001),nrow=3)

X <- matrix(nrow=3,ncol=n)
X[,1] <- c(TN0,TS0,TU0)
Xhat0<- X[,1]

Y<-cbind(Data$sh,Data$nh)

kal <- kalman(Y=Y, A=A, C=C, Sigma_1=Sigma_1, Sigma_2=Sigma_2, Xhat0=Xhat0, V0=V0, n_pred=n_pred)
Reconstructions(kal,0)
Reconstructions(kal,1)
Loglike2(kal)

#sh
Pred<-kal$pred
Pred<-data.frame('x'=Pred[,1],'y'=Pred[,2])
Pred_Var<-kal$Sigma_xx_pred

sd_pred<-sqrt(Pred_Var[1,1,])
t95_minus_pred<-qnorm(0.025)*sd_pred
t95_plus_pred<-qnorm(0.975)*sd_pred

Pred_plus<-Pred$x+t95_plus_pred
Pred_minus<-Pred$x+t95_minus_pred

table <- matrix(cbind(Pred$x[c(171,181,191,201)],Pred_minus[c(171,181,191,201)],
                      Pred_plus[c(171,181,191,201)]),nrow=3,byrow=TRUE)
rownames(table) <- c("Predictions","Lower bounderies","Upper bounderies")
colnames(table) <- c("2020","2030","2040","2050")
table <- as.table(table)
table

#nh
sd_pred<-sqrt(Pred_Var[2,2,])
t95_minus_pred<-qnorm(0.025)*sd_pred
t95_plus_pred<-qnorm(0.975)*sd_pred

Pred_plus<-Pred$y+t95_plus_pred
Pred_minus<-Pred$y+t95_minus_pred

table <- matrix(cbind(Pred$y[c(171,181,191,201)],Pred_minus[c(171,181,191,201)],
                      Pred_plus[c(171,181,191,201)]),nrow=3,byrow=TRUE)
rownames(table) <- c("Predictions","Lower bounderies","Upper bounderies")
colnames(table) <- c("2020","2030","2040","2050")
table <- as.table(table)
table
#--------------------------------------------------------Question 7
xpar = -0.3
ypar = -0.16
upar = 0.2

stdpar = log(0.001)
Sigma1_par = log(0.001)
Sigma1_paru = log(0.02)
Sigma2_par = log(0.001)
VO_paru = log(0.01)
rho=0.92

A=matrix(nrow = 3, ncol = 3,byrow=T,
         c(1,0,1,
           0,1,1,
           0,0,1))
C=matrix(nrow = 2, ncol = 3, byrow = T,
         c(1,0,0,
           0,1,0))
PAR <- c(xpar,ypar,upar,stdpar,Sigma1_par,Sigma1_paru,Sigma2_par,VO_paru)

Kmopt <- optim(PAR,my.obj.7,method = "L-BFGS-B")
  
Xhat0 = Kmopt$par[1:3]
V0=diag(c(exp(Kmopt$par[4]),exp(Kmopt$par[4]),exp(Kmopt$par[8])))
  
Sigma_1 = matrix(c(exp(Kmopt$par[5])*exp(Kmopt$par[5]),exp(Kmopt$par[5])*exp(Kmopt$par[5])*rho,0,
                   exp(Kmopt$par[5])*exp(Kmopt$par[5])*rho,exp(Kmopt$par[5])*exp(Kmopt$par[5]),0,
                  0,0,exp(Kmopt$par[6])),nrow=3)
Sigma_2 = diag(exp(Kmopt$par[7]),2)

kal <- kalman(Y=Y, A=A, C=C, Sigma_1=Sigma_1, Sigma_2=Sigma_2, 
              Xhat0=Xhat0, V0=V0, n_pred=n_pred)
Loglike2(kal)

#----------------------sh
Reconstructions(kal,0)
#---------------------- nh
Reconstructions(kal,1)

#sh
Pred<-kal$pred
Pred<-data.frame('x'=Pred[,1],'y'=Pred[,2])
Pred_Var<-kal$Sigma_xx_pred

sd_pred<-sqrt(Pred_Var[1,1,])
t95_minus_pred<-qnorm(0.025)*sd_pred
t95_plus_pred<-qnorm(0.975)*sd_pred

Pred_plus<-Pred$x+t95_plus_pred
Pred_minus<-Pred$x+t95_minus_pred

table <- matrix(cbind(Pred$x[c(171,181,191,201)],Pred_minus[c(171,181,191,201)],
                      Pred_plus[c(171,181,191,201)]),nrow=3,byrow=TRUE)
rownames(table) <- c("Predictions","Lower bounderies","Upper bounderies")
colnames(table) <- c("2020","2030","2040","2050")
table <- as.table(table)
table

#nh
sd_pred<-sqrt(Pred_Var[2,2,])
t95_minus_pred<-qnorm(0.025)*sd_pred
t95_plus_pred<-qnorm(0.975)*sd_pred

Pred_plus<-Pred$y+t95_plus_pred
Pred_minus<-Pred$y+t95_minus_pred

table <- matrix(cbind(Pred$y[c(171,181,191,201)],Pred_minus[c(171,181,191,201)],
                      Pred_plus[c(171,181,191,201)]),nrow=3,byrow=TRUE)
rownames(table) <- c("Predictions","Lower bounderies","Upper bounderies")
colnames(table) <- c("2020","2030","2040","2050")
table <- as.table(table)
table
#--------------------------------------------------------Question 8
#-------------------------MARIMA MODEL
#--------------------MARIMA(0,1,0) with 2 variables

#-PREDICTION over the reconstructed data
#check help(arma.forecast) and dif.poly argument : "(most often) output from the function define.dif holding the ar-representation of the differencing 
#polynomial (define.dif$dif.poly). If a differenced timeseries was analysed by marima the forecast-variance/covariance matrices are calculated for the 
#aggregated (original) timeseries if 'dif.poly' is specified. If not, the forecast-variance/covariance matrices are calculated for the differenced time 
#series. If forecasting is wanted for the original (not differenced) time series the 'dif.poly' created by define.dif must be specified."
#followed those slides : https://learn.inside.dtu.dk/d2l/le/content/13013/viewContent/3849/View 

Rec<-kal$rec[1:169,]
Rec<-data.frame('x'=Rec[,1],'y'=Rec[,2])
series=Rec
ar=c(1)
ma=c(0)
kvar=2 #number of variables
struct_1<-define.model(kvar=kvar,ar=ar,ma=ma)

difference=matrix(c(1,1,2,1),ncol=2) #variable 1 is differenced 1 and variable 2 as well 
dif.poly<- define.dif(series=series,difference=difference)
series_diff<-t(dif.poly$dif.series)

model<-marima(series_diff, ar.pattern=NULL, ma.pattern=NULL,  
              Plot="none",Check=FALSE, penalty=0)

data_to_predict<-data.frame('x'=matrix(NA, nrow = n_pred, ncol = 1),
                            'y'=matrix(NA, nrow = n_pred, ncol = 1))
#creating a dataframe with all the data (observed and to predict): necessary for the forecast function...
series_full_diff<-rbind(series_diff,data_to_predict)
series_full<-rbind(series,data_to_predict)
Forecast<-arma.forecast(series = series_full , marima = model, nstart=n,
                        nstep = n_pred, dif.poly=dif.poly$dif.poly)

#not sure the output is "differenciated"! En effet, la fonction arma.forecast dit pour dif.poly" If forecasting is wanted for the original 
#(not differenced) time series the 'dif.poly' created by define.dif must be specified."
#----------------------store the forecasts for the model we want and plot 
x_f<-Forecast$forecasts[1,(n+1):(n+n_pred)]
y_f<-Forecast$forecasts[2,(n+1):(n+n_pred)]
#plot the forecasts+CI
#x
Predictions(Forecast,0,kal)

Predictions(Forecast,1,kal)

#sh
sd<-sqrt(Forecast$pred.var[1,1,])
t95_minus<-qnorm(0.025)*sd
t95_plus<-qnorm(0.975)*sd

Pred_plus<-x_f+t95_plus
Pred_minus<-x_f+t95_minus

table <- matrix(cbind(Pred$x[c(171,181,191,201)],Pred_minus[c(171,181,191,201)],
                      Pred_plus[c(171,181,191,201)]),nrow=3,byrow=TRUE)
rownames(table) <- c("Predictions","Lower bounderies","Upper bounderies")
colnames(table) <- c("2020","2030","2040","2050")
table <- as.table(table)
table

#nh
sd_pred<-sqrt(Pred_Var[2,2,])
t95_minus_pred<-qnorm(0.025)*sd_pred
t95_plus_pred<-qnorm(0.975)*sd_pred

Pred_plus<-Pred$y+t95_plus_pred
Pred_minus<-Pred$y+t95_minus_pred

table <- matrix(cbind(Pred$y[c(171,181,191,201)],Pred_minus[c(171,181,191,201)],
                      Pred_plus[c(171,181,191,201)]),nrow=3,byrow=TRUE)
rownames(table) <- c("Predictions","Lower bounderies","Upper bounderies")
colnames(table) <- c("2020","2030","2040","2050")
table <- as.table(table)
table
