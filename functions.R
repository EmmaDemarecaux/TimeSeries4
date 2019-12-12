kalman <- function(Y,A,C,Sigma_1,Sigma_2,V0,Xhat0,n_pred)
{
  ## DEFINITION OF THE INPUT:
  
  ## Y      : Matrix with one observation per row. May be a vector if only one value
  ##          is observed at each time point.
  ## A, C, Sigma_1, Sigma_2 : Matrices as described in the definition of the Kalman
  ##          filter. 
  ## V0     : Initial covariance of Sigma_xx_{1|0}
  ## Xhat0  : Initial estimate of the state X_{1|0}
  ## n_pred: For making forecasts beyond the data
  
  ## DIMENSION OF THE SYSTEM
  n_obs <- dim(Y) #Number of observed parameters
  n_state <- dim(A)[1] #Number of state variables
  
  ## INITIALISATION
  X_hat <- Xhat0
  Sigma_xx <- V0
  Sigma_yy <- C%*%Sigma_xx%*%t(C)+Sigma_2
  
  ## PREPARATION OF THE ARRAYS USED TO SAVE THE RESULTS
  ## For saving the reconstructed and predicted state vectors
  
  X_rec <- array(dim=c(n_obs[1]+n_pred-1,n_state))
  X_pred <- array(dim=c(n_obs[1]+n_pred,n_state))
  
  ## For saving K, Sigmas
  K.out <- array(dim=c(dim(Sigma_xx%*%t(C)%*%solve(Sigma_yy)),n_obs[1]))
  Sigma_xx_rec <- array(dim=c(dim(Sigma_xx),n_obs[1]+n_pred-1))
  Sigma_yy_rec <- array(dim=c(dim(Sigma_yy),n_obs[1]+n_pred-1))
  Sigma_xx_pred <- array(dim=c(dim(Sigma_xx),n_obs[1]+n_pred))
  Sigma_yy_pred <- array(dim=c(dim(Sigma_yy),n_obs[1]+n_pred))
  ## RECONSTRUCTION AND PREDICTION OVER THE DATA SET
  for(tt in 1:n_obs[1])
  {
    K <- Sigma_xx%*%t(C)%*%solve(Sigma_yy)
    ## Reconstruction
    X_hat <- X_hat+K%*%(Y[tt,]-C %*% as.matrix(X_hat))
    X_rec[tt,] <- X_hat
    #Sigma_xx <- Sigma_xx-K%*%C%*%Sigma_xx
    Sigma_xx <- Sigma_xx-K%*%Sigma_yy%*%t(K)
    ## Saving of Reconstruction error variance-covariances matrices
    Sigma_xx_rec [,,tt] <- Sigma_xx
    Sigma_yy_rec [,,tt] <- Sigma_yy
    ## Prediction
    X_hat <- A%*%X_hat
    X_pred [tt+1,] <- X_hat
    Sigma_xx <- A%*%Sigma_xx%*%t(A)+Sigma_1
    Sigma_yy <- C%*%Sigma_xx%*%t(C)+Sigma_2
    ###Saving of Prediction error variance-covariances matrices
    K.out[,,tt] <- K
    Sigma_xx_pred [,,tt+1] <- Sigma_xx
    Sigma_yy_pred [,,tt+1] <- Sigma_yy
  }
  
  ## RECONSTRUCTION AND PREDICTION OVER THE PREDICTION SET
  for(tt in (n_obs[1]+1):(n_obs[1]+n_pred-1))
  {
    ## Reconstruction
    X_hat <- X_hat
    X_rec[tt,] <- X_hat
    #Sigma_xx <- Sigma_xx-K%*%C%*%Sigma_xx
    Sigma_xx <- Sigma_xx
    ## Saving of Reconstruction error variance-covariances matrices
    Sigma_xx_rec [,,tt] <- Sigma_xx
    Sigma_yy_rec [,,tt] <- Sigma_yy
    ## Prediction
    X_hat <- A%*%X_hat
    X_pred [tt+1,] <- X_hat
    Sigma_xx <- A%*%Sigma_xx%*%t(A)+Sigma_1
    Sigma_yy <- C%*%Sigma_xx%*%t(C)+Sigma_2
    ###Saving of Prediction error variance-covariances matrices
    
    Sigma_xx_pred [,,tt+1] <- Sigma_xx
    Sigma_yy_pred [,,tt+1] <- Sigma_yy
  }
  
  ## MAKING THE LIST OF THE OUTPUTS
  outputs <- list(rec=X_rec,pred=X_pred ,K=K.out,Sigma_xx_rec =Sigma_xx_rec ,
                  Sigma_yy_rec =Sigma_yy_rec ,Sigma_xx_pred =Sigma_xx_pred ,
                  Sigma_yy_pred =Sigma_yy_pred )
  return(outputs)
}

Reconstructions <- function(kal,x){
  
  #l'attribut rec de kal est X_rec
  Rec<-kal$rec
  Rec<-data.frame('x'=Rec[,1],'y'=Rec[,2])
  Rec_Var<-kal$Sigma_xx_rec
  Pred<-kal$pred
  Pred<-data.frame('x'=Pred[,1],'y'=Pred[,2])
  Pred_Var<-kal$Sigma_xx_pred
  
  if (x==0){
    
    #----------------------sh
    sd<-sqrt(Rec_Var[1,1,])
    t95_minus<-qnorm(0.025)*sd
    t95_plus<-qnorm(0.975)*sd
    
    Rec_plus<-Rec$x+t95_plus
    Rec_minus<-Rec$x+t95_minus
    
    sd_pred<-sqrt(Pred_Var[1,1,])
    t95_minus_pred<-qnorm(0.025)*sd_pred
    t95_plus_pred<-qnorm(0.975)*sd_pred
    
    Pred_plus<-Pred$x+t95_plus_pred
    Pred_minus<-Pred$x+t95_minus_pred
    
    plot(Year[1:169],Data$sh,type='l', col='orange',lwd=1.5, cex=0.85
         ,xlab='Time [Year]', cex.axis = 1,xlim=c(1850,2050),ylab='',
         main='Temperature anomalies in the Southern hemisphere',
         ylim=c(min(Rec_minus[1:169]),max(Rec_plus[1:169])))
    lines(Year[1:169],Rec_plus[1:169], lty=5, col= 'red', lwd=1.5)
    lines(Year[1:169],Rec_minus[1:169], lty=5, col= 'red', lwd=1.5)
    lines(Year[1:169],Rec$x[1:169],type = 'o', pch=1, lty=1, bty='l', col=blues[2])
    lines(Year[170:201],Pred$x[170:201], lty=1, col= 'purple', pch=1, lwd=1.5)
    lines(Year[170:201],Pred_plus[170:201], lty=1, col= 'red', lwd=1.5)
    lines(Year[170:201],Pred_minus[170:201], lty=1, col= 'red', lwd=1.5)
    legend('topleft',legend = c('Reconstruction','95% confidence interval reconstruction'
                                ,'Southern hemisphere','Predictions','95% confidence interval predictions'),
           col = c(blues[2],'red','orange','purple','red'),lty = c(1,5,1,1,1), pch = c(1,NA,NA,NA,NA),cex=1,
           bty = 'n')}
  
  if (x==1){
    
    #----------------------sh
    sd<-sqrt(Rec_Var[2,2,])
    t95_minus<-qnorm(0.025)*sd
    t95_plus<-qnorm(0.975)*sd
    
    Rec_plus<-Rec$y+t95_plus
    Rec_minus<-Rec$y+t95_minus
    
    sd_pred<-sqrt(Pred_Var[2,2,])
    t95_minus_pred<-qnorm(0.025)*sd_pred
    t95_plus_pred<-qnorm(0.975)*sd_pred
    
    Pred_plus<-Pred$y+t95_plus_pred
    Pred_minus<-Pred$y+t95_minus_pred
    
    plot(Year[1:169],Data$nh,type='l', col='orange',lwd=1.5, cex=0.85
         ,xlab='Time [Year]', cex.axis = 1,xlim=c(1850,2050),ylab='',
         main='Temperature anomalies in the Northern hemisphere',
         ylim=c(min(Rec_minus[1:169]),max(Rec_plus[1:169])))
    lines(Year[1:169],Rec_plus[1:169], lty=5, col= 'red', lwd=1.5)
    lines(Year[1:169],Rec_minus[1:169], lty=5, col= 'red', lwd=1.5)
    lines(Year[1:169],Rec$y[1:169], type = 'o', pch=1, lty=1, bty='l', col=blues[2])
    lines(Year[170:201],Pred$y[170:201], lty=1, col= 'purple', pch=1, lwd=1.5)
    lines(Year[170:201],Pred_plus[170:201], lty=1, col= 'red', lwd=1.5)
    lines(Year[170:201],Pred_minus[170:201], lty=1, col= 'red', lwd=1.5)
    legend('topleft',legend = c('Reconstruction','95% confidence interval reconstruction'
                                ,'Northern hemisphere','Predictions','95% confidence interval predictions'),
           col = c(blues[2],'red','orange','purple','red'),lty = c(1,5,1,1,1), pch = c(1,NA,NA,NA,NA),cex=1,
           bty = 'n')}
}

Loglike <- function(kal){
  Y_tilde <-Y[-1,]-kal$pred[2:n,]
  R = kal$Sigma_yy_pred[,,2:n] 
  
  P1 = matrix(NA, nrow = 168, ncol = 1)
  P2 = matrix(NA, nrow = 168, ncol = 1)
  
  for (i in 1:168){
    P1[i] = log(det(R[,,i]))
    P2[i] = t(Y_tilde[i,]) %*% solve(R[,,i]) %*% Y_tilde[i,]
  }
  return(-0.5 * sum(P1 + P2))}

Loglike2 <- function(kal){
  Y_tilde <-Y[-1,]-kal$pred[2:n,1:2]
  R = kal$Sigma_yy_pred[,,2:n] 
  
  P1 = matrix(NA, nrow = 168, ncol = 1)
  P2 = matrix(NA, nrow = 168, ncol = 1)
  
  for (i in 1:168){
    P1[i] = log(det(R[,,i]))
    P2[i] = t(Y_tilde[i,]) %*% solve(R[,,i]) %*% Y_tilde[i,]
  }
  return(-0.5 * sum(P1 + P2))}

Predictions <- function(Forecast,x,kal){
  Rec<-kal$rec
  Rec<-data.frame('x'=Rec[,1],'y'=Rec[,2])
  Rec_Var<-kal$Sigma_xx_rec
  Pred<-kal$pred
  Pred<-data.frame('x'=Pred[,1],'y'=Pred[,2])
  Pred_Var<-kal$Sigma_xx_pred
  
  if (x==0){
    sdM<-sqrt(Forecast$pred.var[1,1,])
    t95_minusM<-qnorm(0.025)*sdM
    t95_plusM<-qnorm(0.975)*sdM
    
    Pred_plusM<-x_f+t95_plusM
    Pred_minusM<-x_f+t95_minusM
    
    sd_pred<-sqrt(Pred_Var[1,1,])
    t95_minus_pred<-qnorm(0.025)*sd_pred
    t95_plus_pred<-qnorm(0.975)*sd_pred
    
    Pred_plus<-Pred$x+t95_plus_pred
    Pred_minus<-Pred$x+t95_minus_pred
    
    plot(Year[1:169],Data$sh,xlab='Year', ylab='', main="Temperature anomalies - Southern
         hemisphere", type = 'o',pch=1,lty=5,bty='l',cex=0.5,col=blues[9], cex.axis =1,
         xlim=c(1850,2050),ylim=c(min(Pred_minus[-c(1)]),max(Pred_plus[-c(1)])))
    lines(Year[170:201],Pred$x[170:201], lty=1, col= 'purple', pch=1, lwd=1.5)
    lines(Year[170:201],Pred_plus[170:201], lty=5, col= 'purple', lwd=1.5)
    lines(Year[170:201],Pred_minus[170:201], lty=5, col= 'purple', lwd=1.5)
    lines(Year[170:201],Pred_plusM, lty=5, col= 'red', lwd=1.5)
    lines(Year[170:201],Pred_minusM, lty=5, col= 'red', lwd=1.5)
    lines(Year[170:201],x_f,col='red', lwd = 1.5)
    
    
    legend('topleft',legend = c('Southern hemisphere','Predictions MARIMA', '95% confidence interval MARIMA'
                                ,'Predictions Kalman', '95% confidence interval Kalman'),
           col = c(blues[9],'red','red','purple','purple'),lty = c(1,1,5,1,5), pch = c(1,NA,NA,NA,NA),cex=0.7,
           bty = 'n')}    
  if (x==1){
    
    sdM<-sqrt(Forecast$pred.var[2,2,])
    t95_minusM<-qnorm(0.025)*sdM
    t95_plusM<-qnorm(0.975)*sdM
    
    Pred_plusM<-y_f+t95_plusM
    Pred_minusM<-y_f+t95_minusM
    
    sd_pred<-sqrt(Pred_Var[2,2,])
    t95_minus_pred<-qnorm(0.025)*sd_pred
    t95_plus_pred<-qnorm(0.975)*sd_pred
    
    Pred_plus<-Pred$y+t95_plus_pred
    Pred_minus<-Pred$y+t95_minus_pred
    
    plot(Year[1:169],Data$nh,xlab='Year', ylab='', main="Temperature anomalies - Northern
         hemisphere", type = 'o',pch=1,lty=5,bty='l',cex=0.5,col=blues[2], cex.axis =1,
         xlim=c(1850,2050),ylim=c(min(Pred_minus[-c(1)]),max(Pred_plus[-c(1)])))
    lines(Year[170:201],Pred$y[170:201], lty=1, col= 'purple', pch=1, lwd=1.5)
    lines(Year[170:201],Pred_plus[170:201], lty=5, col= 'purple', lwd=1.5)
    lines(Year[170:201],Pred_minus[170:201], lty=5, col= 'purple', lwd=1.5)
    lines(Year[170:201],Pred_plusM, lty=5, col= 'red', lwd=1.5)
    lines(Year[170:201],Pred_minusM, lty=5, col= 'red', lwd=1.5)
    lines(Year[170:201],y_f,col='red', lwd = 1.5)
    legend('topleft',legend = c('Norhtern hemisphere','Predictions MARIMA', '95% confidence interval MARIMA'
                                ,'Predictions Kalman', '95% confidence interval Kalman'),
           col = c(blues[2],'red','red','purple','purple'),lty = c(1,1,5,1,5), pch = c(1,NA,NA,NA,NA),cex=0.7,
           bty = 'n')} 
}

my.obj.4 <- function(PAR){
  Kro <- kalman(Y=Y, A=A, C=C, Sigma_1=diag(exp(PAR[4]),2), Sigma_2=diag(exp(PAR[5]),2), 
                Xhat0=c(PAR[1:2]), V0=diag(exp(PAR[3]),2),n_pred=n_pred)
  Y_tilde <-Y[-1,]-Kro$pred[2:n,]
  R = Kro$Sigma_yy_pred[,,2:n] 
  
  P1 = matrix(NA, nrow = 168, ncol = 1)
  P2 = matrix(NA, nrow = 168, ncol = 1)
  
  for (i in 1:168){
    P1[i] = log(det(R[,,i]))
    P2[i] = t(Y_tilde[i,]) %*% solve(R[,,i]) %*% Y_tilde[i,]
  }
  return(0.5 * sum(P1 + P2))}

my.obj.5 <- function(PAR){
  Kro <- kalman(Y=Y, A=A, C=C, 
                Sigma_1=matrix(c(exp(PAR[4])*exp(PAR[4]),exp(PAR[4])*exp(PAR[4])*PAR[6],
                                 exp(PAR[4])*exp(PAR[4])*PAR[6],exp(PAR[4])*exp(PAR[4])),nrow=2),
                Sigma_2=diag(exp(PAR[5]),2),Xhat0=c(PAR[1:2]), V0=diag(exp(PAR[3]),2),n_pred=n_pred)
  Y_tilde <-Y[-1,]-Kro$pred[2:n,]
  R = Kro$Sigma_yy_pred[,,2:n] 
  
  P1 = matrix(NA, nrow = 168, ncol = 1)
  P2 = matrix(NA, nrow = 168, ncol = 1)
  
  for (i in 1:168){
    P1[i] = log(abs(det(R[,,i])))
    P2[i] = t(Y_tilde[i,]) %*% solve(R[,,i]) %*% Y_tilde[i,]
  }
  return(0.5 * sum(P1 + P2))}

my.obj.7 <- function(PAR){
  Kro <- kalman(Y=Y, A=A, C=C, 
                Sigma_1=matrix(c(exp(PAR[5])*exp(PAR[5]),exp(PAR[5])*exp(PAR[5])*0.92,0,
                                 exp(PAR[5])*exp(PAR[5])*0.92,exp(PAR[5])*exp(PAR[5]),0,
                                 0,0,exp(PAR[6])),nrow=3),
                
                Sigma_2=diag(exp(PAR[7]),2),Xhat0=c(PAR[1:3]), 
                V0=diag(c(exp(PAR[4]),exp(PAR[4]),exp(PAR[8]))),n_pred=n_pred)
  
  Y_tilde <-Y[-1,]-Kro$pred[2:n,1:2]
  R = Kro$Sigma_yy_pred[,,2:n] 
  
  P1 = matrix(NA, nrow = 168, ncol = 1)
  P2 = matrix(NA, nrow = 168, ncol = 1)
  
  for (i in 1:168){
    P1[i] = log(det(R[,,i]))
    P2[i] = t(Y_tilde[i,]) %*% solve(R[,,i]) %*% Y_tilde[i,]
  }
  return(0.5 * sum(P1 + P2))}