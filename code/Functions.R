rm(list=ls())

#A.Jach, A general comovement measure for times series,
#Mathematical and Statistical Methods for Actuarial Sciences and Finance, 
#M. Corazza, et.al., (Eds.), Springer, 279-284, 2021 

#P. Fryzlewicz, H.-S.Oh, Thick-pen transformation for time series,
#Journal of the Royal Statistical Society Series B, 73, 2011

#a set of functions to define TPT, TPMA, MTPT, MTTPMA (plus some auxiliary functions);
#they operate on 'lists'

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FORWARD-LOOKING!!!
#for TPT calculation
Boundary<-function(y,x,xout=NULL,tau=2,stat,...)#forward-looking
  #y: [num] y-axis vector of the same length as x
  #x: [num] x-axis vector of the same length as y (not necessarily equally-spaced)
  #xout: [num] optional (desired) time grid that is possibly different from x
  #tau: [num] scalar, pen thickness
  #stat: min, max, etc (no quotes)
  #output:
  #numerical vector with the boundary values 
  #corresponding to the x-values either a) as in x or b) as in xout (if xout is not null)
{
  if(is.null(xout))
  {  
    I<-length(x)
    result<-rep(NA,I)
    
    for(i in 1:I)
    {
      logical<-(x>=(x[i]) & (x<=x[i]+tau) )#original definition for now!!!
      #if pen is not thick enough, the result will be the same as y (and the boundaries L,U will be identical)
      ysub<-y[logical]
      result[i]<-stat(ysub,...)
    }
  } else{
    
    I<-length(xout)
    result<-rep(NA,I)
    
    for(i in 1:I)
    {
      logical<-(x>=(xout[i]) & (x<=xout[i]+tau) )#original definition for now!!!
      #here 'logical' can consist of FALSEs only, eg, when there are no time stamps of 
      #x that fall inside [xout[i],xout[i]+tau]
      ysub<-y[logical]
      result[i]<-stat(ysub,...)
    }
  }
  result[!is.finite(result)]<-NA
  return(result)
  
}  
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BACKWARD-LOOKING!!!
Boundary<-function(y,x,xout=NULL,tau=2,stat,...)#backward-looking
  #y: [num] y-axis vector of the same length as x
  #x: [num] x-axis vector of the same length as y (not necessarily equally-spaced)
  #xout: [num] optional (desired) time grid that is possibly different from x
  #tau: [num] scalar, pen thickness
  #stat: min, max, etc (no quotes)
  #output:
  #numerical vector with the boundary values 
  #corresponding to the x-values either a) as in x or b) as in xout (if xout is not null)
{
  if(is.null(xout))
  {  
    I<-length(x)
    result<-rep(NA,I)
    
    for(i in 1:I)
    {
      logical<-(x>=(x[i]-tau) & (x<=x[i]) )#modified definition
      #if pen is not thick enough, the result will be the same as y (and the boundaries L,U will be identical)
      ysub<-y[logical]
      result[i]<-stat(ysub,...)
    }
  } else{
    
    I<-length(xout)
    result<-rep(NA,I)
    
    for(i in 1:I)
    {
      logical<-(x>=(xout[i]-tau) & (x<=xout[i]) )#modified definition
      #here 'logical' can consist of FALSEs only, eg, when there are no time stamps of 
      #x that fall inside [xout[i],xout[i]+tau]
      ysub<-y[logical]
      result[i]<-stat(ysub,...)
    }
  }
  
  result[!is.finite(result)]<-NA
  return(result)
  
}  
TPT<-function(X,TIME,xout=NULL,tau)
  #purpose: calculate TPT, for a given tau, for K time series;
  #time series can have different lengths and each can be observed at different time stamps collected in TIME
  #input:
  #X: [list] K time series 
  #TIME: [list] corresponding list of (possibly) irregularly spaced time stamps (possibly, of different lengths)
  #xout: [num] optional time grid that is possibly different from TIME that is common to all K time series
  #tau: [num] scalar/number, pen thickness  
  #output:
  #tpt: [list] with L, U, TIME, xout, tau
  #L,U are lists with boundaries computed for TIME or xout
{
  K<-length(X)
  L<-U<-vector('list',length=K)
  
  #xout can be NULL because 'Boundary' accomodates that
  for(k in 1:K)
  {
    L[[k]]<-Boundary(y=X[[k]],x=TIME[[k]],xout=xout,tau=tau,stat=min,na.rm=TRUE)
    U[[k]]<-Boundary(y=X[[k]],x=TIME[[k]],xout=xout,tau=tau,stat=max,na.rm=TRUE)
  }
  
  mylist<-list(L=L,U=U,TIME=TIME,xout=xout,tau=tau)
  #if xout is null, then original time stamps of X^k [TIME] was used to compute TPT
  return(mylist)
}  

MTPT<-function(X,TIME,xout=NULL,tau)
  #purpose: calculate multi-thickness TPT (as in multi-thickness TPMA, Jach 2021),
  #for K time series for a vector of tau tickness values;
  #time series can have different lengths and each can be observed at different time stamps collected in TIME
  #input:
  #X: [list] K time series 
  #TIME: [list] corresponding list of (possibly) irregularly spaced time stamps (possibly, of different lengths)
  #xout: [num] optional time grid that is possibly different from TIME that is common to all K time series
  #tau: [num] vector with pen thicknesses of length K (NOT scalar)  
  #output:
  #tpt: [list] with L, U, TIME, xout, tau
  #L,U are lists with boundaries computed for TIME or xout
#each boundary is obtained with different pen thickness,
#ie, tau[k], k=1,2,...,K
{
  K<-length(X)
  L<-U<-vector('list',length=K)
  
  for(k in 1:K)
  {
    L[[k]]<-Boundary(y=X[[k]],x=TIME[[k]],xout=xout,tau=tau[k],stat=min,na.rm=TRUE)
    U[[k]]<-Boundary(y=X[[k]],x=TIME[[k]],xout=xout,tau=tau[k],stat=max,na.rm=TRUE)
  }
  
  mtpt<-list(L=L,U=U,TIME=TIME,xout=xout,tau=tau)
  #if xout is null, then original time stamps of X^k [TIME] was used to compute TPT
  return(mtpt)
}  

#TPMA and MTTPMA can be combined

TPMA<-function(X,TIME,xout=NULL,tau)
  #purpose: calculate TPMA, for a given tau, for K time series;
  #input:
  #X: [list] K time series 
  #TIME: [list] corresponding list of (possibly) irregularly spaced time stamps (possibly, of different lengths)
  #xout: [num] optional time grid that is possibly different from TIME that is common to all K time series
  #tau: [num] scalar/number, pen thickness  
  #output:
  #df: [df] data frame quantifying co-movement
  #on a time scale of (tau+1) of K time series over time given in xout or 'xout' computed as the combined grid
{
  #here xout cannot be NULL because we need to align U,L to measure comovement
  #if it's NULL, we combine all time grids and use that
  
  K<-length(X)
  if(is.null(xout))
  {
    temp<-NULL
    for(k in 1:K)
      temp<-c(temp,TIME[[k]])
    
    xout<-unique(sort(temp))#combined time grid
  }  
  
  temp<-TPT(X=X,TIME=TIME,xout=xout,tau=tau)
  xout<-temp$xout#adjustment here in case xout was NULL
  
  #now we CAN combine the slots of the list into a data frame because they have the same nr of elements
  L<-as.data.frame(temp$L)
  U<-as.data.frame(temp$U)
  
  MinU<-apply(U,MARGIN=1,min,na.rm=TRUE)
  MaxL<-apply(L,MARGIN=1,max,na.rm=TRUE)
  MaxU<-apply(U,MARGIN=1,max,na.rm=TRUE)
  MinL<-apply(L,MARGIN=1,min,na.rm=TRUE)
  
  #TPMA
  rho<-(MinU-MaxL)/(MaxU-MinL)
  
  #TPMA cannot be equal to -1, so replace such values by NAs
  rho[rho==-1]<-NA
  rho[!is.finite(rho)]<-NA
  df<-data.frame(rho=rho,xout=xout)
  return(df)
}  

MTTPMA<-function(X,TIME,xout=NULL,tau)
  #purpose: calculate MTTPMA (multi-thickness TPMA, Jach 2021),
  #for K time series for a vector of tau thickness values;
  #input:
  #X: [list] K time series 
  #TIME: [list] corresponding list of (possibly) irregularly spaced time stamps (possibly, of different lengths)
  #xout: [num] optional time grid that is possibly different from TIME that is common to all K time series
  #tau: [num] vector with pen thicknesses of length K (NOT scalar)  
  #output:
  #MTTPMA: [df] data frame quantifying cross-term co-movement
  #on a time scales of (tau[k]+1) of K time series over time given in xout (as given or if not, as calculated below)
{
  K<-length(X)
  
  #here xout cannot be NULL because we need to align U,L to measure comovement
  #if it's NULL, we combine all time grids and use that
  K<-length(X)
  if(is.null(xout))
  {
    temp<-NULL
    for(k in 1:K)
      temp<-c(temp,TIME[[k]])
    
    xout<-unique(sort(temp))#combined time grid
  }  
  
  temp<-MTPT(X=X,TIME=TIME,xout=xout,tau=tau)#we compute MTPT outside
  
  xout<-temp$xout#adjustment here in case xout was NULL
  
  #now we CAN combine the slots of the list into a data frame because they have the same nr of elements
  L<-as.data.frame(temp$L)
  U<-as.data.frame(temp$U)
  
  MinU<-apply(U,MARGIN=1,min,na.rm=TRUE)
  MaxL<-apply(L,MARGIN=1,max,na.rm=TRUE)
  MaxU<-apply(U,MARGIN=1,max,na.rm=TRUE)
  MinL<-apply(L,MARGIN=1,min,na.rm=TRUE)
  
  #MTTPMA
  rho<-(MinU-MaxL)/(MaxU-MinL)
  
  rho[rho==-1]<-NA
  rho[!is.finite(rho)]<-NA
  
  df<-data.frame(rho=rho,time=xout)
  return(df)
}  
