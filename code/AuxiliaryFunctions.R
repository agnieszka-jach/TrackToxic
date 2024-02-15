#auxiliary functions

sim1mindata<-function(days)
#purpose: to simulate 1-min price&volume data for several consecutive days
#input:
#days [int] the number of days, >=2
#output: data.frame with 4 columns, Day, Time, Price, Volume
#in chr, chr, num, num
{
  t0<-as.POSIXct(paste("2050-01-01","09:30:00"),format="%Y-%m-%d %H:%M:%S",tz="EST")
  T<-as.POSIXct(paste("2050-01-01","15:59:00"),format="%Y-%m-%d %H:%M:%S",tz="EST")
  time<-seq(t0,T,"1 min")
  time<-format(time,"%H:%M:%S")#get the time only

  D<-Sys.Date()#eg, "2024-02-01"
  dates<-seq(from=D-days+1,to=D,by="1 day")
  
  vtrend<-0.05*((1:length(time))-length(time)/2)^2+20000
  #U-shaped trend for the volume
  
  df<-NULL
  for(i in 1:days)
  {
    #price follows an AR(1) process, shifted
    price<-100+arima.sim(list(order = c(1,0,0), ar = 0.9), n = length(time))#in USD
    volume<-vtrend+arima.sim(list(order = c(1,0,0), ar = 0.5), n = length(time),sd=1000)
    dftemp<-data.frame(Date=as.character(dates[i]),Time=as.character(time),
                       Price=round(as.vector(price),2),Volume=round(as.vector(volume)))
    df<-rbind(df,dftemp)
  }  
  
  return(df)  
}  

BVrule<-function(pricediff,sigmadiff=0)
#purpose: to calculate BulkVolume proportion
#input:
#pricediff [num] vector of (dollar) price differences \Delta P 
#for intraday data  
#if sigmadiff is not given (sigmadiff=0) then it is computed below  
#sigmadiff [num] scalar with the estimate of \sigma_{\Delta P}
#output: [num] vector of the same lenght as pricediff
#with the normal CDF of the normalized price-differences
{

  if(sigmadiff==0 || is.na(sigmadiff))
  {  
    s<-sd(pricediff,na.rm=TRUE)
  } else s<-sigmadiff
  
  z<-pricediff/s#vector of length n
  prop<-pnorm(z)
  return(prop)
}  

mycdf<-function(x)
#ECDF updated in realtime, allowing for NAs
{
  y<-x#store a copy of x in y
  if(sum(!is.na(x))>=1)#if at least one non-NA
  {
  x<-x[!is.na(x)]
  
  n<-length(x)
  Fy<-rep(NA,length(y))#cdf for the original vector
  Fx<-rep(NA,n)#cdf for the vector without the NA's

  #loop over non-NA values
  for(i in 1:n)
  {
    xx<-x[1:i]
      Fn<-ecdf(xx)
      Fx[i]<-Fn(xx)[i]  #output vector with the element 'i' updated
  }
  Fy[!is.na(y)]<-Fx#fillout the slots with non-NA's by the updated ecdf
  return(Fy)
  } else return(x)#all NAs
  
}  

#this function should be fixed so that it works for incomplete days and little trading
mintoind<-function(time=as.POSIXct(c("2030-01-01 17:30:00",
                                     "2030-01-01 21:30:00",
                                     "2030-01-01 22:59:00",
                                     "2030-01-01 21:32:00"),
                                   format="%Y-%m-%d %H:%M:%S",tz="EST"))
#purpose: to convert calendar time (in 1-min resolution, POSIXct, trading hours)
#to an integer between 1 and 390
#input: time in POSIXct 09:30-16:00h
#output: [int] corresponding integer tags  
{  
  tt<-format(time,"%H:%M:%S")#get the time only
  
  t0<-as.POSIXct(paste("2050-01-01","09:30:00"),format="%Y-%m-%d %H:%M:%S",tz="EST")
  T<-as.POSIXct(paste("2050-01-01","15:59:00"),format="%Y-%m-%d %H:%M:%S",tz="EST")

  myminutes<-seq(t0,T,"1 min")
  myindices<-1:390

  mytime<-as.POSIXct(paste("2050-01-01",tt),format="%Y-%m-%d %H:%M:%S",tz="EST")
  
  out<-myindices[myminutes %in% mytime]
  
  return(out)
}

myscale<-function(x)
#purpose: to normalize a vector to [0,1] in real-time  
{
  y<-x#save a copy of x in y
  x<-x[!is.na(x)]
  
  mi<-cummin(x)
  ma<-cummax(x)
  
  xs<-(x-mi)/(ma-mi)
  
  y[!is.na(y)]<-xs
  return(y)#return normalize y with the NAs
}  

