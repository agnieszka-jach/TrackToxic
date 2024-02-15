rm(list=ls())
graphics.off()


library('data.table')
library('ggplot2')

#set working directory to the location of the following
source("Functions.R")
source("AuxiliaryFunctions.R")

set.seed(7)

#specifications for data simulation
##########
STOCKS<-5
DAYS<-7

#specification for data analysis
#############
VUSERDENOM<-100#used to compute bucket size
TAU<-30#thickness of a pen
SMOOTH<-25 #vpin smoothing window size

#simulated HF data and place it in subfolder simdata

for(i in 1:STOCKS)
{  
  df<-sim1mindata(days=DAYS)
  filename<-paste0("../simdata/Asset",i,".txt")
  write.table(x=df,file=filename,row.names=FALSE)
}

setwd("../simdata/")

myfiles<-list.files(pattern=".txt")

L<-length(myfiles)

DT<-NULL

selected<-1:L#select all assets

for(i in 1:length(selected))
{

#load the data for one asset
l<-selected[i]  
myfile<-myfiles[l]
hf<-read.table(myfile,header=TRUE)
names(hf)<-c("Date","Time","Price","Volume")
#str(hf)
#head(hf)
#tail(hf)

hf$Time<-as.POSIXct(paste("2050-01-01",hf$Time),format="%Y-%m-%d %H:%M:%S",tz="EST")#EST ie New York
#add a 'fake' date here ("2050-01-01") - the interest is only in 'Time'

hf$DateTime<-as.POSIXct(paste(hf$Date,format(hf$Time,"%H:%M:%S")) ,format="%Y-%m-%d %H:%M:%S",tz="EST")
#combine Date and Time columns into a POSIXct column that indicates date and time

#str(hf)
hf$Volume<-round(hf$Volume/1000)#volume in 1000s

df<-data.frame(DateTime=(hf$DateTime),Time=(hf$Time),Date=(hf$Date),
               P=(hf$Price),V=(hf$Volume))#reversed the order if necessary


df$yday<-yday(df$Date)
df$ydayminus1<-df$yday-1

df<-data.table(df)#switch to data.table format to perform Step2
#this is better than excluding 09:30:00 cases, because 
#it could happen that the first transaction happens later (unlikely)

#day-wise computations
#add tagging with index i of 1-min interval (per day);
#compute (dollar) price changes, omitting the overnight return
#add normalized price as well as Vuser (asset-specific)
df<-df[,.(i=mintoind(Time),DateTime,Time,yday,ydayminus1,V,P,
          Pnorm=myscale(P),Vnorm=myscale(V),
          deltaP=c(NA,diff(P)),Vuser=round(sum(V)/VUSERDENOM)),by=.(Date)]


#previous-day info: Vuser
setkey(df,Date)

#create a temporary/auxiliary data.table with V=Vuser=bucket size (day-specific)
temp<-df[,.(Vuser=mean(Vuser,na.rm=TRUE)),by=.(Date)]#Vuser is constant within the same Date

temp[,`:=`(Vuserprev=shift(Vuser, n=1, fill=Vuser[1L], type="lag"))]#replace the first unavailable value by the first available

temp[,c("Vuser"):=NULL]#drop Vuser 


df<-merge(df,temp)#merge so that Vuserprev is concatenated with the 'large' data.table 
#it is important that both data.tables have 'Date' as the key

#Step 2: extend all the variables by repetition (volume-weighted)
#in rep, each element of 'x' is repeated according to each element of 'times'
dfex<-data.table(deltaP=rep(df$deltaP,times=df$V),
                 P=rep(df$P,times=df$V),
                 Pnorm=rep(df$Pnorm,times=df$V),
                 Vnorm=rep(df$Vnorm,times=df$V),
                 Vuser=rep(df$Vuser,times=df$V),
                 Vuserprev=rep(df$Vuserprev,times=df$V),
                 i=rep(df$i,times=df$V),DateTime=rep(df$DateTime,times=df$V),
                 Time=rep(df$Time,times=df$V),
                 Date=rep(df$Date,times=df$V),
                 yday=rep(df$yday,times=df$V),
                 ydayminus1=rep(df$ydayminus1,times=df$V))

#Step3: to add the running index j within a given day, use data.table format

dfex[,`:=`(j=1:length(i)),by=.(Date)]#j is the running index, done within a Date

#add tagging
dfex[,`:=`(TAU=ceiling(length(i)/Vuserprev[1])),by=.(Date)]#TAU is the total number of buckets of a given day

#TAU has identical elements per day, so we take the first value from the vector;
#then we have to truncate the tau vector (within each day) to have full buckets of a given size
dfex[,`:=`(tau=rep(1:TAU[1],each=Vuserprev[1])[1:length(j)]),by=.(Date)]#TAU is the total number of buckets per day

#the remaining steps are done using data.table, not in loops
#some variation in computing VB are done

#simplify the name
dt<-dfex

#previous-day calculations: 
#sd: date-wise (one value per day), volume weighted

setkey(dt,Date)
temp<-dt[,.(sdday=sd(deltaP,na.rm=TRUE)),by=.(Date)]#Date because we get one estimate per day
#add the previous available day's values; repeat as with Vuser above
temp[,`:=`(sddayprev=shift(sdday, n=1, fill=sdday[1L], type="lag"))]
dt<-merge(dt,temp) 

#other calculations:
#sd: bucket-wise (over value per bucket, hence many per day)
setkey(dt,tau,Date)
temp<-dt[,.(sdbucket=sd(deltaP,na.rm=TRUE)),by=.(tau,Date)]
#merge accordingly with the overall info
dt<-merge(dt,temp)#this combines the two data.tables wrt tau and Date


#previous-day calculations: 
#sd: average of bucket-wise (one value per day); also find this value for the previous day
#up to this point, dt had bucket tagging and Date tagging
setkey(dt,Date)
temp<-dt[,.(sdbucketmean=mean(sdbucket,na.rm=TRUE)),by=.(Date)]#Date because we get one estimate per day
temp[,`:=`(sdbucketmeanprev=shift(sdbucketmean, n=1, fill=sdbucketmean[1L], type="lag"))]
dt<-merge(dt,temp)

#Buy (and Sell) volume calculations using different estimates for sigma_{\delta {}}

#use different sigma_{\Delta P} to normalize the dollar price changes and hence different VB

#sdday is a vector so take the first element to scale the price changes
#VB1=sddayprev[1] - normalize by the previous day's sd
#VB2=sdbucketmeanprev[1] - normalizing by the average bucket-wise sd of the previous day
#VB3= no argument; normalizing by the current bucket SD, by default
#VB4=sdday[1] - normalize by the current day's sd - impossible (info unavailable in real time)
#VB5=sdbucketmean[1], normalizing by the average bucket-wise sd of the corresponding day - impossible (info unavailable in real time)

#retain the tagging
dt<-dt[,.(j=tail(j,1),i=tail(i,1),DateTime=tail(DateTime,1),
          Time=tail(Time,1),yday=tail(yday,1),ydayminus1=tail(ydayminus1,1),Vuser=tail(Vuser,1),
          Vuserprev=tail(Vuserprev,1),
          P=tail(P,1),Pnorm=tail(Pnorm,1),Vnorm=tail(Vnorm,1),
          sdday=tail(sdday,1),
          sdbucket=tail(sdbucket,1),sdbucketmean=tail(sdbucketmean,1),
          VB1=sum(BVrule(pricediff=deltaP,sigmadiff=sddayprev[1]),na.rm=TRUE),#prev, one day sd
          VB2=sum(BVrule(pricediff=deltaP,sigmadiff=sdbucketmeanprev[1]),na.rm=TRUE),#prev, average of bucket-sd
          VB3=sum(BVrule(pricediff=deltaP),na.rm=TRUE),#current bucket sd
          VB4=sum(BVrule(pricediff=deltaP,sigmadiff=sdday[1]),na.rm=TRUE),#current, one day sd, infeasible
          VB5=sum(BVrule(pricediff=deltaP,sigmadiff=sdbucketmean[1]),na.rm=TRUE)##prev, average of bucket-sd, infeasible
),by=.(tau,Date)]#by tau index, ie, volume bucket of size Vuser AND by Date

#some time-bars are spread over several volume-buckets, because of 
#a lot of volume being traded in that particular time bar 

#find sell volume accordingly, for each bucket tau
dt[,`:=`(VS1=Vuserprev-VB1,VS2=Vuserprev-VB2,VS3=Vuserprev-VB3,
         VS4=Vuserprev-VB4,VS5=Vuserprev-VB5)]

#find (proportinal) order imbalance, for each bucket tau
dt[,`:=`(OI1=abs(VS1-VB1)/Vuserprev,OI2=abs(VS2-VB2)/Vuserprev,OI3=abs(VS3-VB3)/Vuserprev
         ,OI4=abs(VS4-VB4)/Vuserprev,OI5=abs(VS5-VB5)/Vuserprev)]#per bucket

#obtain VPIN by applying moving average smoothing to OI (withing a given date)
dt<-dt[,.(tau,j,i,DateTime,Time,yday,P,Pnorm,Vnorm,Vuser,Vuserprev,sdday,sdbucket,
          sdbucketmean,
          VB1,VS1,OI1,OI2,OI3,OI4,OI5,
          VPIN1=frollmean(x=OI1,n=SMOOTH),
          VPIN2=frollmean(x=OI2,n=SMOOTH),
          VPIN3=frollmean(x=OI3,n=SMOOTH),
          VPIN4=frollmean(x=OI4,n=SMOOTH),
          VPIN5=frollmean(x=OI5,n=SMOOTH)
          ),by=.(Date)]#by tau index, ie, volume bucket of size Vuser

#find CDF(VPIN) (withing a given date)

#add scaled sdbucket!!! this is done wrt to tau???


dt<-dt[,.(tau,j,i,DateTime,Time,yday,P,Pnorm,Vnorm,Vuser,sdday,sdbucket,sdbucketmean,
          VB1,VS1,OI1,OI2,OI3,OI4,OI5,
          VPIN1,VPIN2,VPIN3,VPIN4,VPIN5,
          cdfVPIN1=mycdf(VPIN1),
          cdfVPIN2=mycdf(VPIN2),
          cdfVPIN3=mycdf(VPIN3),
          cdfVPIN4=mycdf(VPIN4),
          cdfVPIN5=mycdf(VPIN5),
          sdbucketnorm=myscale(sdbucket)),by=.(Date)]#

#add index between 1-390, for computing comovement


myticker<-substr(myfile,1,nchar(myfile)-4)
dt[,`:=`(ticker=myticker)]#tag with ticker

DT<-rbind(DT,dt)

}#end for l loop over files/assets

#plotting

mydays<-levels(as.factor(DT$Date))
t0<-as.POSIXct(paste("2050-01-01","12:00:00"),format="%Y-%m-%d %H:%M:%S",tz="EST")
#ggplot(data=DT[Date %in% mydays[-1],],aes(x=Time,y=OI1,col=ticker))+geom_smooth()+facet_wrap(~Date,nrow=2)
#ggplot(data=DT[Date %in% mydays[-1],],aes(x=Time,y=OI1,col=ticker))+geom_line()+geom_vline(xintercept=t0)+facet_wrap(~Date,nrow=2)
#ggplot(data=DT[Date %in% mydays[-1],],aes(x=i,y=sdbucketnorm,col=ticker))+geom_line()+ facet_wrap(~Date)+coord_cartesian(ylim=c(0,1))
#ggplot(data=DT[Date %in% mydays[-1],],aes(x=tau,y=Pnorm,col=ticker))+geom_line()+ facet_wrap(~Date)#+coord_cartesian(ylim=c(0,1))
#ggplot(data=DT[Date %in% mydays[-1],],aes(x=ticker,y=sdday))+geom_point()+ facet_wrap(~Date,nrow=2)#+coord_cartesian(ylim=c(0,1))
#ggplot(data=DT[Date %in% mydays[-1],],aes(x=i,y=VPIN1,col=ticker))+geom_line()+facet_wrap(~Date,nrow=2)
#ggplot(data=DT,aes(x=yday,y=Vuser,col=ticker))+geom_point()


#customized graphs 
#customized graphs 
#customized graphs 

variable.list<-c("sdbucketnorm",
                 "OI1","OI2","OI3","OI4","OI5",
                 "VPIN1","VPIN2","VPIN3","VPIN4","VPIN5",
                 "cdfVPIN1","cdfVPIN2","cdfVPIN3","cdfVPIN4","cdfVPIN5",
                 "Pnorm","Vnorm")

DF<-NULL


for(k in 1:length(variable.list))
{  
  
  VARIABLE<-variable.list[k]
  
  mydays<-levels(as.factor(DT$Date) )
  tpma<-tpt<-XX<-TT<-vector('list',length=length(mydays))
  
  
  #no need to have tpma here because it only has one slot (x- and y-values)
  for(iii in 1:length(mydays))
    tpt[[iii]]<-XX[[iii]]<-TT[[iii]]<-vector('list',length=length(selected))
  
  xout<-1:390
  
  #loop over tickers and days
  
  for(ii in 1:length(selected))
  {
    myfile<-myfiles[selected[ii]]
    myticker<-substr(myfile,1,nchar(myfile)-4)
    
    temp<-DT[ticker %in% myticker,]
    
    for(iii in 1:length(mydays))
    {    
      XX[[iii]][[ii]]<-as.vector(unlist(temp[Date %in% mydays[iii],][[VARIABLE]]))
      TT[[iii]][[ii]]<-as.vector(unlist(temp[Date %in% mydays[iii],"i"]))
    }
  }  
  
  
  ######################TPMA
  
  for(iii in 1:length(mydays))
  {    
    tpma[[iii]]<-TPMA(X=XX[[iii]],TIME=TT[[iii]],xout=xout,tau=TAU)#TAU declared above
  }
  
  df<-NULL
  
  for(iii in 1:length(mydays))
  {  
    dftemp<-data.frame(TPMA=tpma[[iii]]$rho,
                       xout=tpma[[iii]]$xout,
                       date=mydays[iii])
    df<-rbind(df,dftemp)
  }
  
  df$variable<-VARIABLE
  DF<-rbind(DF,df)#combine various comovement results into one DF
  
  pp<-ggplot(data=df[df$date %in% mydays[-1],],aes(x=xout,y=TPMA))+geom_line()+theme_bw()+facet_wrap(~date,nrow=2)
  pp<-pp+scale_x_continuous(name='Index of a 1-min interval in a trading day, 09:30h-16:00h')
  pp<-pp+coord_cartesian(ylim=c(-1,1))
  pp<-pp+geom_vline(xintercept=c(1,151,301,390),linetype=2)
  pp<-pp+geom_hline(yintercept=c(0),linetype=2)
  pp<-pp+geom_text(x=1,y=1,label='09:30h',size=3)
  pp<-pp+geom_text(x=151,y=1,label='12:00h',size=3)
  pp<-pp+geom_text(x=301,y=1,label='14:30h',size=3)
  pp<-pp+geom_text(x=390,y=1,label='15:59h',size=3)
  
  #pp<-pp+geom_vline(xintercept=c(303),linetype=1)
  #pp<-pp+geom_text(x=303,y=-1,label='14:32h',size=3)#Flash Crash
  
  #pp<-pp+geom_text(x=303,y=0.85,label='Flash Crash')
  pp<-pp+scale_y_continuous(name=paste0('TPMA of ',VARIABLE,' (thickness=',TAU,')'))
  #pp
  
  ###configuration that works well:
  #cdfvipin (smoothing=10,thickness=30)
  
  
  #even better: sdbucket, scaled in real-time using min-max
  #tau=15 (100 buckets per day); issue: Vuser is forward-looking 
  
  print(pp)
  
}#end of k loop over VARIABLE in variable list


# make one plot of all comovement results
#all together
selected.variables<-c("sdbucketnorm",
                      "Pnorm","Vnorm",
                      "OI1","OI2","OI3"
)

#drop the first day too

temp<-DF[DF$variable %in% selected.variables & DF$date %in% mydays[-1],]

temp$Variable<-factor(temp$variable)
temp$Variable<-factor(temp$Variable,levels=levels(temp$Variable),
                      labels=c("OI1","OI2","OI3","Price","SD bucket","Volume") )

pp<-ggplot(data=temp,aes(x=xout,y=TPMA,col=Variable))+geom_line()+facet_wrap(~date,nrow=2)
pp<-pp+scale_x_continuous(name='Index of a 1-min interval in a trading day, 09:30h-16:00h')
pp<-pp+coord_cartesian(ylim=c(-1,1))
pp<-pp+geom_vline(xintercept=c(1,151,301,390),linetype=2)
pp<-pp+geom_hline(yintercept=c(0),linetype=2)
pp<-pp+geom_text(x=1,y=1,label='09:30h',size=3,col=1)
pp<-pp+geom_text(x=151,y=1,label='12:00h',size=3,col=1)
pp<-pp+geom_text(x=301,y=1,label='14:30h',size=3,col=1)
pp<-pp+geom_text(x=390,y=1,label='15:59h',size=3,col=1)
pp<-pp+theme_bw()
pp3<-pp
print(pp3)

