# Model 1
Model1 = function(testET) {
  Model1_ET = rep(NA,nrow(testET))
  par(mfrow=c(4,4),oma=c(2,2,2,2),mar=c(2,2,1,1))
  
  for(p in 1:length(unique(testET$date))) {
    midpt_daytime = sum(range(testET$hour[which(testET$date==unique(testET$date)[p])]))/2 
    sunrise = min(testET$hour[which(testET$date==unique(testET$date)[p])])
    sunset = max(testET$hour[which(testET$date==unique(testET$date)[p])])
    peakdiff = 13.5-midpt_daytime # peak occurs 1:30pm - shift from daytime centre
    cycday = (midpt_daytime-sunrise+peakdiff)*4 # update sinusoidal period based on peak time
    
    ETcal = testET$ET[which(testET$date==unique(testET$date)[p])]
    hourcal = testET$hour[which(testET$date==unique(testET$date)[p])]-
      min(testET$hour[which(testET$date==unique(testET$date)[p])]) # hour starting from sunrise each day
    
    plot(ETcal~hourcal,
         main=unique(testET$date)[p],ylim=c(0,0.7),xaxt="n",pch=20,col="blue",xlab="Hour")
    axis(1,at=seq(min(testET$hour[which(testET$date==unique(testET$date)[p])]),
                  max(testET$hour[which(testET$date==unique(testET$date)[p])]),by=1)-
           min(testET$hour[which(testET$date==unique(testET$date)[p])])+1,
         label=seq(min(testET$hour[which(testET$date==unique(testET$date)[p])]),
                   max(testET$hour[which(testET$date==unique(testET$date)[p])]),by=1))
    
    sine_eachday <- sin(2*pi*hourcal/cycday)
    
    fit_eachday = lm(ETcal ~ sine_eachday+0) # fit sinusoidal function to available daytime ETa records 
    ETa_pred = predict(fit_eachday,data.frame(sine_eachday=sin(2*pi*hourcal/cycday))) # predict all daytime ETa records
    
    points(y=ETa_pred,x=hourcal,col="red",pch=20,type="o")
    
    Model1_ET[which(testET$date==unique(testET$date)[p])] = ETa_pred
    
    if (p == length(unique(testET$date))) {
      legend("topleft",pch=20,col=c("blue","red"),legend=c("obs. ETa","fitted ETa"))
      mtext(side=1,"Hour",outer=T,line=0);mtext(side=2,"ETa (mm)",outer=T,line=0)
      mtext(side=3,"Infill Method 1 - Testing",outer=T,font=2,cex=1.4)
    }
    
  }
  return(Model1_ET)
}

# Model 2
Model2 = function(testET) {
  Model2_ET = rep(NA,nrow(testET))
  par(mfrow=c(4,4),oma=c(2,2,2,2),mar=c(2,2,1,1))
  
  for(p in 1:length(unique(testET$date))) {
    ETcal = testET$ET[which(testET$date==unique(testET$date)[p])]
    hourcal = testET$hour[which(testET$date==unique(testET$date)[p])]-
      min(testET$hour[which(testET$date==unique(testET$date)[p])]) # hour starting from sunrise each day
    
    plot(ETcal~hourcal,
         main=unique(testET$date)[p],ylim=c(0,0.7),xaxt="n",pch=20,col="blue",xlab="Hour")
    axis(1,at=seq(min(testET$hour[which(testET$date==unique(testET$date)[p])]),
                  max(testET$hour[which(testET$date==unique(testET$date)[p])]),by=1)-
           min(testET$hour[which(testET$date==unique(testET$date)[p])])+1,
         label=seq(min(testET$hour[which(testET$date==unique(testET$date)[p])]),
                   max(testET$hour[which(testET$date==unique(testET$date)[p])]),by=1))
    
    
    fit_eachday = lm(ETcal ~ poly(hourcal,2)) # fit 2-deg polynomial to available daytime ETa records 
    ETa_pred = predict(fit_eachday, data.frame(hourcal = hourcal)) # predict all daytime ETa records
    
    points(y=ETa_pred,x=hourcal,col="red",pch=20,type="o")
    
    Model2_ET[which(testET$date==unique(testET$date)[p])] = ETa_pred
    
    if (p == length(unique(testET$date))) {
      legend("topleft",pch=20,col=c("blue","red"),legend=c("obs. ETa","fitted ETa"))
      mtext(side=1,"Hour",outer=T,line=0);mtext(side=2,"ETa (mm)",outer=T,line=0)
      mtext(side=3,"Infill Method 2 - Testing",outer=T,font=2,cex=1.4)
    }
    
  }
  return(Model2_ET)
}

# Model 3
Model3 = function(testET) {
  Model3_ET = rep(NA,nrow(testET))
  par(mfrow=c(4,4),oma=c(2,2,2,2),mar=c(2,2,1,1))
  
  for(p in 1:length(unique(testET$date))) {
    ETcal = testET$ET[which(testET$date==unique(testET$date)[p])]
    ET0cal = testET$ET0[which(testET$date==unique(testET$date)[p])]
    
    hourcal = testET$hour[which(testET$date==unique(testET$date)[p])]-
      min(testET$hour[which(testET$date==unique(testET$date)[p])]) # hour starting from sunrise each day
    
    plot(ETcal~hourcal,
         main=unique(testET$date)[p],ylim=c(0,0.7),xaxt="n",pch=20,col="blue",xlab="Hour")
    axis(1,at=seq(min(testET$hour[which(testET$date==unique(testET$date)[p])]),
                  max(testET$hour[which(testET$date==unique(testET$date)[p])]),by=1)-
           min(testET$hour[which(testET$date==unique(testET$date)[p])])+1,
         label=seq(min(testET$hour[which(testET$date==unique(testET$date)[p])]),
                   max(testET$hour[which(testET$date==unique(testET$date)[p])]),by=1))
    
    ETratiocal =  ETcal/ET0cal
    
    if (!all(is.na(ETratiocal))) {
      if(any(ETratiocal[which(!is.na(ETratiocal))]>5)) {
        ETratiocal[which(!is.na(ETratiocal))][which(ETratiocal[!is.na(ETratiocal)]>5)] = NA
      }
      
      if(any(ETratiocal[which(!is.na(ETratiocal))]<0)) {
        ETratiocal[which(!is.na(ETratiocal))][which(ETratiocal[!is.na(ETratiocal)]<0)] = NA
      }
      
      fit_eachday = lm(ETratiocal ~ poly(hourcal,2))
      
      ETa_pred = predict(fit_eachday, data.frame(hourcal = hourcal))#xc=cos(2*pi*(seq(1,11,by=0.5))/15-7),
      #ETa_pred[which(ETa_pred<0)]=0
      Model3_ET[which(testET$date==unique(testET$date)[p])] = ETa_pred*ET0cal
      
    }
    points(y=ETa_pred*ET0cal,x=hourcal,col="red",pch=20,type="o")
    
    if (p == length(unique(testET$date))) {
      legend("topleft",pch=20,col=c("blue","red"),legend=c("obs. ETa","fitted ETa"))
      mtext(side=1,"Hour",outer=T,line=0);mtext(side=2,"ETa (mm)",outer=T,line=0)
      mtext(side=3,"Infill Method 3 - Testing",outer=T,font=2,cex=1.4)
    }
    
  }
  return(Model3_ET)
}

# Model 4
Model4 = function(testET,trainET) {
  Model4_ET = rep(NA,nrow(testET))
  par(mfrow=c(4,4),oma=c(2,2,2,2),mar=c(2,2,1,1))
  
  for(p in 1:length(unique(testET$date))) {
    ETcal = testET$ET[which(testET$date==unique(testET$date)[p])]
    hourcal = testET$hour[which(testET$date==unique(testET$date)[p])]-
      min(testET$hour[which(testET$date==unique(testET$date)[p])]) # hour starting from sunrise each day
    
    plot(ETcal~hourcal,
         main=unique(testET$date)[p],ylim=c(0,0.7),xaxt="n",pch=20,col="blue",xlab="Hour")
    axis(1,at=seq(min(testET$hour[which(testET$date==unique(testET$date)[p])]),
                  max(testET$hour[which(testET$date==unique(testET$date)[p])]),by=1)-
           min(testET$hour[which(testET$date==unique(testET$date)[p])])+1,
         label=seq(min(testET$hour[which(testET$date==unique(testET$date)[p])]),
                   max(testET$hour[which(testET$date==unique(testET$date)[p])]),by=1))
    
    cortrain = rep(NA,length(unique(trainET$date)))
    for (q in c(1:length(unique(trainET$date)))) { # calculate all pearson correlations between this day and each day with complete data
      if (unique(testET$date)[p]!=unique(trainET$date)[q]) { 
        trainsel = which(trainET$date==unique(trainET$date)[q])
        if (all(!is.na(trainET$ET[trainsel]))) {
          commtesttrain = intersect(testET$hour[which(testET$date==unique(testET$date)[p])],
                                    trainET$hour[trainET$date==unique(trainET$date)[q]])
          whichtest = which(testET$hour[which(testET$date==unique(testET$date)[p])]%in%commtesttrain)
          whichtrain = which(trainET$hour[trainET$date==unique(trainET$date)[q]]%in%commtesttrain)
          cortrain[q] = cor(testET$ET[which(testET$date==unique(testET$date)[p])][whichtest], 
                            trainET$ET[trainsel][whichtrain],
                            use="pairwise.complete.obs")
        }
        
      }
    }
    
    matchdate = unique(trainET$date)[which(cortrain==max(cortrain,na.rm=T))] # find the day with complete data that has the highest correlation with the current day
    commtesttrain = intersect(testET$hour[which(testET$date==unique(testET$date)[p])],
                              trainET$hour[trainET$date==matchdate]) # find the common hours of the two days (since individual days have different daytime lengths included in the dataset)
    whichtest = which(testET$hour[which(testET$date==unique(testET$date)[p])]%in%commtesttrain) # get hour index of the day to fill
    whichtrain = which(trainET$hour[trainET$date==matchdate]%in%commtesttrain)   # get hour index of the matching day (with complete data)
    naind = which(is.na(testET$ET[which(testET$date==unique(testET$date)[p])][whichtest]))
    
    # find which hours in train matches test
    ETa_frac = trainET$ET[which(trainET$date==matchdate)][whichtrain]/sum(trainET$ET[which(trainET$date==matchdate)][whichtrain],na.rm=T) # get the ratio of 30-min ETa to total daily ETa
    ETa_frac_total = sum(trainET$ET[which(trainET$date==matchdate)][-naind])/sum(trainET$ET[which(trainET$date==matchdate)],na.rm=T) # calculate the expected proportion of daily total taken by the non-NA slots 
    ETa_sum_day = sum(testET$ET[which(testET$date==unique(testET$date)[p])][whichtest],na.rm=T)/ETa_frac_total # estimated the expected daily total for the day to infill
    
    
    ETa_pred = ETa_sum_day*ETa_frac 
    
    Model4_ET[which(testET$date==unique(testET$date)[p])][whichtest] = ETa_pred
    points(y=ETa_pred,x=hourcal[whichtest],
           col="red",pch=20,type="o")
    
    if (p == length(unique(testET$date))) {
      legend("topleft",pch=20,col=c("blue","red"),legend=c("obs. ETa","fitted ETa"))
      mtext(side=1,"Hour",outer=T,line=0);mtext(side=2,"ETa (mm)",outer=T,line=0)
      mtext(side=3,"Infill Method 4 - Testing",outer=T,font=2,cex=1.4)
    }
  }
  return(Model4_ET)
}
