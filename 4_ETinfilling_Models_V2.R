# Model 1 - Daily sinusoidal functions of ETa (Sinusoidal)
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
      mtext(side=3,"Infill Method 1 (Sinusoidal) - Testing",outer=T,font=2,cex=1.4)
    }
    
  }
  return(Model1_ET)
}

# Model 2 - Daily smoothing functions of ETa (Smoothing)
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
      mtext(side=3,"Infill Method 2 (Smoothing) - Testing",outer=T,font=2,cex=1.4)
    }
    
  }
  return(Model2_ET)
}

# Model 3 - Daily temporal pattern matching for ETa (MaxCor)
Model3 = function(testET,trainET) {
  Model3_ET = rep(NA,nrow(testET))
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
    
    Model3_ET[which(testET$date==unique(testET$date)[p])][whichtest] = ETa_pred
    points(y=ETa_pred,x=hourcal[whichtest],
           col="red",pch=20,type="o")
    
    if (p == length(unique(testET$date))) {
      legend("topleft",pch=20,col=c("blue","red"),legend=c("obs. ETa","fitted ETa"))
      mtext(side=1,"Hour",outer=T,line=0);mtext(side=2,"ETa (mm)",outer=T,line=0)
      mtext(side=3,"Infill Method 3 (MaxCor) - Testing",outer=T,font=2,cex=1.4)
    }
  }
  return(Model3_ET)
}

# Model 4 - Mean Diurrnal Variation (MDV)
Model4 = function(testET,trainET) {
  Model4_ET = rep(NA,nrow(testET))
  par(mfrow=c(4,4),oma=c(2,2,2,2),mar=c(2,2,1,1))
  
  for(p in 1:length(unique(c(testET$date,trainET$date)))) { # Find the 14-day window starting from each day in the season
    datecut = seq(unique(c(testET$date,trainET$date))[p],unique(c(testET$date,trainET$date))[p]+13,
                  "days")
    fillind = which(is.na(Model4_ET)&!is.na(testET_true$ET)&testET$date%in%datecut) # find which data points are gaps (to infill)
    fillhour = unique(testET$hour[fillind]) # get the hours of gaps to infill
    for (h in 1:length(fillhour)) { # for each hour to infill, get the mean value within the 14-day window to infill gaps
      Model4_ET[which(testET$date%in%datecut&testET$hour==fillhour[h])] = 
        mean(c(trainET$ET[which(trainET$date%in%datecut&trainET$hour==fillhour[h])]),na.rm=T)
    }
    
  }
  Model4_ET[which(is.na(Model4_ET))]=0 # For a few data points at the start or end of the day there is no matching hour within 14 days - zeroed
  
  for(p in 1:length(unique(testET$date))) { # plot the available and infilled data for each day
    ETuse = testET$ET[which(testET$date==unique(testET$date)[p])]
    Model4_ET_sel = Model4_ET[which(testET$date==unique(testET$date)[p]&
                                is.na(testET$ET)&!is.na(testET_true$ET))]
    thouruse = testET$hour[which(testET$date==unique(testET$date)[p])]-
      min(testET$hour[which(testET$date==unique(testET$date)[p])]) # as we start from sunrise each day
    
    plot(ETuse~
           thouruse,
         main=unique(testET$date)[p],ylim=c(0,0.7),xaxt="n",pch=20,col="blue",xlab="Hour")
    axis(1,at=seq(min(testET$hour[which(testET$date==unique(testET$date)[p])]),
                  max(testET$hour[which(testET$date==unique(testET$date)[p])]),by=1)-
           min(testET$hour[which(testET$date==unique(testET$date)[p])])+1,
         label=seq(min(testET$hour[which(testET$date==unique(testET$date)[p])]),
                   max(testET$hour[which(testET$date==unique(testET$date)[p])]),by=1))
    
    points(y=Model4_ET_sel,x=thouruse[which(is.na(ETuse)&
                                           !is.na(testET_true$ET[which(testET$date==unique(testET$date)[p])]))],col="red",pch=20,type="o")
    
    if (p == length(unique(testET$date))) {
      legend("topleft",pch=20,col=c("blue","red"),legend=c("obs. ETa","fitted ETa"))
      mtext(side=1,"Hour",outer=T,line=0);mtext(side=2,"ETa (mm)",outer=T,line=0)
      mtext(side=3,"Infill Method 4 (MDV) - Testing",outer=T,font=2,cex=1.4)
    }
  }
  
  
  return(Model4_ET)
}
