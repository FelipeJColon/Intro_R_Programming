
rm(list=ls(all=TRUE))

# Load packages
require(data.table)
require(dplyr)
require(tidyr)
require(surveillance)
require(lubridate)
require(zoo)

# Specify directory
myDir <- "/local/zck07apu/Documents/GitLab/Bayesian_factor/2x/"

# -------------------------
# Define global vars
# -------------------------
nweeks           <- 52
nweeks2          <- 49
ndays5 <- days5  <- 5
ndays7 <- days7  <- 7
nsim             <- 100
narrays          <- 17
nyears <- years  <- 7

# ----------------------
# Produce data
# ----------------------

# 5-day systems

h1 <- function(N, k, k2, alpha, beta, gama1,
               gama2, gama3, gama4, shift, shift2){
     t=1:N
     if(k==0 & k2==0){h1=alpha+beta*t}
     else{
          if(k==0){
               l=1:k2
               h1=rep(0,N)
               for(i in 1:N){
                    h1[i]=alpha + beta * (t[i] + shift) + 
                         sum(gama3 * cos((2 * pi * l * (t[i] + shift)) / 5) + 
                                  gama4 * sin((2 * pi * l * (t[i] + shift)) / 5))
               }
          }
          else{
               j=1:k
               l=1:k2
               h1=rep(0,N)
               for(i in 1:N){
                    h1[i]=alpha + beta * (t[i] + shift) + 
                         sum(gama1 * cos((2 * pi * j * (t[i] + shift)) / 
                                              (nweeks *5)) + gama2 * 
                                  sin((2 * pi * j * (t[i]+shift2))/(nweeks * 5))) + 
                         sum(gama3 * cos((2 * pi * l * (t[i] + shift)) / 5) + 
                                  gama4 * sin((2 * pi * l * (t[i] + shift))/5))
               }
          }
     }
     h1
}

negbinNoise1 <- function(N, k, k2, alpha, beta, gama1, 
                         gama2, gama3, gama4, phi, shift, shift2){
     mu=exp(h1(N, k, k2, alpha, beta, gama1, gama2, gama3,
               gama4, shift, shift2))
     if(phi==1){yi=rpois(N, mu)}
     else{
          prob = 1/phi 
          size = mu/(phi-1) 
          yi   = rnbinom(N,size=size,prob=prob)
     }
     yi
}


outbreak5 <- function(currentday, weeklength, wtime, yi, 
                      interval, k, k2, alpha, beta, gama1,
                      gama2, gama3, gama4, shift, shift2, 
                      phi, numoutbk, peakoutbk, meanlog, sdlog){
     N  = length(yi)
     t  = 1:N
     mu = exp(h1(N, k, k2, alpha, beta, gama1, gama2,
                 gama3, gama4, shift, shift2))
     s  = sqrt(mu*phi)
     
     ## -------------------
     ## Generate outbreaks
     ## -------------------
     
     # Start time
     startoutbk = sample(wtime, numoutbk, replace = FALSE)
     
     # Outbreak size
     sizeoutbk = rep(0, numoutbk)
     for(i in 1:numoutbk){
          set.seed(i)
          soutbk=1
          sou=1
          while(soutbk<2){
               set.seed(sou)
               soutbk=rpois(1, s[startoutbk[i]] * peakoutbk)
               sou=sou+1
          }
          sizeoutbk[i]=soutbk
     }
     
     # Distribute cases over time using lognormal function
     outbreak=rep(0,2*N)
     for(j in 1:numoutbk){
          set.seed(j)
          outbk    = rlnorm(sizeoutbk[j], meanlog = meanlog, sdlog = sdlog) 
          h        = hist(outbk,breaks=seq(0, ceiling(max(outbk)), interval),
                          plot=FALSE)
          cases    = h$counts
          weight   = rep(0,length(cases))
          duration = startoutbk:(startoutbk+length(cases)-1)
          dayofweek<-duration%%5 # 0 is friday; 1 is monday; 2 is tuesday etc.
          for(i in 1:length(cases)){
               if(dayofweek[i]==0){weight[i]=1.1}
               if(dayofweek[i]==1){weight[i]=1.5}
               if(dayofweek[i]==2){weight[i]=1.1}
               if(dayofweek[i]==3){weight[i]=1}
               if(dayofweek[i]==4){weight[i]=1}
          }
          cases2 = cases * weight
          for (l in 1:(length(cases2))){
               outbreak[startoutbk[j] + (l-1)]=cases2[l] + 
                    outbreak[startoutbk[j] + (l-1)]
          }
     }
     
     for(v in currentday:(currentday+100)){if(outbreak[v] > 0){outbreak[v]=0}}
     outbreak=outbreak[1:N]
     
     # Add outbreaks and noise
     yitot  = yi + outbreak
     result = list(yitot=yitot, outbreak=outbreak, startoutbk=startoutbk,
                   sizeoutbk=sizeoutbk, sd=s, mean=mu)
}

# 7-day systems
h2 <- function(N, k, k2, alpha, beta, gama1, gama2,
               gama3, gama4, shift){
     t=1:N
     if(k==0 & k2==0){h2=alpha + beta * t}
     else{
          if(k==0)
          {
               l=1:k2
               h2=rep(0, N)
               for(i in 1:N){
                    h2[i]=alpha + beta * (t[i] + shift) + 
                         sum(gama3 * cos((2 * pi * l * (t[i] + shift)) / 7) + 
                                  gama4 * sin((2 * pi * l * (t[i] + shift)) /
                                                   7))
               }
          }
          else{
               j=1:k
               l=1:k2
               h2=rep(0,N)
               for(i in 1:N){
                    h2[i]=alpha + beta * (t[i] + shift) + 
                         sum(gama1 * cos((2 * pi * j * (t[i] + shift)) / 
                                              (nweeks * 7)) + gama2 * 
                                  sin((2 * pi * j * (t[i] + shift)) / 
                                           (nweeks * 7))) + 
                         sum(gama3 * cos((2 * pi * l * (t[i] + shift)) / 7) + 
                                  gama4 * sin((2 * pi * l * (t[i]+shift)) / 7))
               }
          }
     }
     h2
}

negbinNoise2 <- function(N, k, k2, alpha, beta, gama1, gama2,
                         gama3, gama4, phi, shift){
     mu=exp(h2(N, k, k2, alpha, beta, gama1, gama2, gama3, gama4, shift))
     if(phi==1){yi=rpois(N,mu)}
     else{
          prob=1/phi 
          size=mu/(phi-1) 
          yi=rnbinom(N,size=size,prob=prob)
     }
     yi
}

outbreak7 <- function(currentday, weeklength, wtime, yi, 
                      interval, k, k2, alpha, beta, gama1,
                      gama2, gama3, gama4, shift, phi,
                      numoutbk, peakoutbk, meanlog, sdlog){
     N=length(yi)
     t=1:N
     mu=exp(h2(N, k, k2, alpha, beta, gama1, gama2, 
               gama3, gama4, shift))
     s=sqrt(mu*phi)
     
     startoutbk=sample(wtime, numoutbk, replace = FALSE)
     
     sizeoutbk=rep(0,numoutbk)
     for(i in 1:numoutbk){
          set.seed(i)
          soutbk=1
          sou=1
          while(soutbk < 2){
               set.seed(sou)
               soutbk=rpois(1, s[startoutbk[i]] * peakoutbk)
               sou=sou+1
          }
          sizeoutbk[i]=soutbk
     }
     
     outbreak=rep(0, 2 * N)
     for( j in 1:numoutbk){
          set.seed(j)
          outbk     =rlnorm(sizeoutbk[j], meanlog = meanlog, sdlog = sdlog) 
          h         =hist(outbk, breaks=seq(0, ceiling(max(outbk)), interval),
                          plot=FALSE)
          cases     =h$counts
          weight    =rep(0,length(cases))
          duration  =startoutbk:(startoutbk+length(cases)-1)
          dayofweek = duration%%7 # 0 is sunday; 1 is monday; 2 is tuesday etc.
          for(i in 1:length(cases)){
               if(dayofweek[i]==0){weight[i]=2}
               if(dayofweek[i]==1){weight[i]=1}
               if(dayofweek[i]==2){weight[i]=1}
               if(dayofweek[i]==3){weight[i]=1}
               if(dayofweek[i]==4){weight[i]=1}
               if(dayofweek[i]==5){weight[i]=1}
               if(dayofweek[i]==6){weight[i]=2}
          }
          cases2 <- cases*weight
          for (l in 1:(length(cases2))){
               outbreak[startoutbk[j]+(l-1)]= cases2[l]+outbreak[startoutbk[j]+(l-1)]
          }
     }
     
     for(v in currentday:(currentday+100)){if(outbreak[v]>0){outbreak[v]=0}}
     outbreak=outbreak[1:N]
     
     yitot=yi+outbreak
     result=list(yitot=yitot, outbreak=outbreak, startoutbk=startoutbk,
                 sizeoutbk=sizeoutbk, sd=s, mean=mu)
}

# --------------------------------
# Specify bank holidays
# --------------------------------
years        <- 7
bankholidays <- fread(file.path(myDir, "Bankholidays.csv"))

# Fix bug on dow (thur and thu)
bankholidays$dayofweek <- as.character(bankholidays$dayofweek)
bankholidays$dayofweek[bankholidays$dayofweek=="Thur"] <- "Thu"

# Sort factor levels
bankholidays$dayofweek <- factor(bankholidays$dayofweek,
                                 levels=c("Mon", "Tue", "Wed",
                                          "Thu", "Fri", "Sat",
                                          "Sun"))

# Fix bug on month (mar instead of "march")
bankholidays$month <- as.character(bankholidays$month)
bankholidays$month[bankholidays$month=="March"] <- "Mar"

# Order factor levels of month
bankholidays$month <- factor(bankholidays$month, 
                             levels=c("Jan", "Feb", "Mar", "Apr",
                                      "May", "Jun", "Jul", "Aug",
                                      "Sep", "Oct", "Nov", "Dec"))

bankhols7 <- bankholidays$bankhol
bankhols7 <- as.numeric(bankhols7)
length(bankhols7)


bankhols5 <- bankhols7[-seq(6, length(bankhols7), 7)] 
bankhols5 <- bankhols5[-seq(6, length(bankhols5), 6)]
bankhols5 <- as.numeric(bankhols5)
length(bankhols5)

# ----------------------------
# Define the data frames
# for input data
# ----------------------------
# Create lsits of empty matrices
lstData  <- replicate(narrays, matrix(nrow=nweeks * ndays7 * nyears, 
                                      ncol=nsim), 
                      simplify=FALSE)
lsOutbr <- replicate(narrays, matrix(nrow=1, ncol=nsim), 
                     simplify=FALSE)

lstzOutbreak <- lstOutbreak <- lstTotal <- lstData

probOut <- lsOutbr

# Assign names to elements of lists
names(lstData)      <- paste0("simulateddata", 1:narrays)
names(lstTotal)     <- paste0("simulatedtotals", 1:narrays)
names(lstOutbreak)  <- paste0("simulatedoutbreak", 1:narrays)
names(lstzOutbreak) <- paste0("simulatedzseasoutbreak", c(6, 7, 16))
names(lsOutbr)      <- paste0("outbreaklength", 1:narrays)
names(probOut)      <- paste0("outbreakprior", 1:narrays)

# Split lists into single matrices
list2env(lstData, .GlobalEnv)
list2env(lstTotal, .GlobalEnv)
list2env(lstOutbreak, .GlobalEnv)
list2env(lstzOutbreak, .GlobalEnv)
list2env(lsOutbr, .GlobalEnv)
list2env(probOut, .GlobalEnv)

# --------------------------------------
# Simulate syndromes and outbreaks
# --------------------------------------

# 5-day week syndromes

days5 <- 5
N     <- nweeks * days5 * years

# sigid6
for(i in 1:nsim){
     set.seed(i)
     
     # Generate syndromic data
     yt = round(negbinNoise1(N=N, k=1, k2=1, alpha=6, beta=0, gama1=0.3,
                             gama2=2, gama3=0.3, gama4=0.5, phi=1.5,
                             shift=-50, shift2=-50)/10)
     out1=rep(0, N)
     
     # Compute seasonal outbreak data
     for(j in 1:years){
          set.seed(j + years * i)
          out=outbreak5(currentday=days5 * nweeks * years, 
                        weeklength=nweeks * days5 * years,
                        wtime=((1 + (j - 1) * days5 * 
                                     nweeks):(20 +(j - 1) * days5 * nweeks)),
                        yi=yt, interval=0.02, k=1,
                        k2=1, alpha=6, beta=0, gama1=0.3, gama2=2,
                        gama3=0.3, gama4=0.5, phi=1.5, shift=-50, shift2=-50, 
                        numoutbk=1, peakoutbk=3 * days5 * 80, meanlog=0, 
                        sdlog=0.5)
          out1=out1 + out$outbreak
     }
     
     out1 = round(out1)
     
     # Compute outbreak2 - non-seasonal
     set.seed(i)
     out2=outbreak5(currentday=days5 * nweeks * years, 
                    weeklength=nweeks * days5 * years,
                    wtime=(length(yt) - 49 * days5 + 1):length(yt),
                    yi=yt, interval=0.25, k=1, k2=1, alpha=6,
                    beta=0, gama1=0.3, gama2=2, gama3=0.3, gama4=0.5,
                    phi=1.5, shift=-50, shift2=-50, numoutbk=1,
                    peakoutbk=2 * days5, meanlog=0, sdlog=0.5)
     
     # Non-seasonal outbreak
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength6[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior6[,i] = outbreaklength6[,i] / (364 * nyears)
     
     # Seasonal outbreak
     zseasoutbreak=zoutbreak + out1
     
     # Seasonal syndromic data + seasonal outbreak
     zt = yt + out1
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     # Set bank holidays to zero
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     # Append weekends to 5-day time series
     zeros=rep(0, 2)
     weekend=seq(days5, days5 * years * nweeks, days5)
     
     for(s in 1:length(weekend)){
          zt=append(zt, zeros, after=2 * (s - 1) + weekend[s])
          zitot=append(zitot, zeros, after=2 * (s - 1) + weekend[s])
          zoutbreak=append(zoutbreak, zeros, after=2 * (s - 1) + weekend[s])
          zseasoutbreak=append(zseasoutbreak, zeros, after=2 * (s - 1) + 
                                    weekend[s])
     }
     
     ## Output
     
     # Data with just seasonal outbreak
     simulateddata6[,i] = zt
     
     # Data plus non-seasonal outbreak
     simulatedtotals6[,i] = zitot
     
     # Non-seasonal outbreak data
     simulatedoutbreak6[,i] = zoutbreak
     
     simulatedzseasoutbreak6[,i]=zseasoutbreak
     
}

#sigid7
for(i in 1:nsim){
     set.seed(i)
     yt=round(negbinNoise1(N=N, k=1, k2=1, alpha=1, beta=0,
                           gama1=0.1, gama2=2, gama3=0.05,
                           gama4=0.05, phi=1, shift=-50, shift2=-50))
     out1=rep(0, N)
     for(j in 1:years){
          set.seed(j+years*i)
          out=outbreak5(currentday=days5 * nweeks * years,
                        weeklength=nweeks * days5 * years,
                        wtime=((1+(j-1) * days5 * nweeks):
                                    (20+(j-1) * days5 * nweeks)),
                        yi=yt, interval=0.02, k=1, k2=1, alpha=1, 
                        beta=0, gama1=0.1, gama2=2, gama3=0.05,
                        gama4=0.05, phi=1, shift=-50, shift2=-50,
                        numoutbk=1, peakoutbk=3 * days5 * 50,
                        meanlog=0, sdlog=0.5)
          out1=out1+out$outbreak
     }
     
     out1 = round(out1)
     
     set.seed(i)
     out2=outbreak5(currentday=days5 * nweeks * years,
                    weeklength=nweeks * days5 * years,
                    wtime=(length(yt)-49 * days5 + 1):length(yt),
                    yi=yt, interval=0.25, k=1, k2=1, alpha=1,
                    beta=0, gama1=0.1, gama2=2, gama3=0.05, 
                    gama4=0.05, phi=1, shift=-50, shift2=-50,
                    numoutbk=1, peakoutbk=2 * days5, meanlog=0,
                    sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength7[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior7[,i] = outbreaklength7[,i] / (364 * nyears)
     
     zseasoutbreak=zoutbreak+out1
     
     # Seasonal syndromic data + seasonal outbreak
     zt = yt + out1
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     # Set bank holidays to zero
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     zeros=rep(0,2)
     weekend=seq(days5, days5 * years * nweeks, days5)
     for(s in 1:length(weekend)){
          yt=append(yt, zeros, after=2 * (s-1) + weekend[s])
          zt=append(zt, zeros, after=2 * (s-1) + weekend[s])
          zitot=append(zitot, zeros, after=2 * (s - 1) + weekend[s])
          zoutbreak=append(zoutbreak, zeros, after=2 * (s - 1) + weekend[s])
          zseasoutbreak=append(zseasoutbreak, zeros, after=2 * (s-1) + 
                                    weekend[s])
     }
     
     simulateddata7[,i]=zt
     simulatedtotals7[,i]=zitot
     simulatedoutbreak7[,i]=zoutbreak
     simulatedzseasoutbreak7[,i]=zseasoutbreak
}

#sigid8
for(i in 1:nsim){
     set.seed(i)
     yt=round(negbinNoise1(N=N, k=0, k2=1, alpha=6, beta=0.0001,
                           gama1=0, gama2=0, gama3=0.6, gama4=0.9,
                           phi=1.5, shift=0, shift2=0)/10)
     
     set.seed(i)
     out2=outbreak5(currentday=days5 * nweeks * years,
                    weeklength=nweeks * days5 * years,
                    wtime=(length(yt)-49 * days5 + 1):length(yt),
                    yi=yt, interval=0.25, k=0, k2=1, alpha=6,
                    beta=0, gama1=0, gama2=0, gama3=0.6, 
                    gama4=0.9, phi=1.5, shift=0, shift2=0,
                    numoutbk=1, peakoutbk=2 * days5, meanlog=0,
                    sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength8[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior8[,i] = outbreaklength8[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     zeros=rep(0, 2)
     weekend=seq(days5, days5 * years * nweeks, days5)
     for(s in 1:length(weekend)){
          zt=append(zt, zeros, after=2 * (s-1) + weekend[s])
          zitot=append(zitot, zeros, after=2 * (s-1) + weekend[s])
          zoutbreak=append(zoutbreak, zeros, after=2 * (s-1) + weekend[s])
     }
     
     simulateddata8[,i]=zt
     simulatedtotals8[,i]=zitot
     simulatedoutbreak8[,i]=zoutbreak
}

#sigid9
for(i in 1:nsim){
     
     set.seed(i)
     yt=round(negbinNoise1(N=N, k=1, k2=1, alpha=3, beta=0,
                           gama1=1.5, gama2=0.1, gama3=0.2,
                           gama4=0.3, phi=1, shift=-150,
                           shift2=-150))
     
     mu=exp(h1(N=N, k=1, k2=1, alpha=3, beta=0, gama1=1.5, 
               gama2=0.1, gama3=0.6, gama4=0.8, shift=-150,
               shift2=-150))
     
     set.seed(i)
     out2=outbreak5(currentday=days5 * nweeks * years, 
                    weeklength=nweeks * days5 * years,
                    wtime=(length(yt)-49 * days5+1):length(yt),
                    interval=0.25, yi=yt, k=1, k2=1, alpha=3,
                    beta=0, gama1=1.5, gama2=0.1, gama3=0.2,
                    gama4=0.3, phi=1, shift=-150, shift2=-150,
                    numoutbk=1, peakoutbk=2 * days5, meanlog=0,
                    sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength9[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior9[,i] = outbreaklength9[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     zeros=rep(0, 2)
     weekend=seq(days5, days5 * years * nweeks, days5)
     for(s in 1:length(weekend)){
          zt=append(zt, zeros, after=2 * (s - 1) + weekend[s])
          zitot=append(zitot, zeros, after=2 * (s - 1) + weekend[s])
          zoutbreak=append(zoutbreak, zeros, after=2 * (s - 1) + weekend[s])
     }
     
     simulateddata9[,i]=zt
     simulatedtotals9[,i]=zitot
     simulatedoutbreak9[,i]=zoutbreak
}

#sigid10
for(i in 1:nsim){
     set.seed(i)
     yt=round(negbinNoise1(N=N, k=1, k2=1, alpha=3, beta=0,
                           gama1=0.2, gama2=0.1, gama3=0.05,
                           gama4=0.15, phi=1, shift=-200,
                           shift2=-200))
     
     set.seed(i)
     out2=outbreak5(currentday=days5 * nweeks * years, 
                    weeklength=nweeks * days5 * years,
                    wtime=(length(yt)-49 * days5 + 1):length(yt),
                    yi=yt, interval=0.25, k=1, k2=1, alpha=3,
                    beta=0, gama1=0.2, gama2=0.1, gama3=0.05,
                    gama4=0.15, phi=1, shift=-200, shift2=-200,
                    numoutbk=1, peakoutbk=2 * days5, meanlog=0,
                    sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength10[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior10[,i] = outbreaklength10[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     zeros=rep(0,2)
     weekend=seq(days5, days5 * years * nweeks, days5)
     for(s in 1:length(weekend)){
          zt=append(zt, zeros, after=2 * (s - 1) + weekend[s])
          zitot=append(zitot, zeros, after=2 * (s - 1) + weekend[s])
          zoutbreak=append(zoutbreak, zeros, after=2 * (s - 1) + weekend[s])
     }
     
     simulateddata10[,i]=zt
     simulatedtotals10[,i]=zitot
     simulatedoutbreak10[,i]=zoutbreak
}

#sigid11
for(i in 1:nsim){
     set.seed(i)
     yt=round(negbinNoise1(N=N, k=1, k2=1, alpha=5, beta=0,
                           gama1=0.2, gama2=0.1, gama3=0.05,
                           gama4=0.1, phi=1, shift=0, shift2=0))
     
     mu=exp(h1(N=N, k=1, k2=1, alpha=5, beta=0, gama1=0.2,
               gama2=0.1, gama3=0.05, gama4=0.1, shift=0,
               shift2=0))
     
     set.seed(i)
     out2=outbreak5(currentday=days5 * nweeks * years,
                    weeklength=nweeks * days5 * years,
                    wtime=(length(yt)-49 * days5 + 1):length(yt),
                    interval=0.25, yi=yt, k=1, k2=1, alpha=5,
                    beta=0, gama1=0.2, gama2=0.1, gama3=0.05,
                    gama4=0.1, phi=1, shift=0, shift2=0,
                    numoutbk=1, peakoutbk=2 * days5, meanlog=0,
                    sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength11[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior11[,i] = outbreaklength11[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     zeros=rep(0, 2)
     weekend=seq(days5, days5 * years * nweeks, days5)
     for(s in 1:length(weekend)){
          zt=append(zt, zeros, after=2 * (s - 1) + weekend[s])
          zitot=append(zitot, zeros, after=2 * (s - 1) + weekend[s])
          zoutbreak=append(zoutbreak, zeros, after=2 * (s - 1) + weekend[s])
     }
     
     simulateddata11[,i]=zt
     simulatedtotals11[,i]=zitot
     simulatedoutbreak11[,i]=zoutbreak
}


#sigid12
for(i in 1:nsim){
     set.seed(i)
     yt=round(negbinNoise1(N=N, k=2, k2=1, alpha=0.5, beta=0,
                           gama1=0.4, gama2=0, gama3=0.05, 
                           gama4=0.15, phi=1, shift=0, shift2=0))
     
     set.seed(i)
     out2=outbreak5(currentday=days5 * nweeks * years,
                    weeklength=nweeks * days5 * years,
                    wtime=(length(yt)-49 * days5 + 1):length(yt),
                    yi=yt, interval=0.25, k=2, k2=1, alpha=0.5,
                    beta=0, gama1=0.4, gama2=0, gama3=0.05,
                    gama4=0.15, phi=1, shift=0, shift2=0,
                    numoutbk=1, peakoutbk=2 * days5, meanlog=0,
                    sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength12[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior12[,i] = outbreaklength12[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     zeros=rep(0,2)
     weekend=seq(days5, days5 * years * nweeks, days5)
     for(s in 1:length(weekend)){
          zt=append(zt, zeros, after=2 * (s - 1) + weekend[s])
          zitot=append(zitot, zeros, after=2 * (s - 1) + weekend[s])
          zoutbreak=append(zoutbreak, zeros, after=2 * (s - 1) + weekend[s])
     }
     
     simulateddata12[,i]=zt
     simulatedtotals12[,i]=zitot
     simulatedoutbreak12[,i]=zoutbreak
}


#sigid13
for(i in 1:nsim){
     
     set.seed(i)
     yt=round(negbinNoise1(N=N, k=1, k2=1, alpha=9, beta=0,
                           gama1=0.5, gama2=0.2, gama3=0.2,
                           gama4=0.5, phi=1, shift=0, 
                           shift2=0)/100)
     
     set.seed(i)
     out2=outbreak5(currentday=days5 * nweeks * years, 
                    weeklength=nweeks * days5 * years,
                    wtime=(length(yt)-49 * days5 + 1):length(yt),
                    yi=yt, interval=0.25, k=1, k2=1, alpha=9,
                    beta=0, gama1=0.5, gama2=0.2, gama3=0.2,
                    gama4=0.5, phi=1, shift=0, shift2=0,
                    numoutbk=1, peakoutbk=2 * days5, meanlog=0,
                    sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength13[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior13[,i] = outbreaklength13[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     zeros=rep(0,2)
     weekend=seq(days5,days5*years * nweeks, days5)
     for(s in 1:length(weekend)){
          zt=append(zt, zeros, after=2 * (s - 1) + weekend[s])
          zitot=append(zitot, zeros, after=2 * (s - 1) + weekend[s])
          zoutbreak=append(zoutbreak, zeros, after=2 * (s - 1) + weekend[s])
     }
     
     simulateddata13[,i]=zt
     simulatedtotals13[,i]=zitot
     simulatedoutbreak13[,i]=zoutbreak
}

# 7-day week data
years <- days7 <- 7
N     <- nweeks * days7 * years

#sigid1
for(i in 1:nsim){
     set.seed(i)
     yt=round(negbinNoise2(N=N, k=1, k2=2, alpha=6, beta=0,
                           gama1=0.2, gama2=0.2, gama3=0.5, 
                           gama4=0.4, phi=2, shift=29)) 
     
     set.seed(i)
     out2=outbreak7(currentday=N, weeklength=nweeks * days7 * years,
                    wtime=(length(yt)-49 * days7 + 1):length(yt),
                    yi=yt, interval=0.25, k=1, k2=2, alpha=6, beta=0,
                    gama1=0.2, gama2=0.2, gama3=0.5, gama4=0.4, phi=2,
                    shift=29, numoutbk=1, peakoutbk=2 * days7, meanlog=0,
                    sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength1[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior1[,i] = outbreaklength1[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     simulateddata1[,i]=zt
     simulatedtotals1[,i]=zitot
     simulatedoutbreak1[,i]=zoutbreak
}

#sigid3
for(i in 1:nsim){
     set.seed(i)
     yt=round(negbinNoise2(N=N, k=1, k2=2, alpha=0.5, beta=0,
                           gama1=1.5, gama2=1.4, gama3=0.5, 
                           gama4=0.4, phi=1, shift=-167))
     
     set.seed(i)
     out2=outbreak7(currentday=N, weeklength=nweeks * days7 * years,
                    wtime=(length(yt)-49 * 7 + 1):length(yt),
                    yi=yt, interval=0.25, k=1, k2=2, alpha=0.5,
                    beta=0, gama1=1.5, gama2=1.4, gama3=0.5,
                    gama4=0.4, phi=1, shift=-167, numoutbk=1,
                    peakoutbk=2 * days7, meanlog=0, sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength3[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior3[,i] = outbreaklength3[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     simulateddata3[,i]=round(zt)
     simulatedtotals3[,i]=round(zitot)
     simulatedoutbreak3[,i]=round(zoutbreak)
}


#sigid4
for(i in 1:nsim){
     set.seed(i)
     yt=round(negbinNoise2(N=N, k=0, k2=2, alpha=5.5, beta=0, 
                           gama1=0, gama2=0, gama3=0.3, gama4=0.25,
                           phi=1, shift=1))
     
     set.seed(i)
     out2=outbreak7(currentday=N, weeklength=nweeks * days7 * 12,
                    wtime=(length(yt)-49 * 7 + 1):length(yt),
                    yi=yt, interval=0.25, k=0, k2=2, alpha=5.5,
                    beta=0, gama1=0, gama2=0, gama3=0.3, gama4=0.25,
                    phi=1, shift=1, numoutbk=1, peakoutbk=2 * days7,
                    meanlog=0, sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength4[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior4[,i] = outbreaklength4[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     simulateddata4[,i]=zt
     simulatedtotals4[,i]=zitot
     simulatedoutbreak4[,i]=zoutbreak
}

#sigid5
for(i in 1:nsim){
     set.seed(i)
     yt=round(negbinNoise2(N=N, k=0, k2=2, alpha=2, beta=0,
                           gama1=0, gama2=0, gama3=0.3, 
                           gama4=0.25, phi=1, shift=1))
     
     set.seed(i)
     out2=outbreak7(currentday=N, weeklength=nweeks * days7 * years,
                    wtime=(length(yt)-49 * days7 + 1):length(yt),
                    yi=yt, interval=0.25, k=0, k2=2, alpha=2, beta=0,
                    gama1=0, gama2=0, gama3=0.3, gama4=0.25, phi=1,
                    shift=1, numoutbk=1, peakoutbk=2 * days7, meanlog=0,
                    sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength5[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior5[,i] = outbreaklength5[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     simulateddata5[,i]=zt
     simulatedtotals5[,i]=zitot
     simulatedoutbreak5[,i]=zoutbreak
}

#sigid14
for(i in 1:nsim){
     set.seed(i)
     yt=round(negbinNoise2(N=N, k=1, k2=2, alpha=2, beta=0.0005,
                           gama1=0.8, gama2=0.8, gama3=0.8, gama4=0.4,
                           phi=4, shift=57))
     
     set.seed(i)
     out2=outbreak7(currentday=N, weeklength=nweeks * days7 * years,
                    wtime=(length(yt)-49 * days7 + 1):length(yt),
                    yi=yt, interval=0.25, k=1, k2=2, alpha=6,
                    beta=0, gama1=0.2, gama2=0.2, gama3=0.5, 
                    gama4=0.4, phi=2, shift=29, numoutbk=1, 
                    peakoutbk=2 * days7, meanlog=0, sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength14[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior14[,i] = outbreaklength14[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     simulateddata14[,i]=zt
     simulatedtotals14[,i]=zitot
     simulatedoutbreak14[,i]=zoutbreak
}


#sigid15
for(i in 1:nsim){
     set.seed(i)
     yt=round(0.1*(negbinNoise2(N=N, k=4, k2=1, alpha=1.5, beta=0,
                                gama1=0.1, gama2=0.1, gama3=1.8,
                                gama4=0.1, phi=1, shift=-85)+2))
     
     set.seed(i)
     out2=outbreak7(currentday=N, weeklength=nweeks * days7 * years,
                    wtime=(length(yt)-49 * days7+1):length(yt),
                    yi=yt, interval=0.25, k=1, k2=2, alpha=2,
                    beta=0, gama1=0.8, gama2=0.8, gama3=0.8,
                    gama4=0.4, phi=4, shift=57, numoutbk=1,
                    peakoutbk=2 * days7, meanlog=0, sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength15[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior15[,i] = outbreaklength15[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     simulateddata15[,i]=zt
     simulatedtotals15[,i]=zitot
     simulatedoutbreak15[,i]=zoutbreak
}

#sigid16
for(i in 1:nsim){
     set.seed(i)
     yt=round(negbinNoise2(N=N, k=1, k2=2, alpha=3, beta=0,
                           gama1=0.8, gama2=0.6, gama3=0.8,
                           gama4=0.4, phi=4, shift=29))
     
     out1=rep(0, N)
     for(j in 1:years){
          set.seed(j + years * i)
          out=
               outbreak5(currentday=days7 * nweeks * years,
                         weeklength=nweeks * days7 * years,
                         wtime=((210 + (j-1) * days7 * 
                                      nweeks):(230+(j-1) * days7 * nweeks)),
                         yi=yt, interval=0.02, k=1, k2=1, alpha=1,
                         beta=0, gama1=0.1, gama2=2, gama3=0.05, 
                         gama4=0.05, phi=1, shift=-50, shift2=-50,
                         numoutbk=1, peakoutbk=3 * days7 * 150,
                         meanlog=0, sdlog=0.5)
          out1=out1 + out$outbreak
     }
     out1 = round(out1)
     
     set.seed(i)
     out2=outbreak7(currentday=N, weeklength=nweeks * days7 * years,
                    wtime=(length(yt)-49 * days7 + 1):length(yt),
                    yi=yt, interval=0.25, k=1, k2=2, alpha=3, beta=0,
                    gama1=0.8, gama2=0.6, gama3=0.8, gama4=0.4, phi=4,
                    shift=29, numoutbk=1, peakoutbk=2 * days7, meanlog=0,
                    sdlog=0.5)
     # Non-seasonal outbreak
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength16[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior16[,i] = outbreaklength16[,i] / (364 * nyears)
     
     # Seasonal outbreak
     zseasoutbreak=zoutbreak + out1
     
     # Seasonal syndromic data + seasonal outbreak
     zt = yt + out1
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     # Set bank holidays to zero
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     simulateddata16[,i]=zt
     simulatedtotals16[,i]=zitot
     simulatedoutbreak16[,i]=zoutbreak
     simulatedzseasoutbreak16[,i]=zseasoutbreak
}


#sigid17
for(i in 1:nsim){
     set.seed(i)
     yt=round(negbinNoise2(N=N, k=0, k2=2, alpha=6, beta=0,
                           gama1=0, gama2=0, gama3=0.8, gama4=0.4,
                           phi=4, shift=1))
     
     set.seed(i)
     out2=outbreak7(currentday=N, weeklength=nweeks * days7 * 12,
                    wtime=(length(yt)-49 * days7 + 1):length(yt),
                    yi=yt, interval=0.25, k=0, k2=2, alpha=6, 
                    beta=0, gama1=0, gama2=0, gama3=0.8, gama4=0.4,
                    phi=4, shift=1, numoutbk=1, peakoutbk=2 * days7,
                    meanlog=0, sdlog=0.5)
     
     zoutbreak=round(out2$outbreak)
     
     # Non-seasonal outbreak length
     outbreaklength17[,i] = sum(zoutbreak>0)
     
     # Probability of non-seasonal outbreak occurrence
     outbreakprior17[,i] = outbreaklength17[,i] / (364 * nyears)
     
     # Syndromic data + seasonal outbreak
     zt = yt
     
     # Add non-seasonal outbreak on top
     zitot = zt + zoutbreak
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zt[b]=0
               zitot[b]=0
               zoutbreak[b]=0
          } 
     }
     
     simulateddata17[,i]=zt
     simulatedtotals17[,i]=zitot
     simulatedoutbreak17[,i]=zoutbreak
}


# ------------------------------------
# Scale up data after a bank holiday
# ------------------------------------

myfiles <- paste0("simulatedtotals", c(1, 3:17))
files5  <- paste0("simulatedtotals", 6:13)
"%ni%"  <- Negate("%in%")

for(file in myfiles){
     myData = get(file)
     
     for(i in 1:nsim){
          
          for(b in 1:nrow(myData)){
               
               if(bankhols7[b]==1 & file %in% files5){
                    # print(paste(file, "YES"))
                    myData[b + 1, i] = round(1.5 * myData[b + 1, i])
                    assign(file, myData)
               } else { if(bankhols7[b]==1 & file %ni% files5){
                    # print(paste(file, "NO"))
                    myData[b, i] = round(2 * myData[b, i])
                    assign(file, myData)
               }
               }
          }
     }
}


# -----------------------------------
# Define alarm output data frames
# -----------------------------------
# 
# # Create lsits of empty matrices
# lstall <- replicate(narrays*narrays, matrix(nrow=nweeks2 * ndays7, ncol=nsim), 
#                     simplify=FALSE)
# 
# names(lstall) <- paste0("alarm", 1:(narrays*narrays))
# 
# list2env(lstall, .GlobalEnv)

# ------------------------------
# Add each possible outbreak 
# to each possible signal
# ------------------------------
nOutbreaks <- nStreams <- 17

# Combine all time series with all outbreaks
for(i in 1:nStreams){
     sim=get(paste0("simulateddata", i))
     
     for(j in 1:nOutbreaks){
          out=get(paste0("simulatedoutbreak", j))
          tot=sim + out
          
          # Assign to object
          assign(paste0("signal", i, "Out", j), tot)
          
          # Save to file
          saveRDS(tot, 
                  file=file.path(myDir, "output", 
                                 paste0("signal", i, "Out", j, ".rds")))
     }
}

# ---------------------
# Create date series
# ---------------------

myDates   <- seq(ymd('2010-01-01'), ymd('2016-12-30'), by = '1 day')

dropDays <- as.POSIXct(c('2010-12-31','2011-12-31', '2012-12-31',
                         '2013-12-31', '2014-12-31', '2015-12-31', 
                         '2016-02-29,', '2012-02-29'))

"%ni%"  <- Negate("%in%")
myDates <- myDates[myDates %ni% dropDays]

# --------------------------------------------
# Convert daily data to 7-day running totals
# --------------------------------------------

# Rolling sum function
rolling <- function(x){
     rollapplyr(x, width=7, FUN=sum, na.rm=T, fill=NA)
} 

# Smooth data
for(i in 1:nStreams){
     
     for(j in 1:nOutbreaks){
          # Retrieve data
          myData  = get(paste0("signal", i, "Out", j))
          
          # Apply smooth function
          newData = apply(myData, 2, rolling)
          
          # Convert data to STS for CUSUM-NB
          newData = sts(newData, start=c(2010, 1), frequency=364,
                        epoch=as.numeric(as.Date(myDates)), 
                        epochAsDate=TRUE)
          
          # Assign to output object
          assign(paste0("weekSignal", i, "Out", j), newData)
          
          # Save output to file
          saveRDS(newData, 
                  file=file.path(myDir, "output", 
                                 paste0("weekSignal", i, "Out", j, ".rds")))
     }
}


# --------------------------------------
# Run outbreak detection algorithm
# --------------------------------------
# Select range of data to monitor, algorithm and prediction interval
in2016  <- 2206:2548
control <- list(range=in2016, alpha=NULL, mu0=list(S=2, trend=TRUE),
                theta=NULL)

for(i in c(1, 3:nStreams)){
     
     for(j in c(1, 3:nOutbreaks)){
          
          out=matrix(nrow=nweeks2 * ndays7, ncol=nsim)
          
          for(k in seq(nsim)){
               
               sim=get(paste0("weekSignal", i, "Out", j))
               det=algo.glrnb(sts2disProg(sim[, k]), control=control)
               ala=as.numeric(as.vector(unlist(det$alarm)))
               
               # Replace missing values with zero (?)
               ala[is.na(ala)]=0
               
               # Get alarms into matrix
               out[, k]=ala
               
               # Plot alarms
               png(file.path(myDir, "plots", "totals", 
                             paste0("signal", i, "out", j,
                                    "sim", k, ".png")),
                   width=8, height=6, units="in", res=300)
               plot(det, main=paste("signal", i, "out", j, "sim", k), 
                    startyear=2010, legend=NULL)
               dev.off()
          }
          
          # Assign output to object
          assign(paste0("alarm", i, "Out", j), out)
          
          # Save output to file
          saveRDS(out, 
                  file=file.path(myDir, "output", 
                                 paste0("alarm", i, "Out", j, ".rds")))
     }
}


# --------------------------
# Retrieve alarm data
# --------------------------
weeklyAll <- list.files(
     path=file.path(myDir, "output"), 
     pattern="week")


alarmAll <- list.files(
     path=file.path(myDir, "output"), 
     pattern="alarm")

# Read files
for(w in weeklyAll){
     
     # Read
     file  <- readRDS(file.path(myDir, "output", w))
     
     # Assign to file
     assign(gsub(".rds", "", w), file)
}

for(a in alarmAll){
     
     # Read
     file  <- readRDS(file.path(myDir, "output", a))
     
     # Assign to file
     assign(gsub(".rds", "", a), file)
}


# --------------------------
# FPR false positive rate
# --------------------------
days  <- 7
ncomb <- nStreams * nOutbreaks
fpr   <- rep(0, (nStreams-1) * (nOutbreaks-1))

# Prepare output object
combinations <- matrix(0, nrow=nStreams, ncol=nOutbreaks)

for(l in c(1, 3:nStreams)){
     
     for(k in c(1, 3:nOutbreaks)){
          
          combinations[l, k]=paste0("alarm", l, "Out", k)
          
     }
}
combinations <- combinations[-2, ]
combined     <- as.vector(combinations)
combined     <- combined[combined != "0"]

for(k in combined){
     
     # Retrieve alarm data
     myData = get(k)
     
     # Get index for simulated outbreak data
     subStr = stringr::str_sub	(k, start= -4)
     getNum = as.numeric(gsub("\\D", "", subStr)) 
     
     # Retrieve outbreak
     myOutb = get(paste0("simulatedoutbreak", getNum))
     
     # Compute FPR
     nu     = 0
     for(j in 1:nsim){
          
          for(i in (49*7):1){
               nu=(myData[nrow(myData)-i + 1, j]==1 & 
                        myOutb[nrow(myOutb)-i+1,j]==0)+nu
          }
     }
     
     index      = which(combined==k)
     fpr[index] = nu/sum(myOutb==0)  
}

#----------------------------
# POD power of detection2
#----------------------------

pod2 <- rep(0, times=(nStreams-1) * (nOutbreaks-1))

for(k in combined){
     
     # Retrieve alarm data
     myData = get(k)
     
     # Get index for simulated outbreak data
     subStr = stringr::str_sub	(k, start= -4)
     getNum = as.numeric(gsub("\\D", "", subStr)) 
     
     # Retrieve outbreak
     myOutb = get(paste0("simulatedoutbreak", getNum))
     
     # Compute POD2
     mu=0
     for(j in 1:nsim){
          nu=0
          
          for(i in (49*days):1){
               nu=nu+(myData[nrow(myData)-i+1,j]==1 & 
                           myOutb[nrow(myOutb)-i+1,j]>0)
          }
          mu=mu+(nu>0)
     }
     
     index       = which(combined==k)
     pod2[index] = mu/nsim
}


#-----------------------
# Sensitivity
#-----------------------
sensit <- rep(0, times=(nStreams-1) * (nOutbreaks-1))

for(k in combined){
     
     # Retrieve alarm data
     myData = get(k)
     
     # Get index for simulated outbreak data
     subStr = stringr::str_sub	(k, start= -4)
     getNum = as.numeric(gsub("\\D", "", subStr)) 
     
     # Retrieve outbreak
     myOutb = get(paste0("simulatedoutbreak", getNum))
     
     # Compute POD2
     nu=0
     for(j in 1:nsim){
          
          for(i in (49*days):1){
               nu=nu+(myData[nrow(myData)-i+1, j]==1 & 
                           myOutb[nrow(myOutb)-i+1, j]>0)
          }
     }
     
     index         = which(combined==k)
     sensit[index] = nu/(nu+fpr[index]*sum(myOutb==0))
}


#-----------------------
# Specificity
#-----------------------
speci <- rep(0, times=(nStreams-1) * (nOutbreaks-1))

for(k in combined){
     
     # Retrieve alarm data
     myData = get(k)
     
     # Get index for simulated outbreak data
     subStr = stringr::str_sub	(k, start= -4)
     getNum = as.numeric(gsub("\\D", "", subStr)) 
     
     # Retrieve outbreak
     myOutb = get(paste0("simulatedoutbreak", getNum))
     
     # Compute POD2
     nu=0
     for(j in 1:nsim){
          for(i in (49*days):1){
               nu=nu+(myData[nrow(myData)-i+1,j]==0 &
                           myOutb[nrow(myOutb)-i+1,j]==0)
          }
     }
     
     index        = which(combined==k)
     speci[index] = nu/(nu+fpr[index]*sum(myOutb==0))
}


#--------------------
# Timeliness
#--------------------
timel <- rep(0, times=(nStreams-1) * (nOutbreaks-1))

for(k in combined){
     
     # Retrieve alarm data
     myData = get(k)
     
     # Get index for simulated outbreak data
     subStr = stringr::str_sub(k, start= -4)
     getNum = as.numeric(gsub("\\D", "", subStr)) 
     
     # Retrieve outbreak
     myOutb = get(paste0("simulatedoutbreak", getNum))
     
     # Compute POD2
     n=0
     ss=0
     for(j in 1:nsim){
          for(i in (nweeks*days*years):(nweeks*days*(years-1)+3*days+1)){
               test=(myOutb[i,j]>0)
               if(test==TRUE){
                    r2=i
                    break
               }
          }
          for(i in (nweeks*(years-1)*days+3*days+1):(nweeks*years*days)){
               test=(myOutb[i,j]>0)
               if(test==TRUE){
                    r1=i
                    break
               }
          }
          for(i in (49*days):1){
               test=(myData[nrow(myData)-i+1,j]==1 &
                          myOutb[nrow(myOutb)-i+1,j]>0)
               if(test==TRUE){
                    ss=ss+(nrow(myOutb)-i+1-r1)/(r2-r1+1)
                    break
               }
          }
          if(i==1 & test!=TRUE){n=n+1}
     }
     
     index        = which(combined==k)
     timel[index] = (ss+n)/nsim
}


#------------------------
# Summary
#------------------------

alarms <- stringr::str_sub(combined, start=1, end=7)
alarms <- as.numeric(gsub("\\D", "", alarms))

outbrs <- stringr::str_sub(combined, start=-4)
outbrs <- as.numeric(gsub("\\D", "", outbrs)) 

outL <- outP <- rep(0, nOutbreaks)

for(i in 1:nOutbreaks){
     
     outL[i] <- median(get(paste0("outbreaklength", i)))
     outP[i] <- median(get(paste0("outbreakprior", i)))
}

outL <- as.vector(na.omit(outL))
outP <- as.vector(na.omit(outP))

summary <- data.table(id=combined, signal=alarms, outbreak=outbrs, 
                      outLength=rep(outL, each=nOutbreaks-1),
                      outPrior=rep(outP, each=nOutbreaks-1),
                      fpr, pod2, sensit, speci, timel)

summary <- summary %>% dplyr::arrange(signal, outbreak) %>% data.table()

meanStats <-  data.table(fpr, pod2, sensit, speci, timel)
meanStats <- colMeans(meanStats)

medStats <-  data.table(fpr, pod2, sensit, speci, timel)
medStats <- apply(medStats, 2, median)


# -------------------------
# Bayesian stats
# -------------------------
ids <- unique(summary$id)

summary$h1         <- 0
summary$h2         <- 0
summary$postH1     <- 0
summary$postH2     <- 0
summary$priorOddsH1<- 0
summary$priorOddsH2<- 0
summary$postOddsH1 <- 0
summary$postOddsH2 <- 0
summary$bayesFac   <- 0
summary$jeffreys   <- 0
summary$kass       <- 0
summary$postProb   <- 0

for(i in seq_along(ids)){
     
     myData=dplyr::filter(summary, id==ids[i])
     
     # --- Prior odds
     # H1 : There is no ongoing outbreak
     # H2 : There is an ongoing outbreak
     
     summary$h1[i] <- 1 - summary$outPrior[i]
     summary$h2[i] <- summary$outPrior[i]     
     
     # Pior odds (ratio of h1 to h2)
     summary$priorOddsH1[i] <- summary$h1[i] / summary$h2[i]
     
     summary$priorOddsH2[i] <- summary$h2[i] / summary$h1[i]
     
     # Posterior probability
     summary$postH2[i] <- (summary$pod2[i] * summary$h2[i]) / 
          ((summary$pod2[i] * summary$h2[i]) + 
                (summary$fpr[i] * summary$h1[i]))
     
     summary$postH1[i] <- 1 - summary$postH2[i]
     # Alternatively:
     # postH1 <- (summary$fpr[1] * h1) /
     #      + ((summary$pod2[1] * h2) + (summary$fpr[1] * h1))
     
     summary$postOddsH1[i] <- summary$postH1[i] / summary$postH2[i]
     
     summary$postOddsH2[i] <- summary$postH2[i] / summary$postH1[i]
     
     # Bayes factor 
     summary$bayesFac[i] <- summary$postOddsH1[i] / summary$priorOddsH1[i]
     
     # Notice that the posterior odds are therefore computed as the
     # product of the Bayes Factor and the prior odds
     # postOdds <- bayesFac * priorOdds
     
     # Interpreting Bayes Factor
     # Jeffreys (1961):
     # BF{H1:H2}    Evidence against H2
     # 1 - 3        Not worth mentioning
     # 3 - 20       Positive
     # 20 - 150     Strong
     # > 150        Very strong
     
     # Notice the value does not even appear on the scale, 
     # so we have to compute the BF for H2 to H1 (reciprocal)
     
     summary$jeffreys[i] <- 1 / summary$bayesFac[i]
     
     # The result gives evidence against H1. So we conclude there is
     # very strong evidence for the presence of an outbreak
     
     # Kass & Raftery (1995):
     # 2log(BF{H2:H1})   Evidence against H1
     # 0 - 2             Not worth mentioning
     # 2 - 6             Positivie
     # 6 - 10            Strong
     # > 10              Very strong
     
     summary$kass[i] <- 2 * log(summary$jeffreys[i]) 
     # We get the same. The evidence against H1 is very strong.
     
     # Compute posterior probability
     summary$postProb[i] <- summary$postOddsH2[i] / (summary$postOddsH2[i] + 1)
}

summPrior <- summary %>% 
     group_by(outbreak) %>%
     summarise(min=min(pod2),
               q25=quantile(pod2, probs=0.25),
               median=median(pod2),
               q75=quantile(pod2, probs=0.75),
               max=max(pod2)) %>%
     t()

summPrior <- summPrior[-1,]

xtable(summPrior)

summPost <- summary %>% 
     group_by(outbreak) %>%
     summarise(min=min(postProb),
               q25=quantile(postProb, probs=0.25),
               median=median(postProb),
               q75=quantile(postProb, probs=0.75),
               max=max(postProb)) %>%
     t()
     
summPost <- summPost[-1,]
xtable(summPost)

# out3 <- dplyr::filter(summary, outbreak==3)

combinations <- unlist(lapply(2:3, 
                              function(i)combn(c(1, 3:nStreams), i, 
                                               simplify=FALSE)),
                       recursive=FALSE)

postProb <- postOdds <- bFcomb <- matrix(0, nrow=length(combinations), 
                                         ncol=length(1:nOutbreaks))

for(j in c(1, 3:nOutbreaks)){
     
     myOut <- dplyr::filter(summary, outbreak==j)
     
     for(i in seq_along(combinations)){
          
          nItems <- length(as.vector(unlist(combinations[i])))
          
          if(nItems == 2){
               
               # Retrieve bayes factors per signam
               ts1 <- myOut$jeffreys[myOut$signal==combinations[[i]][1]]
               ts2 <- myOut$jeffreys[myOut$signal==combinations[[i]][2]]
               
               # Combine bayes factors
               bFcomb[i, j] <-  ts1 * ts2
               
               # Compute posterior odds
               postOdds[i, j] <- bFcomb[i, j] * myOut$priorOddsH2[1]
               
          }else{
               ts1 <- myOut$jeffreys[myOut$signal==combinations[[i]][1]]
               ts2 <- myOut$jeffreys[myOut$signal==combinations[[i]][2]]
               ts3 <- myOut$jeffreys[myOut$signal==combinations[[i]][3]]
               
               # Combine bayes factors
               bFcomb[i, j] <-  ts1 * ts2 * ts3
               
               # Compute posterior odds
               postOdds[i, j] <- bFcomb[i, j] * myOut$priorOddsH2[1]
               
               # Compute posterior probability
               postProb[i, j] <- postOdds[i, j] / (postOdds[i, j] + 1)
               
          }
     }
}

# Get summary stats 2 time series
dfr <- t(data.table(min=apply(postProb[1:120,], 2, min),
                    q25=apply(postProb[1:120,], 2, quantile, probs=0.25),
                    median=apply(postProb[1:120,], 2, median),
                    q75=apply(postProb[1:120,], 2, quantile, probs=0.75),
                    max=apply(postProb[1:120,], 2, max)))
xtable(dfr)

dfr2 <- t(data.table(min=apply(postProb[121:680,], 2, min),
                     q25=apply(postProb[121:680,], 2, quantile, probs=0.25),
                     median=apply(postProb[121:680,], 2, median),
                     q75=apply(postProb[121:680,], 2, quantile, probs=0.75),
                     max=apply(postProb[121:680,], 2, max)))
xtable(dfr2)
