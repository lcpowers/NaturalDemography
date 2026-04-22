setwd("/Volumes/Macintosh HD/Users/clairepowers/Desktop/Research/Fieldwork/soil_temp_data/2023")
library(tidyverse)
library(lubridate)
rm(list = ls())

gp.fns = list.files("./gp",pattern = ".csv",full.names = T)
nr.fns = list.files("./nr",pattern = ".csv",full.names = T)
jp.fns = list.files("./jp",pattern = ".csv",full.names = T)

##### Guanella Pass pendant data 2022-23 #####

gp.temps = NULL

for(fn in gp.fns){
  
  tmp.bn = substr(basename(fn),4,nchar(basename(fn))-4)
  tmp.asp = ifelse(str_detect(tmp.bn,"orth"),"N",
                   ifelse(str_detect(tmp.bn,"outh"),"S","T"))
  tmp.elev = ifelse(str_detect(tmp.bn,"pper"),"upper","lower")
  
  tmp.df = read_csv(fn,skip = 1) %>% 
    select(2,3)
  colnames(tmp.df) = c("dt","T")
  tmp.df = tmp.df %>% 
    mutate(dt = mdy_hms(dt)) %>% 
    separate(dt,sep = " ", into = c("date","time"),remove=F) %>% 
    mutate(site = "gp",aspect = tmp.asp,elev=tmp.elev) %>% 
    filter(date>"2022-08-21")
  
  gp.temps = rbind(gp.temps,tmp.df)   
  rm(tmp.bn,tmp.asp,tmp.df)
  
}

#####

##### Jones Pass pendant data #####

jp.temps = NULL

for(fn in jp.fns){
  
  tmp.bn = substr(basename(fn),4,nchar(basename(fn))-4)
  tmp.asp = ifelse(str_detect(tmp.bn,"orth"),"N",
                   ifelse(str_detect(tmp.bn,"outh"),"S","T"))
  tmp.elev = ifelse(str_detect(tmp.bn,"est"),"west","east")
  
  tmp.df = read_csv(fn,skip = 1) %>% 
    select(2,3)
  colnames(tmp.df) = c("dt","T")
  tmp.df = tmp.df %>% 
    mutate(dt = mdy_hms(dt)) %>% 
    separate(dt,sep = " ", into = c("date","time"),remove=F) %>% 
    mutate(site = "jp",aspect = tmp.asp,elev=tmp.elev) %>% 
    filter(date>"2022-08-21")
  
  jp.temps = rbind(jp.temps,tmp.df)   
  rm(tmp.bn,tmp.asp,tmp.df)
  
}

#####

##### Niwot Ridge pendant data #####
nr.temps = NULL

for(fn in nr.fns){
  
  tmp.bn = substr(basename(fn),4,nchar(basename(fn))-4)
  tmp.asp = ifelse(str_detect(tmp.bn,"orth"),"N",
                   ifelse(str_detect(tmp.bn,"outh"),"S","T"))
  tmp.elev = ifelse(str_detect(tmp.bn,"pper"),"upper","east")
  
  tmp.df = read_csv(fn,skip = 1) %>% 
    select(2,3)
  colnames(tmp.df) = c("dt","T")
  tmp.df = tmp.df %>% 
    mutate(dt = mdy_hms(dt)) %>% 
    separate(dt,sep = " ", into = c("date","time"),remove=F) %>% 
    mutate(site = "nr",aspect = tmp.asp,elev=tmp.elev) %>% 
    filter(date>"2022-08-21")
  
  nr.temps = rbind(nr.temps,tmp.df)   
  rm(tmp.bn,tmp.asp,tmp.df)
  
}

#####

temps = rbind(gp.temps,jp.temps,nr.temps)
write_csv(temps,"~/Desktop/Research/Fieldwork/soil_temp_data/rawT_2223.csv")

daily.temps = temps %>% 
  group_by(date,site,aspect,elev) %>% 
  summarise(meanT = mean(T,na.rm=T),
            varT = var(T),
            minT=min(T),
            maxT=max(T))

write_csv(daily.temps,"~/Desktop/Research/Fieldwork/soil_temp_data/dailyT_2223.csv")

