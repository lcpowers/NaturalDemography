setwd("~/Library/CloudStorage/OneDrive-UCB-O365/Research/Dissertation/Fieldwork/SoilTemps")
library(tidyverse)
library(lubridate)
library(vroom)
rm(list = ls())

gp.fns = list.files("guanella",pattern = ".csv",full.names = T,recursive=T)
# nr.fns = list.files("./nr",pattern = ".csv",full.names = T)
# jp.fns = list.files("./jp",pattern = ".csv",full.names = T)

gp.temps.list <- list()

for(i in 1:length(gp.fns)){
  
  name.parts = strsplit(basename(gp.fns[i]),"_") %>% unlist()
  tmp.plot = name.parts[2]
  tmp.pendID = str_remove(name.parts[3],".csv")
  
  tmp.df = read_csv(gp.fns[i],skip = 1) %>% 
    select(dateTime=2,temp=3)

  tmp.df = tmp.df %>% 
    mutate(dt = mdy_hms(dateTime)) %>% 
    separate(dt,sep = " ", into = c("date","time"),remove=F) %>% 
    mutate(site = "gp",plot = tmp.plot,pendID = tmp.pendID) %>% 
    select(site,plot,pendID,date,time,temp) 
  
  gp.temps.list[[i]] <- tmp.df

}

gp.temps <- data.table::rbindlist(gp.temps.list)

gp.daily.temps <- gp.temps %>% 
  group_by(site,plot,pendID,date) %>% 
  summarise(dailyT = mean(temp,na.rm=T)) %>% 
  group_by(site,plot,date) %>% 
  summarise(dailyT = mean(dailyT,na.rm=T)) %>% 
  mutate(year = year(date),
         yday = yday(date),
         mday = format(as.Date(date,format="%Y-%m-%d"), format = "%m-%d"))

ggplot(gp.daily.temps,aes(x=mday,y=dailyT,color=factor(year)))+
  geom_point()+
  facet_wrap(~plot,ncol=2)+
  #geom_smooth()+
  theme_classic()+
  theme(axis.text.x=element_blank())

