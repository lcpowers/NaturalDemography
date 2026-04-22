library(tidyverse)


setwd("~/Desktop/Doaklab/Research/FieldWork/summer2022tempdata/")

fns = list.files(path = ".",full.names = T,pattern = ".csv")

gp = fns[str_detect(fns,"gp")]
jp = fns[str_detect(fns,"jp")]
niw = fns[str_detect(fns,"niw")]

##### GP ######
gp_temps = NULL
for(fn in gp){
  
  tmp = read_csv(fn,skip = 1) %>% 
    select(contains(c("Date","Temp")))
  colnames(tmp) = c("datetime","t")
  tmp$date = as.Date(tmp$datetime,format = "%m/%d/%y")
  tmp$time = format(as.POSIXct(tmp$datetime,format = "%m/%d/%y %I:%M:%S %p"), format = "%H:%M:%S")
  tmp = select(tmp,date,time,t) %>% 
    mutate(site="gp",
           plot=substr(basename(fn),4,14),
           aspect=substr(plot,7,12),
           elev = substr(plot,1,5))
  
  gp_temps = rbind(gp_temps,tmp)
  
}

gp_daily = gp_temps %>% 
  group_by(date,site,plot,elev,aspect) %>% 
  summarize(t.mean = mean(t,na.rm=T)) %>% 
  filter(date>"2022-07-19")

ggplot(gp_daily,aes(x=date,y=t.mean,color=plot))+
  geom_point()+
  geom_line()+
  facet_wrap(~elev)
######

##### JP ######
jp_temps = NULL
for(fn in jp){
  
  tmp = read_csv(fn,skip = 1) %>% 
    select(contains(c("Date","Temp")))
  colnames(tmp) = c("datetime","t")
  tmp$date = as.Date(tmp$datetime,format = "%m/%d/%y")
  tmp$time = format(as.POSIXct(tmp$datetime,format = "%m/%d/%y %I:%M:%S %p"), format = "%H:%M:%S")
  tmp = select(tmp,date,time,t) %>% 
    mutate(site="jp",
           plot=substr(basename(fn),4,13),
           aspect=substr(plot,6,12),
           east.west = substr(plot,1,4))
  
  jp_temps = rbind(jp_temps,tmp)
  
}

jp_daily = jp_temps %>% 
  group_by(date,site,plot,east.west,aspect) %>% 
  summarize(t.mean = mean(t,na.rm=T)) %>% 
  filter(date>"2022-07-21")

ggplot(jp_daily,aes(x=date,y=t.mean,color=plot))+
  geom_point()+
  geom_line()+
  facet_wrap(~east.west)
######

##### Niwot #####
niw_temps = NULL
for(fn in niw){
  
  tmp = read_csv(fn,skip = 1) %>% 
    select(contains(c("Date","Temp")))
  colnames(tmp) = c("datetime","t")
  tmp$date = as.Date(tmp$datetime,format = "%m/%d/%y")
  tmp$time = format(as.POSIXct(tmp$datetime,format = "%m/%d/%y %I:%M:%S %p"), format = "%H:%M:%S")
  tmp = select(tmp,date,time,t) %>% 
    mutate(site="niw",
           plot=substr(basename(fn),7,17),
           aspect=substr(plot,7,12),
           elev = substr(plot,1,5))
  
  niw_temps = rbind(niw_temps,tmp)
  
}

niw_daily = niw_temps %>% 
  group_by(date,site,plot,elev,aspect) %>% 
  summarize(t.mean = mean(t,na.rm=T))# %>% 
  # filter(date>"2022-07-21")

ggplot(niw_daily,aes(x=date,y=t.mean,color=plot))+
  geom_point()+
  geom_line()

