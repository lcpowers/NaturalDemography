## Compile 2023 rock wall plot iButton data ##

library(tidyverse)
library(vroom)
library(lubridate)
library(reshape2)
setwd("~/Desktop/Research/Fieldwork/soil_temp_data/2023/nr_walls")
rm(list=ls())

# plots = expand.grid("p",1:10,c("north","south","east","west","control")) %>% do.call(paste0,.) %>% sort()
plots = paste0("p",1:10)
temps = NULL
plot_order = paste0("p",1:10) %>% as.factor()

for(plot in plots){
  
  fns = list.files(paste0("./",plot),full.names = T)
  
  df = vroom(fns,skip = 14,id="path") %>% 
    mutate(bsnm = substr(basename(path),1,nchar(basename(path))-4))
  
  df$plot = plot
  
  if(plot=="p10"){
    
    df$aspect = str_sub(df$bsnm,4,nchar(df$bsnm))
  
  }else{
    
    df$aspect = str_sub(df$bsnm,3,nchar(df$bsnm))
    
      }
  
  df_out = select(df,-c(bsnm,path)) %>% 
    separate(`Date/Time`,into = c("date","time"),sep = " ",remove = F) %>% 
    separate(date, into = c("month","day","year")) %>% 
    # fix date messes
    mutate(year=ifelse(nchar(year)==2,paste0("20",year),year),
           date = as_date(paste(year,month,day,sep="-"))) %>% 
    # fix time messes
    separate(time,into=c("h","m","s"),remove=F) %>% 
    mutate(s = ifelse(is.na(s),"00",s),
           h = ifelse(nchar(h)==1,paste0("0",h),h),
           time = paste(h,m,s,sep=":")) %>% 
    select(-c(`Date/Time`,h,m,s,month,day,year,Unit)) %>% 
    mutate(dt = paste(date,time)) %>% 
    rename(temp = Value)
    
    # mutate(date2 = as.Date(date,format="%m/%d/%Y"))
    # mutate(date = as_date(`Date/Time`))
  
  temps = rbind(temps,df_out)
  
  }

temps$plot = factor(temps$plot,levels = plot_order)

daily_temps = temps %>% 
  filter(!str_detect(aspect,"replace")) %>% 
  group_by(date,plot,aspect) %>% 
  summarise(meanT = mean(temp),
            lowT = min(temp),
            highT = max(temp))


ggplot(daily_temps,aes(x=date,y=lowT,color=aspect))+
  geom_point()+
  facet_wrap(~plot)

write_csv(temps,"all_iButtonsTs23.csv")

annual_summary = temps %>%
  filter(!str_detect(aspect,"replace")) %>% 
  group_by(plot,aspect) %>% 
  summarise(meanT = mean(temp),
            lowT = min(temp),
            highT = max(temp)) %>% 
  melt(id.vars = c("plot","aspect"))

ggplot(annual_summary,aes(x=aspect,y=value,fill=aspect))+
  geom_col()+
  facet_grid(plot~variable,scales="free")

winter_summary = temps %>%
  filter(!str_detect(aspect,"replace")) %>% 
  filter(date<"2023-04-01"&date>"2022-11-01") %>% 
  group_by(plot,aspect) %>% 
  summarise(meanT = mean(temp),
            lowT = min(temp),
            highT = max(temp),
            varT = var(temp)) %>% 
  melt(id.vars = c("plot","aspect"))

ggplot(winter_summary,aes(x=aspect,y=value,fill=aspect))+
  geom_col()+
  facet_grid(plot~variable,scales="free")

total_summary = temps %>% 
  filter(!str_detect(aspect,"replace")) %>% 
  group_by(aspect) %>% 
  summarise(meanT = mean(temp),
            lowT = min(temp),
            highT = max(temp),
            varT = var(temp))
