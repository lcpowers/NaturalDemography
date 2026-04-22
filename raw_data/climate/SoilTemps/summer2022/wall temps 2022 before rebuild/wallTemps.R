## Looking at preliminary temp data from decrontructed walls ##
## July 2022

library(readxl)
library(tidyverse)

rm(list=ls())
setwd("~/Desktop/walltemps")

fns = list.files(".",pattern = ".xlsx",full.names = T)

temps = NULL
daily = NULL
for(f in fns){
  
  plot = basename(f) %>% substr(1,2) 
  df = read_xlsx(f, skip = 22) %>%
    tail(nrow(.)-2)
  
  colnames(df) = c("obs","time","temp")
  
    df = df %>% 
      mutate(plot = plot,
             plotnum = substr(plot,2,2),
           aspect=substr(plot,1,1),
           time = as.POSIXct(time),
           temp = as.numeric(temp),
           date = as.Date(time)) %>% 
      filter(date>as.Date("2022-07-08"),
             date<as.Date("2022-07-25"))
  temps = rbind(temps,df)
  
  df2 = df %>% 
    group_by(date,plot,plotnum,aspect) %>% 
    summarise(temp = mean(temp))
  
  daily = rbind(daily,df2)
  rm(df,df2,plot)
  
}

summary(temps)

ggplot(temps,aes(x=time,y=temp,color=plot))+
  #geom_point()+
  geom_line()+
  facet_wrap(plotnum~aspect)

ggplot(daily,aes(x=date,y=temp,color=plot))+
  #geom_point()+
  geom_line()+
  facet_wrap(plot~aspect)

