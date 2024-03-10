library(tidyverse)

library(lavaan)
dat <- read.csv("finalData.csv")
dat
colnames(dat)


dat %>% filter(RS_Condition!="Original") %>%
  filter(CO_Condition!="Original") ->dat


dat %>% select(CO_I1:CO_I9, RS_I1:RS_I10, LOT_I1:LOT_I10, GPA) %>%  na.omit() ->dat
nrow(dat)

dat$GPA <- round((dat$GPA/100)*4.33,2)
head(dat)

write.csv(dat, "zhang2019data.csv")
