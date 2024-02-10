library(lubridate)
library(tidyverse)
#import excel datasets
alpha <- readxl::read_excel("/Users/eyad/Desktop/LIDo LSHTM/LSHTM 1st Rotation/HOCO Alpha Dates.xlsx")
delta <- readxl::read_excel("/Users/eyad/Desktop/LIDo LSHTM/LSHTM 1st Rotation/HOCO Delta Dates.xlsx")

delta_format <- ymd(delta$Date) #set date format as ymd
delta$Date <- decimal_date(delta_format) #change to numeric according to format

alpha_format <- ymd(alpha$Date)
alpha$Date <- decimal_date(alpha$Date)
