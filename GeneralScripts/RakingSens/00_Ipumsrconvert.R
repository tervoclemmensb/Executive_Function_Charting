if (!require("ipumsr")) stop("Reading IPUMS data into R requires the ipumsr package. It can be installed using the following command: install.packages('ipumsr')")
setwd("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/")
ddi <- read_ipums_ddi("usa_00005.xml")
data <- read_ipums_micro(ddi)
dataages<-data[data$AGE>=8 & data$AGE<=35,]
save(dataages,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ACSdataages.Rdata")
