library(cowplot)

load(file= "~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ggNKIacclatwithweights.Rdata")
load(file= "~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ggPNCacclatwithweights.Rdata")

ggNKIweightedunweighted+theme(legend.position = "none")

ggweightedunweighted<-(ggNKIweightedunweighted+theme(legend.position = "none")+ggPNCweightedunweighted+theme(legend.position = "none"))
ggsave(ggweightedunweighted,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Popweight/popweight.pdf",height=3,width=10)

###supdata###
load(file= "~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ggPNCacclatwithweights.suppdata.Rdata")
PNCwanduwfits<-wanduwfits

load(file= "~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_Behavioral/GeneralScripts/RakingSens/ggNKIacclatwithweights.suppdata.Rdata")
NKIwanduwfits<-wanduwfits

allwanduwfits<-plyr::rbind.fill(NKIwanduwfits,PNCwanduwfits)
write.csv(allwanduwfits,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Data/SupportingData/Sup2.csv")
