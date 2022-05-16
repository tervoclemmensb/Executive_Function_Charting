####one off clean function for saving correlation matrices
##specific formatting expectations to cog growth charting paper for now
corsaveclean<-function(cors,groupslong){
  lower<-round(cors,3)
  lower[upper.tri(cors, diag=TRUE)]<-""
  lower<-as.data.frame(lower)
  lower$outcome<-row.names(lower)
  lower$seq<-1:nrow(lower)
  lowermerge<-merge(lower,groupslong[,c("outcome","groupbytype")],by=c("outcome"))
  lowermerge<-lowermerge[order(lowermerge$seq),]
  lowermerge<-lowermerge[,c("groupbytype",names(lowermerge)[!grepl("groupbytype",names(lowermerge))])]
  lowermerge[,c("outcome","seq")]<-NULL
  return(lowermerge)
}

####between####
# source("~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/GeneralScripts/corsave.R")
# corssaveclean<-corsaveclean(cors,groupslong)
# write.table(lowermerge,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/Lunabaselinecorrs.csv",row.names=FALSE,col.names=FALSE,sep=",",quote=FALSE)

# ####within########
# wgcorsave<-data.frame(stabsycog$rwg)
# row.names(wgcorsave)<-gsub(".wg","",row.names(wgcorsave))
# names(wgcorsave)<-gsub(".wg","",names(wgcorsave))
# corssavecleanwithin<-corsaveclean(wgcorsave,groupslong)
# write.table(corssavecleanwithin,file="~/Library/Mobile\ Documents/com~apple~CloudDocs/Projects/R03_behavioral/Figures/Factorfigs/Lunalongitudinalcorrs.csv",row.names=FALSE,col.names=FALSE,sep=",",quote=FALSE)
