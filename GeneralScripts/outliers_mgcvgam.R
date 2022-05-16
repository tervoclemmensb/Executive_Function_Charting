windzthrehold<-function(x,sdcutoff=4){
  x<-scale(x)
  x[abs(x)>sdcutoff]<-NA
  return(x)
}
mgcvgamleveragedetect<-function(outcome,covar,idvar=NULL,sdcutoff=4,df,unioutfirst=FALSE,unicutoff=4){
  #this function does leverage statitiscs based on gam residuals
  #optionally does univariate outlier detection first
  if (unioutfirst){
    outcome[abs(scale(outcome))>unicutoff]<-NA
  }
  
  
  outcome<-unlist(outcome)
 # print(outcome)
  covar<-unlist(covar)
  model<-mgcv::gam(outcome~s(covar))
  m<-model
  if (!is.null(idvar)){
    df$id<-df[,idvar]
    df$id<-as.factor(df$id)
    df$outcome<-outcome
    df$covar<-covar
    ctrl <- list(maxIter = 10000, msMaxIter = 10000)
    #model<-mgcv::gamm(outcome~s(covar),random=list(id=~covar),data=df,control = ctrl)
    model <- tryCatch(
      {
        mgcv::gamm(outcome~s(covar),
                   random=list(id=~covar),data=df,control = ctrl)
      },
      error = function(e){
        print("interceptsonly")
        mgcv::gamm(outcome~s(covar),
                   random=list(id=~1),data=df,control = ctrl)
      }
    )
    m<-model$gam
  }
  zresid<-scale(m$residuals)
  outputoutcome<-outcome
  outputoutcome[abs(as.numeric(zresid))>sdcutoff]<-NA
  return(outputoutcome)
}
# distributeNAsfrommeans<-function(meandf,longdf){
#   substocensor<-meandf[!complete.cases(meandf),]
#   for (i in 1:nrow(substocensor)){
#     subi<-substocensor[i,]
#     NAcols<-names(subi)[which(is.na(subi))]
#     longdf[longdf$id==subi$id,NAcols]<-NA
#     print(subi$id)
#     print(length(which(longdf$id==subi$id)))
#   }
#   return(longdf)
# }

distributeNAsforoutliermeans<-function(meandf,longdf,usevars,idvar,sdcutoff=4){
  ##expects data in usevars to be absolute value scaled
  meandf<-as.data.frame(meandf)
  longdf<-as.data.frame(longdf)
  meandf$id<-as.character(meandf[,idvar])
  longdf$id<-as.character(longdf[,idvar])
  meandfscale<-meandf
  meandfscale[,usevars]<-abs(scale(meandfscale[,usevars]))
  ######
  for (i in 1:nrow(meandfscale)){
    subi<-meandfscale[i,]
    outcols<-names(subi[,usevars])[which(subi[,usevars]>sdcutoff)]
    if (length(outcols)!=0){
    print("censoring based on mean")
    longdf[longdf$id==subi$id,outcols]<-NA
    print(subi$id)
    print(length(which(longdf$id==subi$id)))
    print(outcols)
    }
  }
  return(longdf)
}