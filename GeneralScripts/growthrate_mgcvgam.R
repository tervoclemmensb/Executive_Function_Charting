####convienence functions
ttor<-function(ts,nobs){
  rfromt<-sqrt((t^2)/(t^2+nobs))
}
find_covars_gam <- function(fml, ...) {
  ind <- as.character(fml)[3]
  # s(x1) + x2 +s(x3,"re") -> x1, x2
  vars <- unlist(strsplit(ind, "\\+"))  # formula split by +
  vars <- gsub(" ", "", vars) # remove spaces
  vars <- gsub("\\w+\\((.*?)\\)", "\\1", vars) # remove surrounding s()
  # remove random effect
  no_re <- grep("re[\"']", vars, value=T, invert=T)
  no_re <- gsub(",.*", "", no_re) # nothing after comma
  # remove anything else (likely passed in agevar)
  if (length(list(...)) != 0L) {
    no_re <- no_re[ ! no_re  %in% c(...)]
  }
  return(no_re)
}
too_small <- function(x) abs(x) < 10^-15
clip_on_sig <- function(ci){
  # if confidence interval includes zero
  # signs of x and y will be different, -x * +y  < 0
  # or if both high and low are extremly close to zero
  not_sig <- ci$ci_low * ci$ci_high < 0 |
    (too_small(ci$ci_low) & too_small(ci$ci_high))
  ci$mean_dff_clip <- ci$mean_dff
  ci$mean_dff_clip[not_sig] <- 0
  return(ci)
}
clip_on_siggratia <- function(ci){
  # if confidence interval includes zero
  # signs of x and y will be different, -x * +y  < 0
  # or if both high and low are extremly close to zero
  not_sig <- ci$lower * ci$upper < 0 |
    (too_small(ci$lower) & too_small(ci$upper))
  ci$sig <- 1
  ci$sig[not_sig] <- 0
  return(ci)
}
####gam residual age######
mgcvresidualizecovar<-function(outcome,covar,idvar=NULL,df){
  #this function residualizes covar(e.g., age) via mgcv gam and returns
  returnresids<-rep(NA,length(outcome))
  outcome<-unlist(outcome)
  # print(outcome)
  covar<-unlist(covar)
  outcomecovar<-data.frame(cbind(outcome,covar))
  outcoemvarcomp<-outcomecovar[complete.cases(outcomecovar),]
  completeindices<-which(complete.cases(outcomecovar))
  incompleteindices<-which(!complete.cases(outcomecovar))
  cxmodel<-mgcv::gam(outcoemvarcomp$outcome~s(outcoemvarcomp$covar))
  m<-cxmodel
  returnresids[completeindices]<-cxmodel$residuals
  return(returnresids)
}

########main functions
oldmgcvgampreddata<-function(model,predvar,idvar=NULL,interval_inc=.1){
  if (class(model)[1]=="gamm"){
    model<-model$gam
  }
  modeldata<-data.frame(ydata=model$y, predvar=model$model[, predvar])
  if (!identical(find_covars_gam(model$formula, predvar),character(0))){
    for (cv in find_covars_gam(model$formula,predvar)){
    x <- model$model[, cv]
    if (is.character(x) || is.factor(x) ){
      warning("gam w/factor covar, setting all sim to the first!")
      y <- x[1]
      # TODO: maybe pracma::Mode ?
    } else {
      y <- mean(x, na.rm=T)
    }
    #pp[, cv] <- y
  }
  
  }else{
    y<-"no cov"
  }
  preddata<-data.frame(var=seq(min(modeldata$predvar),max(modeldata$predvar),by=interval_inc),covar=y)
  names(preddata)<-c(predvar,find_covars_gam(model$formula, predvar))
  if (identical(find_covars_gam(model$formula, predvar),character(0))){
    names(preddata)<-c(predvar,"nullcovar") 
  }
  
  yhats <- predict(model,preddata,se.fit=TRUE)
  preddata<-cbind(preddata,yhats$fit,yhats$se.fit)
  names(preddata)<-c(predvar,find_covars_gam(model$formula, predvar),"fit","se")
  if (identical(find_covars_gam(model$formula, predvar),character(0))){
    names(preddata)<-c(predvar,"nullcovar","fit","se") 
  }
  preddata$CI<-2*preddata$se
  return(preddata)
}
mgcvgampreddata<-function(model,predvar,idvar=NULL,interval_inc=.1,varycovarsname=NULL,varycovarlevels=NULL){
  if (class(model)[1]=="gamm"){
    model<-model$gam
  }
  modeldata<-data.frame(ydata=model$y, predvar=model$model[, predvar])
  preddata<-data.frame(var=seq(min(modeldata$predvar),max(modeldata$predvar),by=interval_inc))
  names(preddata)[1]<-predvar
  if (!identical(find_covars_gam(model$formula, predvar),character(0))){
   # if (!(length(find_covars_gam(model$formula, predvar))==1 && find_covars_gam(model$formula, predvar)[1]==varycovarsname)){
    for (cv in find_covars_gam(model$formula,predvar)){
      x <- model$model[, cv]
      if (is.character(x) || is.factor(x)){
        warning("gam w/character or factor covar, setting all sim to the first obs for character/first level if factor")
        y <- x[1]
        if (class(x)=="factor"){y<-levels(x)[1]}
        print(sprintf("covar % set to level %s",cv,y))
      } else {
        y <- mean(x, na.rm=T)
      }
      preddata[, cv] <- y
    }
    #}
  }else if(is.null(varycovarsname)){
    preddata$cov<-"no cov"
  }
  
  #if (!(length(find_covars_gam(model$formula, predvar))==1 && find_covars_gam(model$formula, predvar)[1]==varycovarsname)){
  names(preddata)<-c(predvar,find_covars_gam(model$formula, predvar))
  #}
  if (identical(find_covars_gam(model$formula, predvar),character(0)) && is.null(varycovarsname)){
    names(preddata)<-c(predvar,"nullcovar") 
  }
  
  if (!is.null(varycovarsname)){
    require(reshape)
    orignameorder<-names(preddata)
    preddata[,varycovarsname]<-NULL
    preddata<-reshape::expand.grid.df(preddata,data.frame(varycovar=varycovarlevels))
    names(preddata)[names(preddata)=="varycovar"]<-varycovarsname
    #preddata<-preddata[,orignameorder]
  }
  
  yhats <- predict(model,preddata,se.fit=TRUE)
  preddata<-cbind(preddata,yhats$fit,yhats$se.fit)
  names(preddata)<-c(predvar,find_covars_gam(model$formula, predvar),"fit","se")
  if (identical(find_covars_gam(model$formula, predvar),character(0)) && is.null(varycovarsname)){
    names(preddata)<-c(predvar,"nullcovar","fit","se")
  }else if(identical(find_covars_gam(model$formula, predvar),character(0)) && !is.null(varycovarsname)){
    names(preddata)<-c(predvar,varycovarsname,"fit","se")
  }
  preddata$CI<-2*preddata$se
  return(preddata)
}

smooth_diff <- function(model, newdata, f1, f2, var, alpha = 0.05,
                        unconditional = FALSE) {
  xp <- predict(model, newdata = newdata, type = 'lpmatrix')
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other lochs
  X[, ! (c1 | c2)] <- 0
  ## zero out the parametric cols
  X[, !grepl('^s\\(', colnames(xp))] <- 0
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  data.frame(pair = paste(f1, f2, sep = '-'),
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr)
}

mgcvscalefits<-function(df,outcomevars,predvars,idvar=NULL,mformula,interval_inc=.1,scale=TRUE){
  pairs<-as.data.frame(expand.grid(outcomevars,predvars))
  names(pairs)<-c("outcome","pred")
  allscalefits<-lapply(1:nrow(pairs),function(p){
    if (scale){
      df$outcome<-scale(unlist(df[,as.character(pairs$outcome[p])])) 
    }else{
      df$outcome<-unlist(df[,as.character(pairs$outcome[p])]) 
    }
    df$pred<-unlist(df[,as.character(pairs$pred[p])])
    print(pairs[p,])
    model<-mgcv::gam(mformula,data=df)
    if (!is.null(idvar)){
      df$id<-df[,idvar]
      ctrl <- list(maxIter = 10000, msMaxIter = 10000)
      model <- tryCatch(
        {
          model<-mgcv::gamm(mformula,random=list(id=~1+pred),data=df,control = ctrl)
        },
        error = function(e){
          print("interceptsonly")
          model<-mgcv::gamm(mformula,random=list(id=~1),data=df,control = ctrl)
        }
      )
    }
    modelpred<-mgcvgampreddata(model,predvar="pred",idvar=idvar,interval_inc = interval_inc)
    modelpred$fitscale<-scales::rescale(modelpred$fit,c(0,1))
    modelpred$fitz<-base::scale(modelpred$fit)[,1]
    modelpred$outcome<-as.character(pairs$outcome[p])
    
    if (!is.null(idvar)){
      modelpred$samplesize<-length(unique(df[,idvar]))
    }else{
      modelpred$samplesize<-nrow(df)
    }
    
    if (class(model)[1]=="gamm"){
      nps<-length(summary(model$gam)$s.pv)
      if (nps==1){
        modelpred$modelp<-summary(model$gam)$s.pv 
        modelpred$anovamodelp<-anova.gam(model$gam)$s.table[1,4]
        modelpred$Fstat<-summary(model$gam)$s.table[1,3]
        modelpred$edf<-summary(model$gam)$edf
      } else {  
        modelpred$modelp<-summary(model$gam)$s.pv[which(grepl("pred",row.names(summary(model$gam)$s.table)))]
        modelpred$anovamodelp<-anova.gam(model$gam)$s.table[which(grepl("pred",row.names(anova.gam(model$gam)$s.table))),4]
        modelpred$Fstat<-anova.gam(model$gam)$s.table[which(grepl("pred",row.names(anova.gam(model$gam)$s.table))),3]
        modelpred$edf<-summary(model$gam)$edf[which(grepl("pred",row.names(summary(model$gam)$s.table)))]
      }
    } else {
      nps<-length(summary(model)$s.pv)
      if (nps==1){
        modelpred$modelp<-summary(model)$s.pv
        modelpred$anovamodelp<-anova.gam(model)$s.table[1,4]
        modelpred$Fstat<-summary(model)$s.table[1,3]
        modelpred$edf<-summary(model)$edf
      } else {
        modelpred$modelp<-summary(model)$s.pv[which(grepl("pred",row.names(summary(model)$s.table)))]
        modelpred$anovamodelp<-anova.gam(model)$s.table[which(grepl("pred",row.names(anova.gam(model)$s.table))),4]
        modelpred$Fstat<-anova.gam(model)$s.table[which(grepl("pred",row.names(anova.gam(model)$s.table))),3]
        modelpred$edf<-summary(model)$edf[which(grepl("pred",row.names(summary(model)$s.table)))]
      }}
    return(modelpred)
  })
  allscalefits<-dplyr::bind_rows(allscalefits)
  return(allscalefits) 
}
###bootstrap and jacknife wrappers#####
jackknife_mgcvgamgrowthrate<-function(modelformula,df,jacknifevar,idvar=NULL,mc.cores){
  mclapply(unique(df[,jacknifevar]),mc.cores=mc.cores,function(j){
    str(j)
    jdf<-df[df[,jacknifevar]!=j,]
    jm<-gamm4(modelformula,data=jdf,random=~(1|idvect))
    jf_ci<-gamm4_growthrate_maturationpoint(jm,agevar="pred",idvar="idvect")
    jf_ci$jf_id<-as.character(j)
    return(jf_ci)
  })
}
####proportion of mat point#####
###sim CI proportion of mat point####
proportionofmaxmaturation_mgcvgamsimCI<-function(m,predvar,idvar=NULL,dervproportions=NULL,totalproportions=NULL,interval_inc=.1,scalederivtointerval=TRUE,scaletomaximum=FALSE,n.iterations=10000,ciqnt=c(0.025,0.975)){
  modelpred<-mgcvgampreddata(m,predvar="pred",idvar=idvar,interval_inc=interval_inc)
  signflip=FALSE
  if(sign(modelpred$fit[modelpred$pred==max(modelpred$pred)]-modelpred$fit[modelpred$pred==min(modelpred$pred)])==-1){
    modelpred$fit<-modelpred$fit*-1
    signflip=TRUE ####flip sign so decreases with age (e.g., latency) can use the same % maximum scaling code
  }
  modelpred$fitscale<-rescale(modelpred$fit,c(0,1))
  modelpred$diff<-c(NA,diff(modelpred$fitscale))
  gamsims<-sim_from_mgcvgam(m, predvar, idvar=idvar, n.iterations = n.iterations,interval_inc=interval_inc,scaletomaximum=scaletomaximum)
  derivs<-gamsims$diffs

  derivmatpointsdf<-NULL
  if (!is.null(dervproportions)){
    names(derivs)<-seq(1,length(derivs))
    diffsimsmat <- data.frame(dplyr::bind_rows(derivs))
  #####
  if (scalederivtointerval){scaledderivproportions<-(interval_inc/1)*dervproportions}else{ ###divde by interval (e.g., )
    scaledderivproportions<-dervproportions
  }
  if (scaletomaximum){
    maximumderivproportion<-1/nrow(modelpred)
  scaledderivproportions<-c(scaledderivproportions,maximumderivproportion)
  }
  #########deriv mat points
  derivmatpoints<-lapply(scaledderivproportions,function(p){ ###apply thresholds
    pmats<-unlist(lapply(1:ncol(diffsimsmat),function(simcol){
    if(signflip){
      diffsimsmat[,simcol]<-diffsimsmat[,simcol]*-1
    }
    if (scaletomaximum){
      diffsimsmat[,simcol]<-rescale(diffsimsmat[,simcol],c(0,1))
    }
    matpoint<-min(modelpred$pred[abs(diffsimsmat[,simcol])<=p],na.rm=TRUE)
    return(matpoint)
    }))
    pmatsci<-data.frame(t(quantile(pmats[is.finite(pmats)],ciqnt, na.rm=T)))
    names(pmatsci)<-c("CI95low","CI95high")
    pmatsci$mean<-mean(pmats[is.finite(pmats)])
    pmatsci$sd<-sd(pmats[is.finite(pmats)],na.rm=TRUE)
    pmatsci$thresholdinterval<-p
    pmatsci$thresholdperunit<-p/interval_inc
    pmatsci$thresholdtype<-"entered/default"
    pmatsci$iterwithoutmatpoint<-length(which(is.infinite(pmats)))
    if (scaletomaximum && pmatsci$thresholdinterval==maximumderivproportion){
      pmatsci$thresholdtype<-"max non-linearity"
    }
    return(list(pmatsci=pmatsci,pmatsdf=data.frame(matpoints=pmats,thresholdinterval=p)))
    })
  derivmatpointsdf<-data.frame(dplyr::bind_rows(lapply(1:length(derivmatpoints),function(xi){derivmatpoints[[xi]]$pmatsci})))
  derivmatpointsdf$type<-"deriv"
  derivmatpointsdf$scaletomaximum<-scaletomaximum
  
  allrawmats<-data.frame(dplyr::bind_rows(lapply(1:length(derivmatpoints),function(xi){derivmatpoints[[xi]]$pmatsdf})))
  return(list(allmatpointsCI=derivmatpointsdf,matpointsims=allrawmats))
  }
  ####total mat points####
  if (!is.null(totalproportions)){
  fits<-gamsims$fits
  names(fits)<-seq(1,length(fits))
  fitsimsmat <- data.frame(dplyr::bind_rows(fits))
  
  totalmatpoints<-lapply(totalproportions,function(p){ ###apply thresholds
    pmats<-unlist(lapply(1:ncol(fitsimsmat),function(simcol){
      if(signflip){
        fitsimsmat[,simcol]<-fitsimsmat[,simcol]*-1
        fitsimsmat[,simcol]<-rescale(fitsimsmat[,simcol],c(0,1))
      }
      matpoint<-min(modelpred$pred[abs(fitsimsmat[,simcol])>=p],na.rm=TRUE)
      return(matpoint)
    }))
    pmatsci<-data.frame(t(quantile(pmats,ciqnt, na.rm=T)))
    names(pmatsci)<-c("low","high")
    pmatsci$mean<-mean(pmats)
    pmatsci$thresholdinterval<-p
    pmatsci$thresholdperunit<-p
    pmatsci$thresholdtype<-"entered/default"
    return(list(pmatsci=pmatsci,pmatsdf=data.frame(matpoints=pmats,thresholdinterval=p)))
  })
  totalmatpointsdf<-totalmatpoints[[1]]$pmatsci
  totalmatpointsdf$type<-"prop. total"
  totalmatpointsdf$scaletomaximum<-scaletomaximum
  
  ###pick up and return raw mat points for derivatives#####
  
  allmatpoints<-rbind(totalmatpointsdf,derivmatpointsdf)
  return(list(allmatpointsCI=allmatpoints,matpointsims=totalmatpoints[[1]]$pmatsdf))
}
}
####main gamm wrapper functions##########
gamm4_growthrate<-function(m, agevar, idvar = NULL, n.iterations = 10000, qnt = c(0.025, 
                                                                                 0.975)) 
{
  simdiff <- sim_diff1_from_gamm4(m, agevar, idvar, n.iterations = n.iterations)
  ci <- ci_from_simdiff1(simdiff$pred, simdiff$ages, qnt = qnt)
  ci$fit <- simdiff$fit
  return(ci)
}

mgcvgam_growthrate_maturationpoint<-function(m, agevar, idvar = NULL, n.iterations = 10000, qnt = c(0.025, 
                                                                                  0.975),interval_inc=.1,scaletomaximum=FALSE,usegratia=TRUE) 
{
  if (usegratia){
    require(gratia)
    #deriv_gratia<-fderiv(m)
    v <- m$model[, agevar]
    agevarintervaldata <- mgcvgampreddata(m,agevar,idvar,interval_inc=.1)
    CIlevel<-qnt[2]-qnt[1]
    f1deriv<-gratia::derivatives(m,term=agevar,partial_match=TRUE,interval="simultaneous",level=CIlevel,type="backward",newdata=agevarintervaldata)
    names(f1deriv)<-c("smooth","term","ages","mean_dff","se","crit","ci_low","ci_high")
    f1deriv$ages<-agevarintervaldata[,agevar]
    f1deriv$matpoint <-gam_maturation_pointmaxsig(f1deriv)
    ciclip<-clip_on_sig(f1deriv)
  } else {
  
  simdiff <- sim_from_mgcvgam(m, agevar, idvar, n.iterations = n.iterations,scaletomaximum=scaletomaximum)
  ci <- ci_from_simdiff1(simdiff$diffs, simdiff$ages, qnt = qnt)
  ci$fit <- simdiff$meanfit
  ci$matpoint <-gam_maturation_pointmaxsig(ci)
  ciclip<-clip_on_sig(ci)
  }
  return(ciclip)
}

# gratia_growthrate_maturationpoint<-function(m, agevar, idvar = NULL, n.iterations = 10000,interval_inc=.1, CIlevel=.95,scaletomaximum=FALSE) 
# {
#   require(gratia)
#   v <- m$model[, agevar]
#   agevarintervaldata <- data.frame(pred=seq(min(v), max(v), by=interval_inc))
#   names(agevarintervaldata)<-agevar
#   f1deriv <-gratia::fderiv(m,term=agevar,newdata=agevarintervaldata) 
#   f1derivCI<-confint(f1deriv,param=agevar,type="simultaneous",level=CIlevel,nsims=n.iterations)
#   names(f1derivCI)<-c("term","ci_low","ci_high","mean_dff")
#   f1derivCI$ages<-agevarintervaldata[,agevar]
#   f1derivCI$matpoint <-gam_maturation_pointmaxsig(f1derivCI)
#   ciclip<-clip_on_sig(f1derivCI)
#   return(ciclip)
# }

gam_maturation_pointmaxsig <- function(ci) {
  
  # when ci bounds include 0 (different sign), no longer signficant
  # clip out insignificant derivitive
  if (is.na(ci$ci_low[1])) ci <- ci[-1, ]
  
  # get mean_df_clip column
  if (! "mean_dff_clip" %in% names(ci)) ci <- clip_on_sig(ci)
  
  # find maturation point after the first signficant age
  onset_sig <- ci$ages[ci$mean_dff_clip != 0]
  maturation_pnt <- NA
  if (length(onset_sig)>0L && !all(is.na(onset_sig))) {
    #mat_points_idx <- ci$mean_dff_clip==0 & ci$ages > onset_sig[1]
    mat_points_idx<- which((ci$ages > max(ci$ages[ci$mean_dff_clip!=0]))) ###Ages greater than maximum significant (better defintion of maturation 20191130)
    if (length(mat_points_idx) > 0L && any(mat_points_idx))
      maturation_pnt <- min(ci$ages[mat_points_idx], na.rm=T) 
    if(length(onset_sig)>0L && !any(mat_points_idx))
      maturation_pnt<-max(ci$ages)
  }
  return(maturation_pnt)
}


sim_from_mgcvgam <- function(m, agevar, idvar=NULL,
                               n.iterations=10000, interval_inc=.1,scaletomaximum=FALSE) {
  
  if (class(m)[1]=="gamm"){
    m<-m$gam
    idvar<-NULL ###only use fixed effects portion of GAMM (no REs)
  }
  v <- m$model[, agevar]
  cond_list <- list(seq(min(v), max(v), by=interval_inc))
  pp <- data.frame(a=cond_list[[1]], b=Inf)
  # names should match what went into the model
  names(pp) <- c(agevar, idvar)
  
  # what if idvar is factor (Inf wont work)
  if (is.null(idvar)) {
    # do nothing. no idvar
  } else if (is.factor(m$model[, idvar])){
    # select idvar with the middle most random effect
    # random effects are coefficents like s(idvar).xxxxx
    # where xxxx is the index of the specific idvar factor name
    idvarpatt <- sprintf("s\\(%s\\)", idvar)
    idvarpatt. <- sprintf("s\\(%s\\).", idvar)
    randeff <- m$coefficients[ grep(idvarpatt, names(m$coefficients)) ]
    medval <- sort(randeff)[floor(length(randeff)/2)]
    med_re_name <- names(which(randeff == medval))
    median_idx <- gsub(idvarpatt., "", med_re_name)
    median_subj <- levels(m$model[, idvar])[as.numeric(median_idx)]
    warning("gam w/factor idvar, ",
            "setting the middle most random effect subject: ",
            median_subj)
    pp[, 2] <- median_subj
    
    # alternatively, select the first
    # pp[, 2] <- m$model[1, idvar]
  } else {
    warning("predition with continous (non-factor) idvar will give 'Inf' fit")
    # maybe pick middle value instead?
    # pp[, 2] <- mean(m$model[, idvar], na.rm=T)
  }
  
  # for all covars, pick out the mean
  for (cv in find_covars_gam(m$formula, agevar)) {
    x <- m$model[, cv]
    if (is.character(x) || is.factor(x) ){
      warning("gam w/factor covar, setting all sim to the first!")
      y <- x[1]
      # TODO: maybe pracma::Mode ?
    } else {
      y <- mean(x, na.rm=T)
    }
    pp[, cv] <- y
  }
  
  Xp <- predict(m, pp, type="lpmatrix")
  
  mu_beta <- coef(m)
  sigma_Vb <- vcov(m)
  # variance-covariance matrix of the main parameters  fitted model
  # used as: a positive-definite symmetric matrix specifying
  #  the covariance matrix of the variables.
  
  # set.seed(10)
  mrand <- MASS::mvrnorm(n.iterations, mu_beta, sigma_Vb)
  
  # ilink <- family(m)$linkinv
  # ilink<-m$family$linkinv()
  # only want inetercept and agevar
  keep_cols <- grep(paste0("Intercept|", agevar), dimnames(Xp)[[2]], value=T)
  Xp_agevar <- Xp[, keep_cols]
  mrand_agevar <- mrand[, keep_cols]
  
  # generate a whole bunch of plausable values, get the diff
  sims <- lapply(seq_len(n.iterations), function(i)  {
    fit <- m$family$linkinv((Xp_agevar %*% mrand_agevar[i, ]))
    if (scaletomaximum){
     fit<-rescale(fit,c(0,1))
    }
    dff <- c(NA, diff(fit))
    return(list(dffs=dff,fits=fit))
  })
  
  
  return(list(diffs=lapply(sims,"[[","dffs"),fits=lapply(sims,"[[","fits"),ages=pp[, 1], meanfit=predict(m, pp)))
}



ci_from_simdiff1 <- function(pred, ages, qnt=c(.025, .975)) {
  
  names(pred) <- 1:length(pred)
  mm <- t(dplyr::bind_rows(pred))
  
  # this is the ouptut !
  mean_dff <- apply(mm, 2, mean)
  ci <- apply(mm, 2, quantile, qnt, na.rm=T)
  colnames(ci) <- ages
  out <- data.frame(mean_dff=mean_dff, ages=ages)
  ci_out <- t(ci)
  dimnames(ci_out)[[2]] <- c("ci_low", "ci_high")
  return(cbind(out, ci_out))
  
  # NEVER REACHED -- left as bad documentation
  # old: return just ci and mean_dff
  return(list(ci=ci, mean_dff=mean_dff))
  
  # this is for fun
  ages[which.min(ci[1, ])]
  ages[which.min(ci[2, ])]
  
  plot(ages, mean_dff)
  for (i in 1:10) lines(ages, pred[[i]])
}

proportionofmaturationpointmcgvgam<-function(df,outcomevars,idvar=NULL,predvars,mformula,matfromsig=TRUE,jackknife=FALSE,jackknifesplits,jackknifecores){
  require(scales)
  pairs<-as.data.frame(expand.grid(outcomevars,predvars))
  names(pairs)<-c("outcome","pred")
  
  matfunction<-function(env = parent.frame()){
  allmatpoints<-lapply(1:nrow(pairs),function(p){
    print(pairs[p,])
    df$outcome<-scale(unlist(df[,as.character(pairs$outcome[p])]))
    df$pred<-unlist(df[,as.character(pairs$pred[p])])
    model<-mgcv::gam(mformula,data=df)
    if (!is.null(idvar)){
      model<-mgcv::gamm(mformula,random=list(id=~pred),data=df)
    }
    matpointsim<-proportionofmaxmaturation_mgcvgamsimCI(model,"pred",interval_inc=.1,scaletomaximum=FALSE,totalproportions = NULL,dervproportions=c(.05,.025))
    if(matfromsig){
      ci<-mgcvgam_growthrate_maturationpoint(model,"pred")
      sigmat<-data.frame(matpoint=c(mean(ci$matpoint)))
      sigmatdf<-data.frame(matrix(NA,nrow=1,ncol=ncol(matpointsim$allmatpointsCI)))
      names(sigmatdf)<-names(matpointsim$allmatpointsCI)
      sigmatdf$mean<-mean(ci$matpoint)
      sigmatdf$type<-"last sig"
      matpoints<-rbind(matpointsim$allmatpointsCI,sigmatdf)}
    
    matpoints$outcome<-as.character(pairs$outcome[p])
    matpoints$pred<-as.character(pairs$pred[p])
    rawsims<-matpointsim$matpointsims
    rawsims$outcome<-as.character(pairs$outcome[p])
    rawsims$pred<-as.character(pairs$pred[p])
    
    if (class(model)[1]=="gamm"){
    nps<-length(summary(model$gam)$s.pv)
      if (nps==1){
      matpoints$modelp<-summary(model)$s.pv 
      } else {  
    matpoints$modelp<-summary(model$gam)$s.pv[which(grepl("pred",row.names(summary(model$gam)$s.table)))]
    }
    } else {
    nps<-length(summary(model)$s.pv)
      if (nps==1){
        matpoints$modelp<-summary(model)$s.pv
      } else {
    matpoints$modelp<-summary(model)$s.pv[which(grepl("pred",row.names(summary(model)$s.table)))]
    }}
    return(list(matpointsummary=matpoints,rawmatpoints=rawsims))
  })
  #allmatpoints<-dplyr::bind_rows(allmatpoints)
  return(allmatpoints)
  }
  
  if (!jackknife){
  matpointsout<-matfunction()
  return(matpointsout)
  } else {
  fulldf<-df
  require(parallel)
  require(caret)
  jackknifefolds<-createFolds(row.names(fulldf),k=jackknifesplits,list = TRUE, returnTrain = FALSE)
  matpointlist<-mclapply(1:length(jackknifefolds),mc.cores=jackknifecores,function(foldi){
    fold<-jackknifefolds[[foldi]]
    df<-fulldf[fold,]
    matpointsout<-matfunction()
    matpointsout$fold<-foldi
    return(matpointsout)
  })
  matpointlistmat<-dplyr::bind_rows(matpointlist)  
  
  }
  
}


####plotting functions##
mgcvgam_ggplot<-function(predframe,derivframe,model,d,predvar,idvar,plotxlim=NULL,yvar = as.character(model$formula[2]), 
                         plotsavename = NA, xplotname = "Age", yplotname = yvar, draw_maturation = T, 
                         draw_points = T, show_all_fill = F, ci_plot = T,hexbin=F,hexbinnum=100){
  
  if (class(model)[1]=="gamm"){
    model<-model$gam
  }
  
  ci<-derivframe ###named "ci" for legacy with deriv functions
  
  require(ggplot2)
  if (!"gam" %in% class(model))
    stop("model given must be a gam model!")
  if (!"data.frame" %in% class(d)) 
    stop("d given must be a data.frame!")
  if (!"data.frame" %in% class(ci)) 
    stop("ci is not growthrate_gam() output")
  if (!yvar %in% names(model$model)) 
    stop(yvar, "not in model dataframe!")
  maturation_pnt <- mean(ci$matpoint) ###mat point as repeating column in ci
  if (is.na(maturation_pnt) && draw_maturation) {
    warning("No maturation point!")
    draw_maturation <- F
  }
  fill_column <- ifelse(show_all_fill, "mean_dff", "mean_dff_clip")
  deriv_range <- range(ci$mean_dff, na.rm = T)
  tile <- ggplot(ci[-1, ]) + aes_string(x = "ages", y = 1, 
                                        fill = fill_column) + geom_raster(interpolate = TRUE) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), 
                         limits = deriv_range) + xlab(sprintf("\n%s", xplotname))
  if (!is.null(plotxlim)){
    tile<-tile+xlim(plotxlim)
  }
  
  if (draw_maturation) 
    tile <- tile + geom_segment(linetype = 2, colour = "black", 
                                aes(x = maturation_pnt, xend = maturation_pnt, y = 0.5, 
                                    yend = 1.5))
    tilelegend<-cowplot::get_legend(tile)
  tile_luna <- lunaize_geomraster(tile) + theme(text = element_text(size = 36))
  modeldata <- data.frame(ydata = model$y, predvar = model$model[, 
                                                                predvar])
  condlist <- list(a = ci$ages)
  names(condlist) <- predvar

  ageplot <- ggplot(predframe) + aes_string(x = predvar, y = "fit") + 
    geom_line(colour = "black", size = 2) + ylab(yplotname) + 
    xlab(xplotname)
  
  if (ci_plot) {
    ageplot <- ageplot + geom_ribbon(aes(ymin = fit - CI, 
                                         ymax = fit + CI), alpha = 0.3)
  }
  
  if (draw_points){
    if (hexbin){
      require(hexbin)
      idvar<-NULL
    hexplot<-ggplot(data = modeldata)+geom_hex(aes(y = ydata, 
                                                    x = predvar),bins=hexbinnum)+
      scale_fill_gradientn(colours=c("grey80","black"),name = "Frequency",na.value=NA)
    
    hexlegend<-cowplot::get_legend(hexplot)
    
    ageplot<-hexplot+geom_line(data=predframe,aes_string(x = predvar, y = "fit"),
                                                         colour = "black", size = 2) + 
      ylab(yplotname) + xlab(xplotname)
      if (ci_plot){
        ageplot<-ageplot+geom_ribbon(data=predframe,aes_string(x=predvar,ymin = "fit - CI", 
                                                        ymax = "fit + CI"), colour="black",alpha = 0.3)
      }
    
    ageplot<-ageplot+theme(legend.position = "none")
    
    } else{
    ageplot <- ageplot + geom_point(data = modeldata, aes(y = ydata, 
                                                          x = predvar), alpha = 0.2)
    }
  }
  
  if (!is.null(plotxlim)){
    ageplot<-ageplot+xlim(plotxlim)
  }
  
  
  if (!is.null(idvar) && draw_points) 
    ageplot <- ageplot + geom_line(data = d, aes_string(y = yvar, 
                                                        group = idvar), alpha = 0.2)
  ageplot_luna <- LNCDR::lunaize(ageplot) + theme(text = element_text(size = 36), 
                                                  axis.title.x = element_blank(), axis.text.x = element_blank())
  if (hexbin){
    ageplot_luna<-ageplot_luna+theme(legend.position = "none")
  }
  
  g <- gam_growthrate_plot_combine(ageplot_luna, tile_luna, 
                                   plotsavename)
  list_of_plots <- list(tile = tile_luna, ageplot = ageplot_luna, 
                        both = g)
  if (hexbin){
  return(list(list_of_plots=list_of_plots,tilelegend=tilelegend,hexlegend=hexlegend))
  } else{
    return(list(list_of_plots=list_of_plots,tilelegend=tilelegend))
  }

}



lunaize_geomrasterxkeep<-function(x){
  x+
    theme_bw()+
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y     = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_blank(),
      legend.position  = "none")+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}

too_small <- function(x) abs(x) < 10^-15
clip_on_sig <- function(ci){
  # if confidence interval includes zero
  # signs of x and y will be different, -x * +y  < 0
  # or if both high and low are extremly close to zero
  not_sig <- ci$ci_low * ci$ci_high < 0 |
    (too_small(ci$ci_low) & too_small(ci$ci_high))
  ci$mean_dff_clip <- ci$mean_dff
  ci$mean_dff_clip[not_sig] <- 0
  return(ci)
}
lunaize_geomraster<-function(x){
  x+
    theme_bw()+
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y     = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_blank(),
      axis.title.x     = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.x      = element_blank(),
      legend.position  = "none")+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
lunaize_topplot<-function(x){
  x+
    theme_bw()+theme(text = element_text(size = 36))+
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title.y     = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.y      = element_blank(),
      panel.border = element_blank(),
      legend.position  = "none")
}

multi_mgcvgam_growthrate_multiplot<-function(df,outcomevars,predvars,plotxlim=c(8,34),idvar=NULL,mformula,zscale=TRUE,hexbin=FALSE,plotsavespecifier=NULL,topscatter=NULL,topmodelformula=NULL,totaltiles=NULL,
                                             toppredvar){
  ###wrapper script for gam_growthrate######
  pairs<-as.data.frame(expand.grid(outcomevars,predvars))
  names(pairs)<-c("outcome","pred")
  ggrout<-NULL
  for (p in 1:nrow(pairs)){
    if (zscale){
      df$outcome<-base::scale(unlist(df[,as.character(pairs$outcome[p])]))
    }else{
      df$outcome<-unlist(df[,as.character(pairs$outcome[p])])
    }
    df$pred<-unlist(df[,as.character(pairs$pred[p])])
    model<-mgcv::gam(mformula,data=df)
    gammodel<-model
    if (!is.null(idvar)){
      ctrl <- list(maxIter = 10000, msMaxIter = 10000)
      model <- tryCatch(
        {
          model<-mgcv::gamm(mformula,random=list(id=~1+pred),data=df,control = ctrl)
        },
        error = function(e){
          print("interceptsonly")
          model<-mgcv::gamm(mformula,random=list(id=~1),data=df,control = ctrl)
        }
      )
      gammodel<-model
      model<-model$gam
    }
    ggr<-mgcvgam_growthrate_maturationpoint(model,agevar="pred",idvar=idvar)
    ggr$var<-as.character(pairs$outcome[p])
    ggrout<-rbind(ggrout,ggr)
  }
  ggrouttoplot<-ggrout[!is.na(ggrout$mean_dff),]
  ggrouttoplotclip<-clip_on_sig(ggrouttoplot)
  ggrouttoplotclip$se<-(ggrouttoplotclip$ci_high-ggrouttoplotclip$mean_dff)/1.96
  ggrouttoplotclip$sd<-ggrouttoplotclip$se*sqrt(10000) #n interations from model ###change to sample size
  ggrouttoplotclip$tstat<-(ggrouttoplotclip$mean_dff-0)/(ggrouttoplotclip$sd/sqrt(10000))
  
  #ggrouttoplotclip$percentdiff<-ggrouttoplotclip$mean_dff/ggrouttoplotclip$fit
  ggrouttoplotclip$tstat[ggrouttoplotclip$mean_dff_clip==0]<-0
  
  #ggrouttoplotclip %>% group_by(var) %>% summarise(r=range(tstat)[1],r2=range(tstat)[2])
  deriv_range<-range(ggrouttoplotclip$tstat)
  fillvar="tstat"
  if (zscale){
    deriv_range<-range(ggrouttoplotclip$mean_dff)
    fillvar="mean_dff_clip"
  }
  
  if (sign(deriv_range[1])==sign(deriv_range[2])){
    if (sign(max(deriv_range))==1){
      deriv_range[deriv_range==min(deriv_range)]<-0
    } else {
      deriv_range[deriv_range==max(deriv_range)]<-0
    }
  }
  
  allplots<-list()
  for (vi in 1:length(unique(ggrouttoplotclip$var))){
    print(vi)
    v<-unique(ggrouttoplotclip$var)[vi]
    ci<-ggrouttoplotclip[ggrouttoplotclip$var==v,]
    tile <- ggplot(ci) + aes_string(x = "ages", y = 1, fill = fillvar) + geom_raster(interpolate = TRUE) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range)+xlim(plotxlim)
    #+ggtitle(v)
    if (vi !=length(unique(ggrouttoplotclip$var))){
      tile_luna <- lunaize_geomraster(tile) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())#,axis.text.x=element_text(color="white"),axis.ticks.x=element_line(color="white"))
    }else{
      print("final")
      tile_luna <- lunaize_geomraster(tile) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
    }
    allplots[[vi]]<-tile_luna
    legend <- cowplot::get_legend(tile)
    tilegrid<-cowplot::plot_grid(plotlist=allplots,ncol = 1)
    returnplot<-tilegrid
  }
  if (!is.null(totaltiles)){
    fillplots<-list()
    print(totaltiles)
    addtiles<-totaltiles-length(allplots)
    cifill<-ci
    cifill$tstat<-0
    cifill$mean_dff_clip<-0
    tilefill <- ggplot(cifill) + aes_string(x = "ages", y = 1, fill = fillvar) + geom_raster(interpolate = TRUE) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range)+xlim(plotxlim)
    tilefillluna <- lunaize_geomraster(tilefill) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
    for (ai in seq(1,addtiles)){
      print(ai)
      #ati<-length(allplots)+ai
      fillplots[[ai]]<-tilefillluna
    }
    tilegrid<-cowplot::plot_grid(plotlist=c(allplots,fillplots),ncol = 1)
    returnplot<-tilegrid
  }
  # if (!is.null(interleavetiles)){
  #   print(interleavetiles)
  #   cifill<-ci
  #   cifill$tstat<-0
  #   tilefill <- ggplot(cifill) + aes_string(x = "ages", y = 1, fill = "tstat") + geom_raster(interpolate = TRUE) +
  #     scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range)
  #   tilefillluna <- lunaize_geomraster(tilefill) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
  #   for (ai in seq(1,addtiles)){
  #     print(ai)
  #     ati<-length(allplots)+ai
  #     allplots[[ati]]<-tilefillluna
  #   }
  # }
  if (!is.null(plotsavespecifier)){
    plotname<-sprintf("%s.pdf",plotsavespecifier)
    save_plot(plotname, tilegrid, ncol = 1, base_height=8,base_asp = 1.1)
    dev.off()
    grid.newpage()
    plotname<-sprintf("%s.legend.pdf",plotsavespecifier)
    pdf(plotname,width=2,height=2)
    grid.draw(legend)
    dev.off()
  }
  
  if (!is.null(topscatter) && !is.null(topmodelformula)){
    print("adding top scatter based on top model and predvars[1]")
    topmodel<-mgcv::gam(topmodelformula,data=df)
    df$toppredvar<-df[,toppredvar]
    if (!is.null(idvar)){
      topmodel<-mgcv::gamm(topmodelformula,random=list(id=~toppredvar),data=df)
    }
    ggrtop<-mgcvgam_growthrate_maturationpoint(topmodel,agevar=predvars[1],idvar=idvar)
    toppred<-mgcvgampreddata(topmodel,predvar=predvars[1],idvar=idvar)
    scatter<-mgcvgam_ggplot(predframe=toppred,derivframe=ggrtop,model=topmodel,plotxlim=plotxlim,predvar=predvars[1],idvar=idvar,d=df,hexbin=hexbin)
    scatter<-list(lunaize_topplot(scatter$list_of_plots$ageplot))
    allplotstp<-allplots
    allplotstp2<-c(scatter,allplotstp)
    numbars=length(allplots)
    if (!is.null(totaltiles)){
      allplotstp2<-c(scatter,allplotstp,fillplots)
      numbars=totaltiles
    }
    topplottilegrid<-cowplot::plot_grid(plotlist=allplotstp2,ncol = 1,rel_heights=c(10,rep(1,numbars)))
    returnplot<-topplottilegrid
    if (!is.null(plotsavespecifier)){
      topplotname<-sprintf("%s.withtopplot.pdf",plotsavespecifier)
      save_plot(topplotname, topplottilegrid, ncol = 1, base_height=12,base_width =10)
    }
  }
  return(list(returnplot=returnplot,tilelegend=legend))
  
}
#####
multi_mgcvgam_growthrate_multiplot_rasteronly<-function(df,outcomevars,predvars,plotxlim=c(8,36),breaks=(c(10,15,20,25,30,35)),
                                                        idvar=NULL,mformula,zscale=TRUE,totaltiles=NULL,xscaleplot=TRUE, derivrangemanual=NULL,
                                                        derivcolourmanual=c("blue","red"),
                                                        xrangeticks=TRUE,xrangetickvals=NULL,datarangegrey=TRUE,devageguides=TRUE
                                                        ){ ###derivative raster only (no top plots)
  ###wrapper script for gam_growthrate######
  if(is.null(xrangetickvals)){xrangetickvals<-plotxlim}
  pairs<-as.data.frame(expand.grid(outcomevars,predvars))
  names(pairs)<-c("outcome","pred")
  ggrout<-NULL
  for (p in 1:nrow(pairs)){
    if (zscale){
      df$outcome<-base::scale(unlist(df[,as.character(pairs$outcome[p])]))
    }else{
      df$outcome<-unlist(df[,as.character(pairs$outcome[p])])
    }
    df$pred<-unlist(df[,as.character(pairs$pred[p])])
    model<-mgcv::gam(mformula,data=df)
    gammodel<-model
    if (!is.null(idvar)){
      ctrl <- list(maxIter = 10000, msMaxIter = 10000)
      model <- tryCatch(
        {
          model<-mgcv::gamm(mformula,random=list(id=~1+pred),data=df,control = ctrl)
        },
        error = function(e){
          print("interceptsonly")
          model<-mgcv::gamm(mformula,random=list(id=~1),data=df,control = ctrl)
        }
      )
      gammodel<-model
      model<-model$gam
    }
    ggr<-mgcvgam_growthrate_maturationpoint(model,agevar="pred",idvar=idvar)
    ggr$var<-as.character(pairs$outcome[p])
    ggrout<-rbind(ggrout,ggr)
  }
  ggrouttoplot<-ggrout[!is.na(ggrout$mean_dff),]
  ggrouttoplotclip<-clip_on_sig(ggrouttoplot)
  ggrouttoplotclip$se<-(ggrouttoplotclip$ci_high-ggrouttoplotclip$mean_dff)/2
  ggrouttoplotclip$sd<-ggrouttoplotclip$se*sqrt(10000) #n interations from model ###change to sample size
  ggrouttoplotclip$tstat<-(ggrouttoplotclip$mean_dff-0)/(ggrouttoplotclip$sd/sqrt(10000))
  
  ggrouttoplotclip$tstat[ggrouttoplotclip$mean_dff_clip==0]<-0
  deriv_range<-range(ggrouttoplotclip$tstat)
  fillvar="tstat"
  if (zscale){
    deriv_range<-range(ggrouttoplotclip$mean_dff_clip)
    fillvar="mean_dff_clip"
  }
  
  if (sign(deriv_range[1])==sign(deriv_range[2])){
    if (sign(max(deriv_range))==1){
      deriv_range[deriv_range==min(deriv_range)]<-0
    } else {
      deriv_range[deriv_range==max(deriv_range)]<-0
    }
  }
  
  if(!is.null(derivrangemanual)){
    deriv_range<-derivrangemanual
    print(sprintf("using manual deriv range %s %s",deriv_range[1],deriv_range[2]))
    print("values outside of range will be set to max via oob=squish")
  }
  
  allplots<-list()
  for (vi in 1:length(unique(ggrouttoplotclip$var))){
    print(vi)
    v<-unique(ggrouttoplotclip$var)[vi]
    ci<-ggrouttoplotclip[ggrouttoplotclip$var==v,]
    tile <- ggplot(ci) + aes_string(x = "ages", y = 1, fill = fillvar) + geom_raster(interpolate = TRUE) +
      scale_fill_gradient2(low = derivcolourmanual[1], mid = "white", high = derivcolourmanual[2], midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range,oob=squish)+xlim(plotxlim)
    #+ggtitle(v)
    if (vi !=length(unique(ggrouttoplotclip$var))){
      tile_luna <- lunaize_geomraster(tile) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())#,axis.text.x=element_text(color="white"),axis.ticks.x=element_line(color="white"))
    }else{
      print("final")
      tile_luna <- lunaize_geomraster(tile) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
    }
    if(datarangegrey){
      require(ggpattern)
      tile_luna<-tile_luna+geom_rect(xmin = xrangetickvals[1], xmax = min(df$pred,na.rm=TRUE),   ymin =0.5, ymax = 1.5,   fill = "grey88") 
      tile_luna<-tile_luna+geom_rect(xmin = max(df$pred,na.rm=TRUE), xmax = xrangetickvals[2],   ymin =0.5, ymax = 1.5,   fill = "grey88") 
    }
    if(xrangeticks){
      tile_luna<-tile_luna+geom_segment(linetype = 1, size=.5,colour = "black", aes(x = min(df$pred,na.rm=TRUE), xend =  min(df$pred,na.rm=TRUE), y = 0.5, yend = 1.5))
      tile_luna<-tile_luna+geom_segment(linetype = 1, size=.5,colour = "black", aes(x = max(df$pred,na.rm=TRUE), xend = max(df$pred,na.rm=TRUE), y = 0.5, yend = 1.5))
    }
    
    allplots[[vi]]<-tile_luna
    legend <- cowplot::get_legend(tile)
    tilegrid<-cowplot::plot_grid(plotlist=allplots,ncol = 1)
    returnplot<-tilegrid
  }
  fillplots<-NA
  pairswithfill<-pairs
  if (!is.null(totaltiles)){
    fillplots<-list()
    print(totaltiles)
    addtiles<-totaltiles-length(allplots)
    cifill<-ci
    cifill$tstat<-0
    cifill$mean_dff_clip<-0
    tilefill <- ggplot(cifill) + aes_string(x = "ages", y = 1, fill = fillvar) + geom_raster(interpolate = TRUE) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range)+xlim(plotxlim)
    tilefillluna <- lunaize_geomraster(tilefill) + theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
    for (ai in seq(1,addtiles)){
      print(ai)
      #ati<-length(allplots)+ai
      fillplots[[ai]]<-tilefillluna
    }
    tilegrid<-cowplot::plot_grid(plotlist=c(allplots,fillplots),ncol = 1)
    returnplot<-tilegrid
    fillforpais<-data.frame(matrix(ncol = ncol(pairs), nrow = addtiles))
    names(fillforpais)<-names(pairs)
    pairswithfill<-rbind(pairs,fillforpais)
  }
  tilescaleplot<-NA
  if(xscaleplot){
    ciscaleplot<-ci
    ciscaleplot$tstat<-0
    ciscaleplot$mean_dff_clip<-0
    tilescaleplot <- ggplot(ciscaleplot) + aes_string(x = "ages", y = 1, fill = fillvar) + geom_raster(interpolate = TRUE) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", breaks = sort(c(0, deriv_range)), limits = deriv_range)+
      scale_x_continuous(limits = plotxlim,breaks=breaks)
    tilescaleplot <- lunaize_geomrasterxkeep(tilescaleplot)+theme(text = element_text(size = 20))+theme(axis.title.x=element_blank())
    if(xrangeticks){
      tilescaleplot<-tilescaleplot+geom_segment(linetype = 1, colour = "black", aes(x = min(df$pred,na.rm=TRUE), xend =  min(df$pred,na.rm=TRUE), y = 0.5, yend = 1.5))
      tilescaleplot<-tilescaleplot+geom_segment(linetype = 1, colour = "black", aes(x = max(df$pred,na.rm=TRUE), xend = max(df$pred,na.rm=TRUE), y = 0.5, yend = 1.5))
    }
    if (devageguides){
      tilescaleplot<-tilescaleplot+geom_segment(linetype = 2, colour = "black", aes(x = 12, xend =  12, y = 0.5, yend = 1.5))
      tilescaleplot<-tilescaleplot+geom_segment(linetype = 2, colour = "black", aes(x = 15, xend=15, y = 0.5, yend = 1.5))
      tilescaleplot<-tilescaleplot+geom_segment(linetype = 2, colour = "black", aes(x = 18, xend=18, y = 0.5, yend = 1.5))
    }
    if(datarangegrey){
      tilescaleplot<-tilescaleplot+geom_rect(xmin = xrangetickvals[1], xmax = min(df$pred,na.rm=TRUE),   ymin =0.5, ymax = 1.5,   fill = "grey88") 
      tilescaleplot<-tilescaleplot+geom_rect(xmin = max(df$pred,na.rm=TRUE), xmax = xrangetickvals[2],   ymin =0.5, ymax = 1.5,   fill = "grey88") 
    }
    tilescaleplotlist<-list(tilescaleplot)
    if(!is.null(totaltiles)){
      returnplot<-cowplot::plot_grid(plotlist=c(allplots,fillplots,tilescaleplotlist),ncol = 1)
    }else{
      returnplot<-cowplot::plot_grid(plotlist=c(allplots,tilescaleplotlist),ncol = 1)
    }
  }
  

  return(list(returnplot=returnplot,tilelegend=legend,derivdata=ggrouttoplotclip,pairs=pairs,pairswithfill=pairswithfill,allplots=allplots,fillplots=fillplots,tilescaleplot=tilescaleplot))
  
}
#####



mgcv_plotcolumnwrapper<-function(df,outcomes,primaryvar,covars=NULL,idvar=NULL,correportmethod="kendall"){
  #vars
  df<-as.data.frame(df)
  df$primaryvar<-as.numeric(df[,primaryvar])
  primaryvar<-NULL
  #fml<-as.formula(sprintf("%s~s(%s)+%s",modelpairs$outcome[i],modelpairs$primaryvar[i],modelpairs$covars[i]))
  if (!is.null(covars)){
    modelpairs<-expand.grid(outcome=outcomes,primaryvar="primaryvar",covars=covars)
    modelpairs$formula<-unlist(lapply(1:nrow(modelpairs),function(i){
      fml<-as.formula(sprintf("%s~s(%s)+%s",modelpairs$outcome[i],modelpairs$primaryvar[i],modelpairs$covars[i])) 
    }))
  }else{
    modelpairs<-expand.grid(outcome=outcomes,primaryvar="primaryvar")
    modelpairs$formula<-unlist(lapply(1:nrow(modelpairs),function(i){
      fml<-as.formula(sprintf("%s~s(%s)",modelpairs$outcome[i],modelpairs$primaryvar[i])) 
    }))
  }
  #models
  if (is.null(idvar)){
    modelpreds<-lapply(1:nrow(modelpairs),function(i){
      model<-mgcv::gam(formula=modelpairs$formula[[i]],data=df)
      pred<-mgcvgampreddata(model,predvar="primaryvar")
      pred$outcome<-modelpairs$outcome[i]
      pred$modeltype<-"cross-sectional"
      pred$modelp<-summary(model)$s.table[row.names(summary(model)$s.table)=="s(primaryvar)",4]
      return(pred)
    })
    
  }else{
    print("mixed models")
    df$id<-as.factor(df[,idvar])
    ctrl <- list(maxIter = 10000, msMaxIter = 10000)
    modelpreds<-lapply(1:nrow(modelpairs),function(i){
      m<-NULL
      model<-NULL
      modeltype<-NULL
      print(modelpairs$outcome[i])
      model <- tryCatch(
        {
          m<-mgcv::gamm(formula=modelpairs$formula[[i]],
                     random=list(id=~primaryvar),data=df,control = ctrl)
          #modeltype="ranef slope+int"
          #return(list(m=m,modeltype=modeltype))
        },
        error = function(e){
          print("interceptsonly")
          m<-mgcv::gamm(formula=modelpairs$formula[[i]],
                     random=list(id=~1),data=df,control = ctrl)
          #modeltype="ranef int only"
          #return(list(m=m,modeltype=modeltype))
        }
      )
      m<-model
      pred<-mgcvgampreddata(m,predvar="primaryvar")
      pred$outcome<-modelpairs$outcome[i]
      #pred$modeltype<-modeltype
      pred$modelp<-summary(m$gam)$s.table[row.names(summary(m$gam)$s.table)=="s(primaryvar)",4]
      pred$edf<-summary(m$gam)$s.table[row.names(summary(m$gam)$s.table)=="s(primaryvar)",1]
      pred$cor<-cor(df$primaryvar,as.numeric(df[,as.character(modelpairs$outcome[i])]),method=correportmethod,
                    use="pairwise.complete.obs")
      return(pred)
    })
  }
  allpreds<-do.call(rbind,modelpreds)
}
mgcv_plotcolumnwrapper_fromvardf<-function(df,outcomesprimaryvarsdf,idvar=NULL,correportmethod="kendall",returnrandom=FALSE){
  #vars
  #for (modelrow in 1:nrow(outcomesprimaryvarsdf)){
  modelpreds<-lapply(1:nrow(outcomesprimaryvarsdf),function(modelrow){
    
    varsrow<-outcomesprimaryvarsdf[modelrow,]
    primaryvar<-as.character(varsrow[,"primaryvar"])
    
    df<-as.data.frame(df)
    df$primaryvar<-as.numeric(df[,primaryvar])
    pvarsave<-primaryvar
    primaryvar<-NULL
    
    covars<-grep("covar",names(outcomesprimaryvarsdf),value=TRUE)
    if (length(covars)!=0){
     
      smoothcovars<-grep("scovar",covars,value=TRUE) 
      nonsmoothcovars<-covars[!grepl("scovar",covars)]
      
      if (length(nonsmoothcovars)!=0){
        nonsmoothcovarchars<-paste(unlist(lapply(nonsmoothcovars,function(ci){
        return(as.character(varsrow[,ci]))
      })),collapse="+")
      }else{
        nonsmoothcovarchars<-nonsmoothcovars
      }
      
      if (length(smoothcovars)!=0){
        smoothcovarchars<-paste(unlist(lapply(smoothcovars,function(ci){
        if (length(unique(df[,as.character(varsrow[,ci])]))<10){
          return(sprintf("s(%s,k=%s)",as.character(varsrow[,ci]),length(unique(df[,as.character(varsrow[,ci])]))))
        }else{
        return(sprintf("s(%s)",as.character(varsrow[,ci])))
        }
      })),collapse="+")
      }else{
        smoothcovarchars<-smoothcovars
      }
      
      covarchars<-paste(c(smoothcovarchars,nonsmoothcovarchars),collapse="+")
      
      modelbase<-sprintf("%s~s(%s)+",varsrow$outcome,"primaryvar")
      fml<-as.formula(paste0(modelbase,covarchars)) 
      
    }else{
      fml<-as.formula(sprintf("%s~s(%s)",varsrow$outcome,"primaryvar"))
    }
    print(fml)
    #models
    if (is.null(idvar)){
      model<-mgcv::gam(formula=fml,data=df)
      pred<-mgcvgampreddata(model,predvar="primaryvar")
      pred$outcome<-varsrow$outcome
      pred$modeltype<-"cross-sectional"
      pred$primaryvarname<-pvarsave
      pred$modelp<-summary(model)$s.table[row.names(summary(model)$s.table)=="s(primaryvar)",4]
      pred$edf<-summary(model)$s.table[row.names(summary(model)$s.table)=="s(primaryvar)",1]
      pred$cor<-cor(df$primaryvar,as.numeric(df[,as.character(varsrow$outcome)]),method=correportmethod,
                    use="pairwise.complete.obs")
      ranefs<-NA
    }else{
      print("mixed models")
      df$id<-as.factor(df[,idvar])
      ctrl <- list(maxIter = 10000, msMaxIter = 10000)
      m<-NULL
      model<-NULL
      modeltype<-NULL
      print(as.character(varsrow$outcome))
      model <- tryCatch(
        {
          m<-mgcv::gamm(formula=fml,
                        random=list(id=~primaryvar),data=df,control = ctrl)
          #modeltype="ranef slope+int"
          #return(list(m=m,modeltype=modeltype))
        },
        error = function(e){
          print("interceptsonly")
          m<-mgcv::gamm(formula=fml,
                        random=list(id=~1),data=df,control = ctrl)
          #modeltype="ranef int only"
          #return(list(m=m,modeltype=modeltype))
        }
      )
      m<-model
      pred<-mgcvgampreddata(m,predvar="primaryvar")
      pred$outcome<-varsrow$outcome
      #pred$modeltype<-modeltype
      pred$modelp<-summary(m$gam)$s.table[row.names(summary(m$gam)$s.table)=="s(primaryvar)",4]
      pred$edf<-summary(m$gam)$s.table[row.names(summary(m$gam)$s.table)=="s(primaryvar)",1]
      pred$fstat<-summary(m$gam)$s.table[row.names(summary(m$gam)$s.table)=="s(primaryvar)",3]
      pred$cor<-cor(df$primaryvar,as.numeric(df[,as.character(varsrow$outcome)]),method=correportmethod,
                    use="pairwise.complete.obs")
      pred$primaryvarname<-pvarsave
      
      mranef<-ranef(m$lme)
      ranefdf<-data.frame(mranef$id)
      ranefdf[,idvar]<-row.names(ranefdf)
      if (length(which(names(ranefdf)=="primaryvar"))>0){
      names(ranefdf)[names(ranefdf)=="primaryvar"]<-paste0("ranslope_",pvarsave)
      }else{
        ranefdf[,paste0("ranslope_",pvarsave)]<-NA
      }
      names(ranefdf)[grepl("Intercept",names(ranefdf))]<-"ranintercept"
      ranefdf$outcome<-varsrow$outcome
    }
    return(list(pred=pred,ranefdf=ranefdf))
  })
  
  modelpredstemp<-lapply(modelpreds,"[[","pred")
  allpreds<-do.call(rbind,modelpredstemp)
  if(returnrandom){
  ranefstemp<-lapply(modelpreds,"[[","ranefdf")
  ranefs<-do.call(dplyr::bind_rows,ranefstemp)
  return(list(allpreds=allpreds,ranefs=ranefs))
  }else{
  return(allpreds)
  }
}

