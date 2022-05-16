####LUNA####
allefvars<-c("Anti_CRLat","Anti_CRR","Mix_CRLat","Mix_CRR","SOC.Problems.solved.in.minimum.moves",
             "SOC.Overallmeaninitialthinkingtime","DMS.PC","DMS.Median.correct.latency","SSP.Span.length",
             "nfixbreak","best_acc_m_exclude","first_lat_m_exclude")

accvars<-c("Anti_CRR","Mix_CRR","DMS.Percent.correct","SSP.Span.length","nfixbreak")
accvarsperc<-c("Anti_CRR","Mix_CRR","DMS.PC")
accvarscount<-c("SSP.Span.length","nfixbreak")

###twelve measuers six tasks#####


###NCANDA####
allefvars<-c("cnp_pcet_pcet_acc2","cnp_pcet_pcetrtcr","cnp_sfnb2_sfnb_mcr","cnp_sfnb2_sfnb_mrtc","cnp_spcptnl_scpt_tp","cnp_spcptnl_scpt_tprt","stroop_total_mean")

accvars<-c("cnp_pcet_pcet_acc2","cnp_sfnb2_sfnb_mcr","cnp_spcptnl_scpt_tp")
latvars<-c("cnp_pcet_pcetrtcr","cnp_sfnb2_sfnb_mrtc","cnp_spcptnl_scpt_tprt","stroop_total_mean")
#1 non CNP task and measure 
####NKI#####
accvars<-c("PCET_PCET_ACC2","SLNB2_SLNB_MCR","SPCPTNL_SCPT_TP","TOWacc","DFLacc")
latvars<-c("PCET_PCETRTCR","SLNB2_SLNB_MRTC","SPCPTNL_SCPT_TPRT","CWIlat","TMTlat")

allefvars<-c(accvars,latvars)
#4 Non CNP tasks and 4 Measures

###PNC###
groups<-data.frame(acc=c("Accuracycomposite","PCET_ACC2","PCPT_T_TP","LNB_MCR"),
                   lat=c("Latencycomposite","PCET_RTCR","PCPT_T_TPRT","LNB_MRTC"),
                   
####six meausres three tasks from WEb CNP#####


#total meaurs 6+4+1+12
#total tasks 6+4+1+5

