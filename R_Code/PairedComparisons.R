WLdef <- function(res, base, measure, sig){
  pres <- pairedComparisons(res, base, maxs=rep(TRUE,3), p.value=sig)
  WL <-matrix(0,ncol=5, nrow=(length(workflowNames(res))-1))
  colnames(WL) <- c("Win", "sigWin", "Loss", "SigLoss","Tie")
  rownames(WL) <- setdiff(workflowNames(res), base)
  # for each data set count wins/losses
  for (e in 1:length(taskNames(res))){
    for(nm in setdiff(workflowNames(res), base)){
      if(pres[[measure]]$WilcoxonSignedRank.test[nm,2,e]<0){
        WL[nm,1] <- WL[nm,1]+1
        if(pres[[measure]]$WilcoxonSignedRank.test[nm,3,e]<sig) WL[nm,2] <- WL[nm,2]+1
      } else if(pres[[measure]]$WilcoxonSignedRank.test[nm,2,e]==0){
        WL[nm,5] <- WL[nm,5]+1
      } else {
        WL[nm,3] <- WL[nm,3]+1
        if(pres[[measure]]$WilcoxonSignedRank.test[nm,3,e]<sig) WL[nm,4] <- WL[nm,4]+1
      } 
    }
  }
  WL
}

WLdef(exp,"mc.lm","F1",0.05)[c(1,4,7),]
WLdef(exp,"mc.svm","F1",0.05)[c(11,14,17),]
WLdef(exp,"mc.mars","F1",0.05)[c(21,24,27),]
WLdef(exp,"mc.rf","F1",0.05)[c(31,34,37),]
WLdef(exp,"mc.rpart","F1",0.05)[c(41,44,47),]

WLdef(exp,"mc.lm_UNDERB","F1",0.05)[c(2,3),]
WLdef(exp,"mc.svm_UNDERB","F1",0.05)[c(12,13),]
WLdef(exp,"mc.mars_UNDERB","F1",0.05)[c(22,23),]
WLdef(exp,"mc.rf_UNDERB","F1",0.05)[c(32,33),]
WLdef(exp,"mc.rpart_UNDERB","F1",0.05)[c(42,43),]

WLdef(exp,"mc.lm_OVERB","F1",0.05)[c(5,6),]
WLdef(exp,"mc.svm_OVERB","F1",0.05)[c(15,16),]
WLdef(exp,"mc.mars_OVERB","F1",0.05)[c(25,26),]
WLdef(exp,"mc.rf_OVERB","F1",0.05)[c(35,36),]
WLdef(exp,"mc.rpart_OVERB","F1",0.05)[c(45,46),]

WLdef(exp,"mc.lm_SMOTEB","F1",0.05)[c(8,9),]
WLdef(exp,"mc.svm_SMOTEB","F1",0.05)[c(18,19),]
WLdef(exp,"mc.mars_SMOTEB","F1",0.05)[c(28,29),]
WLdef(exp,"mc.rf_SMOTEB","F1",0.05)[c(38,39),]
WLdef(exp,"mc.rpart_SMOTEB","F1",0.05)[c(48,49),]

WLdef(exp,"mc.arima","F1",0.05)
WLdef(exp,"mc.BDES","F1",0.05)
