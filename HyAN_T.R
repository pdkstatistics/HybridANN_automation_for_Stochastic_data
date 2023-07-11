## PACKAGES 
##############################################################
easypackages::packages(GP)
easypackages::packages(P.Hybrid)
options(max.print=999999)



## DATA SET
##############################################################
#DT=xlsx::read.xlsx(file.choose(),sheetIndex = 1)
#DT1=read.table("clipboard",header=TRUE)
DT=DT1
dim(DT)
variable.names(DT)
TS1=ts(DT, start=c(2001, 10), frequency=12)#
TS1

#=============#
#XXX=1:2
XXX=1:dim(DT)[2]   # Selected Variables
#=============#




## IMPUTATION 
##############################################################
# imputing the missing values
#------------------------------------------
#plot(ts(TS1))
TS=round(na_seadec(TS1,algorithm ="mean"),2)
#plot(ts(TS))



## MODELLING 
#############################################################=

# Inputs
L.HyAN.MDLs=list()
HyAN.D.Tr=c()
HyAN.D.Tt=c()
L.HyAN.FITTED=list()

for (i in XXX) {
    
    # DATA SPLITTING
    #--------------------------------------------------------------------------
    YY =TS[,i]
    T=round(length(YY)*0.85) #####
    Tr=head(YY,T)
    Tt=tail(YY,length(YY)-T)
    
    
    # HYBRID MODELLING
    #--------------------------------------------------------------------------
    L.HyAN.MDLs=c(L.HyAN.MDLs,list(hybridModel(Tr,models = "an")))
    names(L.HyAN.MDLs)[i]=variable.names(TS)[i]
    
    
    ## MODEL DIAGNOSTICS
    #--------------------------------------------------------------------------
    ## Training data ######
    # Tabling the estimates
    X.Tr=tibble(Y =as.vector(Tr),
                Ey=as.vector(L.HyAN.MDLs[[i]]$fitted),
                Residuals=as.vector(L.HyAN.MDLs[[i]]$residuals))
    X.Tr=drop_na(X.Tr)
    
    #  Diagnostic statistics 
    HyAN.D.Tr=rbind(HyAN.D.Tr,
                      as.data.frame(with(X.Tr,
                                         list(V.Name=variable.names(TS)[i],
                      Arima=as.character(L.HyAN.MDLs[[i]]$auto.arima),
                      NN=L.HyAN.MDLs[[i]]$nnetar$method,
                                              MAE=DescTools::MAE(Ey,Y),
                                              MAPE=DescTools::MAPE(Ey,Y),
                                              SMAPE=DescTools::SMAPE(Ey,Y),
                                              MSE=DescTools::MSE(Ey,Y),
                                              RMSE=DescTools::RMSE(Ey,Y),
                                              TheilU=DescTools::TheilU(Y,Ey)))))
    
    
    ## Testing data ######
    # Tabling the estimates
    X.Tt=tibble(Y =as.vector(Tt),
                Ey=as.data.frame(forecast::forecast(L.HyAN.MDLs[[i]],
                                h=length(Tt)))$`Point Forecast`,
                Residuals=Y-Ey)
    
    #  Diagnostic statistics 
    HyAN.D.Tt=rbind(HyAN.D.Tt,
                      as.data.frame(with(X.Tt,
                                         list(V.Name=variable.names(TS)[i],
                                              MAE=DescTools::MAE(Ey,Y),
                                              MAPE=DescTools::MAPE(Ey,Y),
                                              SMAPE=DescTools::SMAPE(Ey,Y),
                                              MSE=DescTools::MSE(Ey,Y),
                                              RMSE=DescTools::RMSE(Ey,Y),
                                              TheilU=DescTools::TheilU(Y,Ey)))))
    
    ## FITTED (E(Y) of selected models for full_data[YY])
    #----------------------------------------------------
    L.HyAN.FITTED=c(L.HyAN.FITTED,
                    list(tibble(Y=YY,
                         Ey=c(L.HyAN.MDLs[[i]]$fitted,X.Tt$Ey))))
    names(L.HyAN.FITTED)[i]=variable.names(TS)[i]
    
}
L.HyAN.MDLs
HyAN.D.Tr
HyAN.D.Tt
L.HyAN.FITTED



## RESIDUAL DIAGNOSTICS ####################################

HyAN.RES.AC=c()
HyAN.RES.NOR=c()
HyAN.RES.RUNS=c()


for (i in XXX) {
    
    RESID=L.HyAN.MDLs[[i]]$residuals
    
    
    # 1. AUTOCORRELATION
    #-------------------------------------------------------
    # 1.1 Box-Pierce test
    # 1.2 Box-Ljung test
    HyAN.RES.AC=rbind(HyAN.RES.AC,
                   data.frame(BP.stat=Box.test(RESID,lag=10,
                                               fitdf=0)$statistic,
                              BP.p   =Box.test(RESID,lag=10,
                                               fitdf=0)$p.value,
                              BP.inf =if(Box.test(RESID,lag=10,
                                                  fitdf=0)$p.value<=0.05){
                                  "Autocorrelation"}else{"No_Autocorrelation"},
                              
                              
                              BL.stat=Box.test(RESID,lag=10,
                                               fitdf=0,type="Lj")$statistic,
                              BL.p   =Box.test(RESID,lag=10,
                                               fitdf=0,type="Lj")$p.value,
                              BL.inf =if(Box.test(RESID,lag=10,fitdf=0,
                                                  type="Lj")$p.value<=0.05){
                                  "Autocorrelation"}else{"No_Autocorrelation"}))
    
    rownames(HyAN.RES.AC)[i]=variable.names(TS)[i]
    
    
    
    # 2. NORMALITY
    #-------------------------------------------------------
    
    # Shapiro-Wilks Test
    HyAN.RES.NOR=rbind(HyAN.RES.NOR,
                    data.frame(W.stat=shapiro.test(RESID)$statistic,
                               p=round(shapiro.test(RESID)$p.value,3),
                               inf=if(shapiro.test(RESID)$p.value<=0.05){
                                   "Non-Normal"}else{"Normal"}))
    rownames(HyAN.RES.NOR)[i]=variable.names(TS)[i]
    
    
    # 3. RANDOMNESS
    #-------------------------------------------------------
    
    # Runs test
    runs.test=randtests::runs.test
    HyAN.RES.RUNS=rbind(HyAN.RES.RUNS,
                     data.frame(U=runs.test(RESID)$statistic,
                                Runs=runs.test(RESID)$runs,
                                p=round(runs.test(RESID)$p.value,3),
                                inf=if(runs.test(RESID)$p.value<=0.05){
                                    "Non-Random"}else{"Random"}))
    rownames(HyAN.RES.RUNS)[i]=variable.names(TS)[i]
    
}

HyAN.RES.AC
HyAN.RES.NOR
HyAN.RES.RUNS

## PLOTS #################################################

#=============#
#VAR.NAMES=as.vector(read.table("clipboard",header=FALSE))
#=============#

for (i in XXX) {
    
    PLOTS=as.data.frame(L.HyAN.MDLs[[i]]$residuals)
    colnames(PLOTS)=variable.names(TS)[i]
    
    ##LOCATION
    png(file=paste("E:/STATISTICS/PhD/THESIS/STUDIES/RESULTS/GRAPHS/",i,".",
                   variable.names(PLOTS)[1], ".png", sep=""),
        width=700, height=175,pointsize=100,res=96)
    
    A=ggPacf(PLOTS[,1], lag.max =12)+
        labs(title = "")+
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5))
    
    
    B=ggAcf(PLOTS[,1], lag.max =12)+
        labs(title = "")+
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5))
    
    
    C=ggdraw() +
        draw_plot(B, x = 0, y = 0, width = 0.5, height = 1) +
        draw_plot(A, x = 0.5, y = 0, width =0.5, height = 1) +
        draw_label(VAR.NAMES[i,],fontface="bold",
                   x = 0.5, y = 0.97, 
                   hjust = 0.5, vjust = 0.5,
                   size = 10)
    
    print(C)
    dev.off()
}

#-----------------------------------------------------------------------------

sink("E:/STATISTICS/PhD/THESIS/STUDIES/RESULTS/GRAPHS/HyAN_FINAL.txt")  # START
sink()  # END

## EXPORT
cat("HYBRID ARIMA+ANN\n======================\n")
cat("TRAINING ACCURACY\n======================\n")
HyAN.D.Tr %>% mutate_if(is.numeric, round, digits=4)
cat("TESTING ACCURACY\n======================\n")
HyAN.D.Tt %>% mutate_if(is.numeric, round, digits=4)
cat("FITTED\n======================\n")
L.HyAN.FITTED %>% 
cat("RESIDUAL DIAGNOSTICS \n ========================\n")
cat("Box-Pierce test \nBox-Ljung test \n ========================\n")
HyAN.RES.AC %>% mutate_if(is.numeric, round, digits=4)
cat("Shapiro-Wilks Test \n========================\n")
HyAN.RES.NOR %>% mutate_if(is.numeric, round, digits=4)
cat("RANDOMNESS \n========================\n")
HyAN.RES.RUNS %>% mutate_if(is.numeric, round, digits=4)
################################################################################