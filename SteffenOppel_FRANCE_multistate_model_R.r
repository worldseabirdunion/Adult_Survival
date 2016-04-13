### Analysis of Yelkouan Shearwater survival rates using RMARK###
### script written by Steffen Oppel, RSPB, <steffen.oppel@rspb.org.uk>####

## PLEASE CITE: 
##Oppel, S., Raine, A.F., Borg, J.J., Raine, H., Bonnaud, E., Bourgeois, K., 
##Breton, A.R., 2011. Is the Yelkouan shearwater Puffinus yelkouan threatened by
##low adult survival probabilities? Biological Conservation 144: 2255-2263.




###################################################################################
#### IMPORT DATA FOR MULTISTATE MODEL, FRANCE 2004-2010 , 7 encounter occasions ###
###################################################################################


setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Malta\\Analysis\\Survival_analysis\\FRANCE_multistate_model")
setwd("C:\\STEFFEN\\WORK\\RSPB\\Malta\\Analysis\\Survival_analysis\\FRANCE_multistate_model")
library(RMark)

YESH<-import.chdata("FRANCE_Encounter_hist_alphabetic.inp", header = F, field.names= c("ch","island","frequency"), use.comments=TRUE)
YESH.proc<-process.data(YESH, model="Multistrata", begin.time=2004, groups=c("island"))
YESH.ddl<-make.design.data(YESH.proc, parameters=list(Psi=list(subtract.stratum=c("B","P"))))
YESH.ddl<-add.design.data(YESH.proc,YESH.ddl,parameter="S",type="age",bins=c(0,.5,7),name="transience")   ## adding a transience parameter
#YESH.ddl<-YESH.ddl$Psi[!(YESH.ddl$Psi$stratum=="B"&YESH.ddl$Psi$tostratum=="P"),]     ## setting Psi (breeder->prospector) = 0
#YESH.ddl$Psi$BtoP<-0
#YESH.ddl$Psi$BtoP[YESH.ddl$Psi$stratum=="P"&YESH.ddl$Psi$tostratum=="B"]=1



###################################################################################
#### SETTING UP MODEL RUN WITH PLAUSIBLE VARIATION IN 'S' AND 'p' #################
###################################################################################

do.analysis<-function()
{
S.island<-list(formula=~transience+island)
S.state<-list(formula=~transience+stratum)
S.state.island<-list(formula=~transience+island+stratum)
S.transience<-list(formula=~transience)
#p.island<-list(formula=~island)             # removed as provided poor model fit
p.stratum<-list(formula=~stratum)
p.island.stratum<-list(formula=~stratum+island)
#p.island.time<-list(formula=~-1+time:island)   # removed as provided poor model fit
p.time.stratum<-list(formula=~stratum+time)
#p.island.stratum.time<-list(formula=~time:stratum:island) # removed as parameters not identifiable
#p.time<-list(formula=~-1+time)                 # removed as provided poor model fit
#p.dot<-list(formula=~1)                        # removed as provided poor model fit
Psi.dot<-list(formula=~-1 + stratum:tostratum)
cml<-create.model.list("Multistrata")
model.list<-mark.wrapper(cml,data=YESH.proc,ddl=YESH.ddl)
return(model.list)
}
YESH_survival<-do.analysis()

#### ADJUST FOR OVERDISPERSION USING OUTPUT FROM U-CARE ############################
c<-62.427/37
YESH_survival<-adjust.chat(c,YESH_survival)


model.table<-model.table(YESH_survival,use.lnl=TRUE) 
model.table$EvidenceRatio<-2.71828182845904523536^(-0.5*model.table$DeltaQAICc) 
model.table




######################################################################################
############## ESTIMATE MODEL-AVERAGED SURVIVAL PARAMETERS FOR EACH GROUP ############
######################################################################################


## creating vector lists that specify the parameters of interest (only survival parameters)

par.index.france.breeders<-YESH.ddl$S[YESH.ddl$S$stratum=="B",]
par.index.france.prospectors<-YESH.ddl$S[YESH.ddl$S$stratum=="P",]
PQbreed<-min(as.numeric(row.names(par.index.france.breeders[par.index.france.breeders$island=="PQ",])))+2
PQprosp<-min(as.numeric(row.names(par.index.france.prospectors[par.index.france.prospectors$island=="PQ",])))+2
PCbreed<-min(as.numeric(row.names(par.index.france.breeders[par.index.france.breeders$island=="PC",])))+2
PCprosp<-min(as.numeric(row.names(par.index.france.prospectors[par.index.france.prospectors$island=="PC",])))+2


## model-average the results and compile the parameters of interest in a Table
## if you remove "drop=F" then models with ill-defined parameter estimates will be dropped, leading to slightly different output

parameter.estimates<-model.average(YESH_survival,parameter="Psi", vcv=T, drop=F)
parameter.estimates$estimates[parameter.estimates$estimates$stratum=='B',2:6]
parameter.estimates$estimates[parameter.estimates$estimates$stratum=='P',2:6]

parameter.estimates<-model.average(YESH_survival,parameter="S", vcv=T, drop=F)

Pooled_YESH_survival_estimates<-data.frame(Group=c("Porquerolles breeders", "Port-Cros breeders","Porquerolles prospectors", "Port-Cros prospectors"))
Pooled_YESH_survival_estimates$mean_survival[1]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PQbreed,2]
Pooled_YESH_survival_estimates$mean_survival[2]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PCbreed,2]
Pooled_YESH_survival_estimates$mean_survival[3]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PQprosp,2]
Pooled_YESH_survival_estimates$mean_survival[4]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PCprosp,2]
Pooled_YESH_survival_estimates$stand_error[1]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PQbreed,3]
Pooled_YESH_survival_estimates$stand_error[2]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PCbreed,3]
Pooled_YESH_survival_estimates$stand_error[3]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PQprosp,3]
Pooled_YESH_survival_estimates$stand_error[4]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PCprosp,3]
Pooled_YESH_survival_estimates$lower95CI[1]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PQbreed,4]
Pooled_YESH_survival_estimates$lower95CI[2]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PCbreed,4]
Pooled_YESH_survival_estimates$lower95CI[3]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PQprosp,4]
Pooled_YESH_survival_estimates$lower95CI[4]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PCprosp,4]
Pooled_YESH_survival_estimates$upper95CI[1]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PQbreed,5]
Pooled_YESH_survival_estimates$upper95CI[2]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PCbreed,5]
Pooled_YESH_survival_estimates$upper95CI[3]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PQprosp,5]
Pooled_YESH_survival_estimates$upper95CI[4]<-parameter.estimates$estimates[parameter.estimates$estimates$par.index==PCprosp,5]
Pooled_YESH_survival_estimates



########## MODEL-AVERAGING RECAPTURE PROBABILITY ###############

par.index.france.breeders<-YESH.ddl$p[YESH.ddl$p$stratum=="B",]
par.index.france.prospectors<-YESH.ddl$p[YESH.ddl$p$stratum=="P",]
PQbreed<-as.numeric(row.names(par.index.france.breeders[par.index.france.breeders$island=="PQ",]))+84
PQprosp<-as.numeric(row.names(par.index.france.prospectors[par.index.france.prospectors$island=="PQ",]))+84
PCbreed<-as.numeric(row.names(par.index.france.breeders[par.index.france.breeders$island=="PC",]))+84
PCprosp<-as.numeric(row.names(par.index.france.prospectors[par.index.france.prospectors$island=="PC",]))+84


## model-average the results and compile the parameters of interest in a Table
## if you remove "drop=F" then models with ill-defined parameter estimates will be dropped, leading to slightly different output

recap.estimates<-model.average(YESH_survival,parameter="p", vcv=T, drop=F)

Pooled_YESH_survival_estimates$mean_recapture[1]<-mean(recap.estimates$estimates[recap.estimates$estimates$par.index %in% PQbreed,2])
Pooled_YESH_survival_estimates$mean_recapture[2]<-mean(recap.estimates$estimates[recap.estimates$estimates$par.index %in% PCbreed,2])
Pooled_YESH_survival_estimates$mean_recapture[3]<-mean(recap.estimates$estimates[recap.estimates$estimates$par.index %in% PQprosp,2])
Pooled_YESH_survival_estimates$mean_recapture[4]<-mean(recap.estimates$estimates[recap.estimates$estimates$par.index %in% PCprosp,2])
Pooled_YESH_survival_estimates$min_recap[1]<-min(recap.estimates$estimates[recap.estimates$estimates$par.index %in% PQbreed,4])
Pooled_YESH_survival_estimates$min_recap[2]<-min(recap.estimates$estimates[recap.estimates$estimates$par.index %in% PCbreed,4])
Pooled_YESH_survival_estimates$min_recap[3]<-min(recap.estimates$estimates[recap.estimates$estimates$par.index %in% PQprosp,4])
Pooled_YESH_survival_estimates$min_recap[4]<-min(recap.estimates$estimates[recap.estimates$estimates$par.index %in% PCprosp,4])
Pooled_YESH_survival_estimates$max_recap[1]<-max(recap.estimates$estimates[recap.estimates$estimates$par.index %in% PQbreed,5])
Pooled_YESH_survival_estimates$max_recap[2]<-max(recap.estimates$estimates[recap.estimates$estimates$par.index %in% PCbreed,5])
Pooled_YESH_survival_estimates$max_recap[3]<-max(recap.estimates$estimates[recap.estimates$estimates$par.index %in% PQprosp,5])
Pooled_YESH_survival_estimates$max_recap[4]<-max(recap.estimates$estimates[recap.estimates$estimates$par.index %in% PCprosp,5])
Pooled_YESH_survival_estimates





######################################################################################
### MULTIPLY SURVIVAL ESTIMATES AND ESTIMATE STANDARD ERRORS USING DELTA METHDOD #####
######################################################################################


### EXTRACTING THE VCV MATRIX FOR SE ESTIMATION WITH DELTA METHOD
parameter.estimates$vcv

# build your own variance covariance matrix
breeder.vcv<-matrix(nrow=2, ncol=2)
breeder.vcv[1,1]<-parameter.estimates$vcv[PCbreed, PCbreed]
breeder.vcv[2,1]<-parameter.estimates$vcv[PCbreed, PQbreed]
breeder.vcv[1,2]<-parameter.estimates$vcv[PQbreed, PCbreed]
breeder.vcv[2,2]<-parameter.estimates$vcv[PQbreed, PQbreed]

prospector.vcv<-matrix(nrow=2, ncol=2)
prospector.vcv[1,1]<-parameter.estimates$vcv[PCprosp, PCprosp]
prospector.vcv[2,1]<-parameter.estimates$vcv[PCprosp, PQprosp]
prospector.vcv[1,2]<-parameter.estimates$vcv[PQprosp, PCprosp]
prospector.vcv[2,2]<-parameter.estimates$vcv[PQprosp, PQprosp]


### MULTIPLY SURVIVAL ESTIMATES AND ESTIMATE STANDARD ERRORS ####
library(msm)

## calculate the MEAN ANNUAL SURVIVAL PROBABILITY as an average across 3 years, with each year composed of 8 non-breeding months
mean_breeder_survival<-(Pooled_YESH_survival_estimates[1,2]+Pooled_YESH_survival_estimates[2,2])/2
mean_prospector_survival<-(Pooled_YESH_survival_estimates[3,2]+Pooled_YESH_survival_estimates[4,2])/2


## calculate the S.E. of the ANNUAL SURVIVAL PROBABILITY using the delta method
se_breed<-deltamethod(~((x1+x2)/2), Pooled_YESH_survival_estimates[1:2,2], breeder.vcv, ses=TRUE)
se_prosp<-deltamethod(~((x1+x2)/2), Pooled_YESH_survival_estimates[3:4,2], prospector.vcv, ses=TRUE)

## combine the output
output<-cbind("mean survival"=mean_breeder_survival, "standard error"=se_breed)
output<-rbind(output, cbind("mean survival"=mean_prospector_survival, "standard error"=se_prosp))

## export the results
write.table(output, "FRANCE_overall_mean_annual_survival_estimate.csv", row.names=F, sep=',')
 




######### TESTING DIFFERENCES IN SURVIVAL RATE USING Z-TESTS ###########

library(PASWR)
x<-Pooled_YESH_survival_estimates[4,2]
y<-output[1,1]
#n<-dim(subset(YESH, island=='PC'))[1]
n1=82   ### 171 breeders (82 PQ, 89 PC), 67 non-breeders (37 PQ, 30 PC)
n2=37
sigma.x<-Pooled_YESH_survival_estimates[4,3]*sqrt(n1)
sigma.y<-output[1,2]*sqrt(n2)
zsum.test(mean.x=x, sigma.x = sigma.x, n.x = n1, mean.y = y,sigma.y = sigma.y, n.y = n2, alternative = "two.sided", mu = 0,conf.level = 0.95)


x<-0.770
y<-0.880
n1=120   
n2=200
sigma.x<-0.019*sqrt(n1)
sigma.y<-0.018*sqrt(n2)
zsum.test(mean.x=x, sigma.x = sigma.x, n.x = n1, mean.y = y,sigma.y = sigma.y, n.y = n2, alternative = "two.sided", mu = 0,conf.level = 0.95)




######################################################################################
############## EXPORT RESULTS AND CLEAN UP ###########################################
######################################################################################


write.table(Pooled_YESH_survival_estimates,"Survival_estimates_France2004-2010_multistate.csv", sep=",", row.names=F)
write.table(model.table,"YESH_model_comparison_France2004-2010_multistate.csv", sep=",", row.names=F)
rm(top_model)
rm(YESH_survival)
cleanup(ask=FALSE)


######################################################################################
############## RUNNING THE TOP MODEL ###########################################
######################################################################################

#top_model<-mark(data=YESH.proc,ddl=YESH.ddl,model="Multistrata",model.parameters=list(S=list(formula=~transience:island:stratum), p=list(formula=~-1 + stratum:island), Psi=list(formula=~-1+stratum:tostratum)))
#summary(top_model)

