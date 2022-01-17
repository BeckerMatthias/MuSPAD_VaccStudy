#Author: Barbora Kessel, 
#Affiliation: Helmholtz Centre for Infection Research, Dept. of Epidemiology, 
#             Braunschweig, Germany
#Contact: barbora.kessel@helmholtz-hzi.de
#Date: 2022-01-12


#This code generates the supplementary tables 4, 5, and 6 for the publication:
#Comparative magnitude and persistence of humoral SARS-CoV-2 vaccination responses
#in the adult population in Germany
#Rstudio Version 1.4.1717 was used with R version 4.1.2


#needed libraries: lawstat, nlme

#read in data for mix-and-match cohort
dat<-read.csv("data_mix_and_match_shared.csv",h=T,sep=",",stringsAsFactors = F)

#create some variables for convenience
dat$male<-as.numeric(dat$sex=="male")

dat$impf<-rep("",nrow(dat))
dat$impf[dat$vacc_type=="Pfizer-Pfizer"]<-"P/P"
dat$impf[dat$vacc_type=="Moderna-Moderna"]<-"M/M"
dat$impf[dat$vacc_type=="Astra-Astra"]<-"A/A"
dat$impf[dat$vacc_type=="Astra-Pfizer"]<-"A/P"
dat$impf[dat$vacc_type=="Astra-Moderna"]<-"A/M"
dat$impf[dat$vacc_type=="Janssen"]<-"J"
table(dat$impf,useNA="ifany")

dat$period<-as.numeric(dat$time_since_full_vacc>=28)

dat$impf_factor<-as.factor(dat$impf)

#separate previously infected probands into dat1
wh<-which(dat$test_reported_positiv==1 | dat$N>1)
dat1<-dat[wh,]
dat<-dat[-wh,]

############Supplementary table 4################
library(nlme)

#remove J and A/M
datsel<-dat[dat$impf!="J" & dat$impf!="A/M",]
datsel$impf_factor<-factor(datsel$impf)

#logit-transform the ACE 2 binding inhibitions
#negative values set to 0.001
datsel$WT2<-pmax(0.001,datsel$WT)
datsel$logitWT2<-log(datsel$WT2/(1-datsel$WT2))
datsel$BETA2<-pmax(0.001,datsel$BETA)
datsel$DELTA2<-pmax(0.001,datsel$DELTA)
datsel$logitBETA2<-log(datsel$BETA2/(1-datsel$BETA2))
datsel$logitDELTA2<-log(datsel$DELTA2/(1-datsel$DELTA2))
datsel$logitALPHA<-log(datsel$ALPHA/(1-datsel$ALPHA))
datsel$logitGAMMA<-log(datsel$GAMMA/(1-datsel$GAMMA))

#merge immunosuppression and cancer into one condition
datsel$other2<-as.numeric(datsel$immune_supp  | 
                            datsel$cancer)

#Due to privacy protection issues, we do not publicly share the exact age of
#our probands. In order to run the code we used for Supplementary table 4,
#the following ARTIFICIAL AGE variable may be created. However, note that 
#the estimates and p-values WON'T BE THE SAME as reported in the Supplement!

set.seed(1002)
datsel$age<-rep(NA,nrow(datsel))
datsel$age[datsel$age_group=="18-25"]<-sample(18:25,sum(datsel$age_group=="18-25"),
                                               replace=T)
datsel$age[datsel$age_group=="26-45"]<-sample(26:45,sum(datsel$age_group=="26-45"),
                                               replace=T)
datsel$age[datsel$age_group=="46-65"]<-sample(46:65,sum(datsel$age_group=="46-65"),
                                               replace=T)
datsel$age[datsel$age_group=="66-79"]<-sample(66:79,sum(datsel$age_group=="66-79"),
                                               replace=T)
datsel$age[datsel$age_group==">79"]<-sample(80:99,sum(datsel$age_group==">79"),
                                               replace=T)
#create indicator for age over 70
datsel$age70<-as.numeric(datsel$age>70)

###fit the model for WT (part A of the Table)
f1<-lme(logitWT2~(male+period+age+cardiovascular+
                          hypertension+
                          diabetes+
                          lung_disease+
                          other2):impf_factor +impf_factor-1,
              random=~1|Plate,data=datsel[!is.na(datsel$cardiovascular),],
              weights=varIdent(form=~1|age70*impf_factor))

intervals(f1,which="var-cov")
aux<-summary(f1)

#print results for chosen vaccination scheme
chosen<-"P/P" # or "M/M" or "A/P" or "A/A"
data.frame(est=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("Value")],digits=2),
            sd=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("Std.Error")],digits=3),
            p_val=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("p-value")],digits=3))

###fit the model for ALPHA (part B of the Table)
f1<-lme(logitALPHA~(male+period+age+cardiovascular+
                    hypertension+
                    diabetes+
                    lung_disease+
                    other2):impf_factor +impf_factor-1,
        random=~1|Plate,data=datsel[!is.na(datsel$cardiovascular),],
        weights=varIdent(form=~1|age70*impf_factor))

intervals(f1,which="var-cov")
aux<-summary(f1)

#print results for chosen vaccination scheme
chosen<-"P/P" # or "M/M" or "A/P" or "A/A"
data.frame(est=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("Value")],digits=2),
           sd=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("Std.Error")],digits=3),
           p_val=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("p-value")],digits=3))

###fit the model for BETA (part C of the Table)
f1<-lme(logitBETA2~(male+period+age+cardiovascular+
                    hypertension+
                    diabetes+
                    lung_disease+
                    other2):impf_factor +impf_factor-1,
        random=~1|Plate,data=datsel[!is.na(datsel$cardiovascular),],
        weights=varIdent(form=~1|age70*impf_factor))

intervals(f1,which="var-cov")
aux<-summary(f1)

#print results for chosen vaccination scheme
chosen<-"P/P" # or "M/M" or "A/P" or "A/A"
data.frame(est=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("Value")],digits=2),
           sd=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("Std.Error")],digits=3),
           p_val=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("p-value")],digits=3))

###fit the model for GAMMA (part D of the Table)
f1<-lme(logitGAMMA~(male+period+age+cardiovascular+
                    hypertension+
                    diabetes+
                    lung_disease+
                    other2):impf_factor +impf_factor-1,
        random=~1|Plate,data=datsel[!is.na(datsel$cardiovascular),],
        weights=varIdent(form=~1|age70*impf_factor))

intervals(f1,which="var-cov")
aux<-summary(f1)

#print results for chosen vaccination scheme
chosen<-"P/P" # or "M/M" or "A/P" or "A/A"
data.frame(est=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("Value")],digits=2),
           sd=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("Std.Error")],digits=3),
           p_val=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("p-value")],digits=3))

###fit the model for DELTA (part E of the Table)
f1<-lme(logitDELTA2~(male+period+age+cardiovascular+
                    hypertension+
                    diabetes+
                    lung_disease+
                    other2):impf_factor +impf_factor-1,
        random=~1|Plate,data=datsel[!is.na(datsel$cardiovascular),],
        weights=varIdent(form=~1|age70*impf_factor))

intervals(f1,which="var-cov")
aux<-summary(f1)

#print results for chosen vaccination scheme
chosen<-"P/P" # or "M/M" or "A/P" or "A/A"
data.frame(est=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("Value")],digits=2),
           sd=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("Std.Error")],digits=3),
           p_val=round(aux$tTable[grepl(chosen,row.names(aux$tTable)),c("p-value")],digits=3))

##fit model with interactions, i.e. allow coefficients vary
##in the peak and the plateau periods

#merge immunosuppression, cancer and lung_disease into one condition
datsel$other3<-as.numeric(datsel$immune_supp  | 
                            datsel$cancer | datsel$lung_disease)

f2<-lme(logitWT2~(male+age+cardiovascular+
                             hypertension+
                             diabetes+
                             other3+period+
                             male*period+age*period+cardiovascular*period+
                             hypertension*period+
                             diabetes*period+
                             other3*period):impf_factor +impf_factor-1,
              random=~1|Plate,data=datsel[!is.na(datsel$cardiovascular),],
              weights=varIdent(form=~1|age70*impf_factor))

aux<-summary(f2)

#look at significance of coefficients for chosen vaccination scheme in the two periods
# the choice below shows the p-values for age reported in the main text

chosen<-"M/M"
wh<-which(grepl(chosen,row.names(aux$tTable)))
L<-cbind(diag(7),diag(7))
stds<-sqrt(diag(L%*%f2$varFix[wh,wh]%*%t(L)))
data.frame(est1=round(aux$tTable[wh[1:7],c("Value")],digits=2),
            sd1=round(aux$tTable[wh[1:7],c("Std.Error")],digits=3),
            p_val1=round(2*(1-pnorm(abs(aux$tTable[wh[1:7],c("Value")]/aux$tTable[wh[1:7],c("Std.Error")]))),
                  digits=3),
            est2=round(aux$tTable[wh[1:7],c("Value")]+aux$tTable[wh[8:14],
                                                            c("Value")],digits=2),
            sd2=round(stds,digits=3),
            p_val2=round(2*(1-pnorm(abs(aux$tTable[wh[1:7],c("Value")]+aux$tTable[wh[8:14],
                                                                           c("Value")])/stds)),digits=3))



############Supplementary table 5################

#fix order of comparison
comp_ids<-cbind(c("J","J","J","J","J","A/A","A/A","A/A","A/A",
                  "A/P","A/P","A/P","A/M","A/M","P/P"),
                c("M/M","P/P","A/M","A/P","A/A","M/M","P/P","A/M","A/P",
                  "M/M","P/P","A/M","M/M","P/P","M/M"))
feature<-"Spike" #or "RBD" or "S1" or "S2"
res<-data.frame(comp=rep("",nrow(comp_ids)),
                pvals=rep(NA,nrow(comp_ids)),
                est=rep(NA,nrow(comp_ids))
                )
k<-1
for (i in 1: nrow(comp_ids)){
  res[k,1]<-paste(comp_ids[i,1]," - ",comp_ids[i,2],sep="")
  aux<-lawstat::brunner.munzel.test(dat[dat$impf==comp_ids[i,1],feature],
                                    dat[dat$impf==comp_ids[i,2],feature])
  
  res[k,2]<-aux$p.value
  res[k,3]<-aux$estimate
  k<-k+1
  
}

#sort the p-values and compare them to the Bonferonni-Holm sequence of levels
aux<-sort(res[,2],index.return=T)
holm_seq<-0.05*1/(nrow(res):1)
signif<-res[aux$ix,2]<holm_seq
#determine the original order
aux2<-sort(aux$ix,index.return=T)
#print the results
data.frame(comp=res[,1],est=round(res[,3],digits=2),p_val=round(res[,2],digits=3),signif=signif[aux2$ix])


############Supplementary table 6################

#compare convalescent and others
table(dat1$impf)
all_groups<-c("M/M","P/P","A/P","A/A","J") #only one observation in A/M

feature<-"RBD" #or "WT"

res<-data.frame(comp=all_groups,res=rep(NA,length(all_groups)),
                stat=rep(NA,length(all_groups)))
for(i in 1: length(all_groups)){
  aux<-lawstat::brunner.munzel.test(dat[dat$impf==all_groups[i],feature],
                                    dat1[dat1$impf==all_groups[i],feature])
  res[i,2]<-aux$p.value
  res[i,3]<-aux$estimate
}

#print the results
data.frame(comp=res[,1],est=round(res[,3],digits=2),p_val=round(res[,2],digits=3))
