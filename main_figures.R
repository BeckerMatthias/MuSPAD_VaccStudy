#Author: Matthias Becker 
#Affiliation: NMI Natural and Medical Sciences Institute at the University of Tuebingen, 
#             Reutlingen, Germany
#Contact: matthias.becker@nmi.de
#Date: 2022-01-13

#This code generates the figures for the publication:
#Comparative magnitude and persistence of humoral SARS-CoV-2 vaccination responses in the adult population in Germany
#Rstudio Version 1.2.5001 was used with R version 3.6.1
#Currently this code displays all figures in RStudio. If "#" is removed from all "#svg(...)" and "#dev.off()",
#the code will export files in the correct height and width as .svg files
#To generate final figures for the publication from the .svg files, InkScape vector graphics editor was used


##### Library and functions ####

library(beeswarm) # Used for stripchart depiction of data
library(gplots)   # Used for col2hex() function, enables use of named colors in custom boxplot


##Custom Boxplot/Stripchart overlay function
BoxMaB <- function(data, #List containing all sets to be displayed in sequential order (left to right)
                   MWU = F, #Either False to disable stat calculation or list of vectors of length 3 for all comparisons to be made
                   #Elements 1 and 2 of the vector give the data sets to be compared 
                   #Element 3 gives the level at which to display the stats in the graph
                   spacing = c(1:length(data)), #Vector, for which length must be = # of sets denoting the spacing on the x-axis
                   color = "gray30", # Vector of colors to use
                   logscale = F,#Whether to display data on log scale
                   main = "",#Title to be displayed above plot
                   ylab = "",#ylab for plot, xlab has to be handeled separately by user
                   points = T,#whether to draw the stripchart
                   ptcex = 0.5,#size of points in stripchart
                   ptcorral = "none",#Beeswarm corralling of points
                   yaxticks = T)#if T, will use default y-axis intervals for non log scale and c(0.01,0.1,1,10,100) for log scale, 
                                #alternatively give vector where to put y axis ticks
{
  require(beeswarm)
  require(gplots)
  
  #log transform data for display if required
  if(logscale==T){for(i in (1:length(data))){data[[i]] <- log10(data[[i]])}}
  
  #set parameters for stats display and set y-axis limit
  if (class(MWU)=="logical"){if (MWU==F){ylim=range(data)}}else{
    MWUlevels <- unlist(MWU)[seq(from=3,to=length(MWU)*3,by=3)]
    sbase=range(data)[2]
    sincr=(range(data)[2]-range(data)[1])*0.1
    ylim=c(range(data)[1],sbase+sincr*max(MWUlevels))
  }
  
  #Draw box backgrounds and set graph limits
  if(points==T){
    boxplot(data,
            at=spacing,
            xlim=c(range(spacing)[1]-0.5,range(spacing)[2]+0.5),
            ylim=ylim,
            outline=F,
            col=paste(col2hex(color),"30",sep=""),
            lwd=0.75,border="#00000000",xaxt="n",yaxt="n",
            main=main,
            ylab=ylab)
  }else if(points==F){
    boxplot(data,
            at=spacing,
            xlim=c(range(spacing)[1]-0.5,range(spacing)[2]+0.5),
            ylim=ylim,
            outline=T,
            col=color,
            xaxt="n",yaxt="n",
            main=main,
            ylab=ylab)
  }
  #Draw points for stripchart
  if(points==T){
    beeswarm(data,
             at=spacing,
             pch=21,col="black",bg=color,
             cex=ptcex,
             corral=ptcorral, 
             xaxt="n",
             add=T)
  }
  
  #Draw boxes
  if(points==T){
    boxplot(data,at=spacing,outline=F,lwd=0.75,col="#00000000",add=T,yaxt="n",xaxt="n")
  }
  
  #Draw y-axis
  if (class(yaxticks)=="logical"){
    if (yaxticks==T&logscale==F) {axis(2)} else if (yaxticks==T&logscale==T) {
      axis(2,at=log10(c(0.01,0.1,1,10,100,1000)),labels=c("0.01","0.1","1","10","100","1,000"))} else if (yaxticks==F){}
  }else if (logscale==F) {axis(2,at=yaxticks)} else if (logscale==T) {
    axis(2,at=log10(yaxticks),labels=format(yaxticks,big.mark=",",trim=T))}
  
  #Perform MWU statistics and display them
  if (class(MWU)=="logical"){if (MWU==F){}}else{
    for (i in c(1:length(MWU))){
      if(length(data[[MWU[[i]][1]]])==0|length(data[[MWU[[i]][2]]])==0){}else{
        lines(c(rep(spacing[MWU[[i]][1]],2),rep(spacing[MWU[[i]][2]],2)),
              c(sbase+sincr*MWUlevels[i]-sincr*0.5,rep(sbase+sincr*MWUlevels[i]-sincr*0.25,2),sbase+sincr*MWUlevels[i]-sincr*0.5))
        text(mean(c(spacing[MWU[[i]][1]],spacing[MWU[[i]][2]])),sbase+sincr*MWUlevels[i]+sincr*0.1,
             paste(formatC(wilcox.test(data[[MWU[[i]][1]]],data[[MWU[[i]][2]]],alternative="two.sided")$p.value,
                           format="e",digits=2)),
             cex = 1)
      }
    }
  }
}


#####
##### Data Read in and setup ######


##Read Data from .csv files
MM <- read.csv("data_mix_and_match_shared.csv",h=T,sep=",")
TP <- read.csv("data_timepoints_shared.csv",h=T,sep=",")
LT <- read.csv("data_longitudinal_shared.csv",h=T,sep=",")

##For Mix/Match Group: Remove donors with reported positive test or Nucleocapsid (N) signal >1, 
##save them in a separate data frame separately
MM_inf <-  MM[!(MM$N<=1&MM$test_reported_positiv==0),]
MM <- MM[MM$N<=1&MM$test_reported_positiv==0,]

## Color scheme ###
mycol <- c("#0086cb","#cbab00","#00cb12","#00cbbb","#96cb00","#727272") #blue to yellow via green

#####
##### Figure 1  #####

AGs <- c("Spike","RBD","S1","S2")
labels <- c(expression("Normalised Spike IgG Signal"),
            expression("Normalised RBD"[wt]~"IgG Signal"),
            expression("Normalised S1 IgG Signal"),
            expression("Normalised S2 IgG Signal"))
for (AG in AGs){
  #svg(paste("Figure1_",AG,".svg",sep=""),7.5,5)
  par(mfrow=c(1,1),mar=c(3,4,2,2),pty="m")
  BoxMaB(data=list(MM[MM$vacc_type=="Moderna-Moderna",AG],
                   MM[MM$vacc_type=="Pfizer-Pfizer",AG],
                   MM[MM$vacc_type=="Astra-Astra",AG],
                   MM[MM$vacc_type=="Astra-Moderna",AG],
                   MM[MM$vacc_type=="Astra-Pfizer",AG],
                   MM[MM$vacc_type=="Janssen",AG]),
         color=mycol,
         MWU=F,
         spacing=c(1:6),
         logscale = F,
         ptcex=0.5,
         ptcorral="wrap")
  axis(1,1:6,labels=c("M/M","P/P","A/A","A/M","A/P","J"))
  mtext(text=paste("n=",c(length(MM[MM$vacc_type=="Moderna-Moderna",AG]),
                          length(MM[MM$vacc_type=="Pfizer-Pfizer",AG]),
                          length(MM[MM$vacc_type=="Astra-Astra",AG]),
                          length(MM[MM$vacc_type=="Astra-Moderna",AG]),
                          length(MM[MM$vacc_type=="Astra-Pfizer",AG]),
                          length(MM[MM$vacc_type=="Janssen",AG])),
                   sep=""),
        side=1,at=1:6,line=2)
  mtext(text=c(round(mean(MM[MM$vacc_type=="Moderna-Moderna","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Pfizer-Pfizer","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Astra-Astra","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Astra-Moderna","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Astra-Pfizer","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Janssen","time_since_full_vacc"]),1)),
        side=3,at=1:6,line=1)
  mtext(text=paste("sd=",c(round(sd(MM[MM$vacc_type=="Moderna-Moderna","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_type=="Pfizer-Pfizer","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_type=="Astra-Astra","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_type=="Astra-Moderna","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_type=="Astra-Pfizer","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_first_type=="Janssen","time_since_full_vacc"]),1)),
                   sep=""),
        side=3,at=1:6,line=0)
  mtext(expression("mean"~Delta~"t"[2]),side=3,line=0.5,at=0)
  mtext(text=labels[which(AGs==AG)],2,line=2)
  #dev.off()
}
rm(AGs,AG,labels)

#####
##### Figure 2  ####

#svg(paste("Figure2.svg",sep=""),5,14)

##Data (Using all antigens and all Samples)
hmdata <- MM
hmdata$Vacc_short <- ""
hmdata[hmdata$vacc_type=="Moderna-Moderna","Vacc_short"] <- "M/M"
hmdata[hmdata$vacc_type=="Pfizer-Pfizer","Vacc_short"] <- "P/P"
hmdata[hmdata$vacc_type=="Astra-Astra","Vacc_short"] <- "A/A"
hmdata[hmdata$vacc_type=="Astra-Moderna","Vacc_short"] <- "A/M"
hmdata[hmdata$vacc_type=="Astra-Pfizer","Vacc_short"] <- "A/P"
hmdata[hmdata$vacc_type=="Janssen","Vacc_short"] <- "J"

##Color
hmcol <- colorRampPalette(c("#cccccc","#888888","#000000","#cc0000","#ff0000"))(50)#Grays into Reds, stretched

##Select/Transform Data
hmdata <- hmdata[order(hmdata$Vacc_short),]#sort after vaccination
z <- hmdata[,c("Spike","RBD","S1","S2")]
z[,c("Spike","RBD","S1","S2")] <- log(z[,c("Spike","RBD","S1","S2")])#log transform data
z <- scale(z)
z[z>2.5] <- 2.5 #remove peaks
z[z<(-2.5)] <- -2.5 #remove peaks


##Clustering
o <- vector(length=nrow(hmdata)) # vector in which order will be collected
c <- 1 # counter to help with writing into vector o during loop
for (v in unique(hmdata$Vacc_short)){#loop for each group
  hc <- hclust(dist(z[hmdata$Vacc_short==v,]),method="centroid") # hclust for each subgroup separately
  # get order of subgroup data after clustering, but with indexes from the large table hmdata
  o[c:(c+nrow(hmdata[hmdata$Vacc_short==v,])-1)] <- which(hmdata$Vacc_short==v)[c(rev(hc$order))] 
  c <- c+nrow(hmdata[hmdata$Vacc_short==v,])
}
z <- z[o,]#re-order according to separate clusterings
hmdata <- hmdata[o,]

##Insert NAs into z and hmdata to create group searation visually (because using "rowsep" leads to data loss)
nNAs <- 20
tmpz <- list()#split and re-bind is easiest
c <- 0
for (i in unique(hmdata$Vacc_short)){
  c <- c+1
  tmpz[[i]] <- z[hmdata$Vacc_short==i,]
  tmpz[[2*c]] <- matrix(nrow=nNAs,ncol=ncol(z))
}
z <- do.call(rbind,tmpz)#re-bind

##Get adjusted row labels
hmnames <- rep("",nrow(z))
hmnames[which(complete.cases(z))] <- hmdata$Vacc_short

##Draw Heatmap
heatmap.2(as.matrix(z),main="",Colv=F,Rowv=F,dendrogram='none',
          labRow=hmnames,labCol=colnames(z),
          scale='none',col=hmcol,na.color="white",trace='none', margins=c(6,4),
          breaks=51,symkey=F,
          cexRow = 0.5, cexCol = 1.5)


#dev.off()
#Rownames and legend are further edited in InkScape for polishing
rm(hmdata,o,v,c,i,hc,z,hmcol,nNAs,tmpz,hmnames)


#####
##### Figure 3  ####

AGs <- c("WT")
labels <- c(expression("ACE2 binding inhibition (%)"))
for (AG in AGs){
  #svg(paste("Figure3.svg",sep=""),7.5,5)
  par(mfrow=c(1,1),mar=c(3,4,2,2),pty="m")
  BoxMaB(data=list(MM[MM$vacc_type=="Moderna-Moderna",AG],
                   MM[MM$vacc_type=="Pfizer-Pfizer",AG],
                   MM[MM$vacc_type=="Astra-Astra",AG],
                   MM[MM$vacc_type=="Astra-Moderna",AG],
                   MM[MM$vacc_type=="Astra-Pfizer",AG],
                   MM[MM$vacc_type=="Janssen",AG]),
         color=mycol,
         MWU=F,
         spacing=c(1:6),
         logscale = F,
         ptcex=0.5,
         ptcorral="wrap")
  axis(1,1:6,labels=c("M/M","P/P","A/A","A/M","A/P","J"))
  mtext(text=paste("n=",c(length(MM[MM$vacc_type=="Moderna-Moderna",AG]),
                          length(MM[MM$vacc_type=="Pfizer-Pfizer",AG]),
                          length(MM[MM$vacc_type=="Astra-Astra",AG]),
                          length(MM[MM$vacc_type=="Astra-Moderna",AG]),
                          length(MM[MM$vacc_type=="Astra-Pfizer",AG]),
                          length(MM[MM$vacc_type=="Janssen",AG])),
                   sep=""),
        side=1,at=1:6,line=2)
  mtext(text=c(round(mean(MM[MM$vacc_type=="Moderna-Moderna","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Pfizer-Pfizer","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Astra-Astra","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Astra-Moderna","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Astra-Pfizer","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Janssen","time_since_full_vacc"]),1)),
        side=3,at=1:6,line=1)
  mtext(text=paste("sd=",c(round(sd(MM[MM$vacc_type=="Moderna-Moderna","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_type=="Pfizer-Pfizer","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_type=="Astra-Astra","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_type=="Astra-Moderna","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_type=="Astra-Pfizer","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_first_type=="Janssen","time_since_full_vacc"]),1)),
                   sep=""),
        side=3,at=1:6,line=0)
  mtext(expression("mean"~Delta~"t"[2]),side=3,line=0.5,at=0)
  mtext(text=labels[which(AGs==AG)],2,line=2)
  #dev.off()
}
rm(AGs,AG,labels)

#####
##### Figure 4  ####


AGs <- c("WT")
labels <- c(expression("ACE2 binding inhibition (%)"))
for (AG in AGs){
  ##a - Split by delta T
  day <- 28 #which day post vaccination tp split groups at
  #svg(paste("Figure4a.svg",sep=""),9,6)
  par(mfrow=c(1,1),mar=c(3,4,2,2),pty="m")
  BoxMaB(data=list(MM[MM$vacc_type=="Moderna-Moderna"&MM$time_since_full_vacc<day,AG],
                   MM[MM$vacc_type=="Moderna-Moderna"&MM$time_since_full_vacc>=day,AG],
                   MM[MM$vacc_type=="Pfizer-Pfizer"&MM$time_since_full_vacc<day,AG],
                   MM[MM$vacc_type=="Pfizer-Pfizer"&MM$time_since_full_vacc>=day,AG],
                   MM[MM$vacc_type=="Astra-Astra"&MM$time_since_full_vacc<day,AG],
                   MM[MM$vacc_type=="Astra-Astra"&MM$time_since_full_vacc>=day,AG],
                   MM[MM$vacc_type=="Astra-Moderna"&MM$time_since_full_vacc<day,AG],
                   MM[MM$vacc_type=="Astra-Moderna"&MM$time_since_full_vacc>=day,AG],
                   MM[MM$vacc_type=="Astra-Pfizer"&MM$time_since_full_vacc<day,AG],
                   MM[MM$vacc_type=="Astra-Pfizer"&MM$time_since_full_vacc>=day,AG],
                   MM[MM$vacc_type=="Janssen"&MM$time_since_full_vacc<day,AG],
                   MM[MM$vacc_type=="Janssen"&MM$time_since_full_vacc>=day,AG]),
         color=rep(mycol,each=2),
         spacing=c(1,2,3.5,4.5,6,7,8.5,9.5,11,12,13.5,14.5),
         logscale = F,
         ptcex=0.5,
         ptcorral="wrap")
  axis(1,c(1,2,3.5,4.5,6,7,8.5,9.5,11,12,13.5,14.5),labels=rep(paste(c("<",">="),day,sep=""),6))
  mtext(side=1,at=-0.5,text=expression(Delta~t[2]~"(days)"),line=1)
  mtext(side=1,at=c(1.5,4,6.5,9,11.5,14),text=c("M/M","P/P","A/A","A/M","A/P","J"),line=2)
  mtext(text=c("n =",
               length(MM[MM$vacc_type=="Moderna-Moderna"&MM$time_since_full_vacc<day,AG]),
               length(MM[MM$vacc_type=="Moderna-Moderna"&MM$time_since_full_vacc>=day,AG]),
               length(MM[MM$vacc_type=="Pfizer-Pfizer"&MM$time_since_full_vacc<day,AG]),
               length(MM[MM$vacc_type=="Pfizer-Pfizer"&MM$time_since_full_vacc>=day,AG]),
               length(MM[MM$vacc_type=="Astra-Astra"&MM$time_since_full_vacc<day,AG]),
               length(MM[MM$vacc_type=="Astra-Astra"&MM$time_since_full_vacc>=day,AG]),
               length(MM[MM$vacc_type=="Astra-Moderna"&MM$time_since_full_vacc<day,AG]),
               length(MM[MM$vacc_type=="Astra-Moderna"&MM$time_since_full_vacc>=day,AG]),
               length(MM[MM$vacc_type=="Astra-Pfizer"&MM$time_since_full_vacc<day,AG]),
               length(MM[MM$vacc_type=="Astra-Pfizer"&MM$time_since_full_vacc>=day,AG]),
               length(MM[MM$vacc_type=="Janssen"&MM$time_since_full_vacc<day,AG]),
               length(MM[MM$vacc_type=="Janssen"&MM$time_since_full_vacc>=day,AG])),
        side=3,at=c(0,1,2,3.5,4.5,6,7,8.5,9.5,11,12,13.5,14.5),line=0.5)
  mtext(text=labels[which(AGs==AG)],2,line=2)
  #dev.off()
  
  ##b - Split by sex
  #svg(paste("Figure4b.svg",sep=""),9,6)
  par(mfrow=c(1,1),mar=c(3,4,2,2),pty="m")
  BoxMaB(data=list(MM[MM$vacc_type=="Moderna-Moderna"&MM$sex=="female",AG],
                   MM[MM$vacc_type=="Moderna-Moderna"&MM$sex=="male",AG],
                   MM[MM$vacc_type=="Pfizer-Pfizer"&MM$sex=="female",AG],
                   MM[MM$vacc_type=="Pfizer-Pfizer"&MM$sex=="male",AG],
                   MM[MM$vacc_type=="Astra-Astra"&MM$sex=="female",AG],
                   MM[MM$vacc_type=="Astra-Astra"&MM$sex=="male",AG],
                   MM[MM$vacc_type=="Astra-Moderna"&MM$sex=="female",AG],
                   MM[MM$vacc_type=="Astra-Moderna"&MM$sex=="male",AG],
                   MM[MM$vacc_type=="Astra-Pfizer"&MM$sex=="female",AG],
                   MM[MM$vacc_type=="Astra-Pfizer"&MM$sex=="male",AG],
                   MM[MM$vacc_type=="Janssen"&MM$sex=="female",AG],
                   MM[MM$vacc_type=="Janssen"&MM$sex=="male",AG]),
         color=rep(mycol,each=2),
         spacing=c(1,2,3.5,4.5,6,7,8.5,9.5,11,12,13.5,14.5),
         logscale = F,
         ptcex=0.5,
         ptcorral="wrap")
  axis(1,c(1,2,3.5,4.5,6,7,8.5,9.5,11,12,13.5,14.5),labels=rep(c("f","m"),6))
  mtext(side=1,at=-0.5,text="sex",line=1)
  mtext(side=1,at=c(1.5,4,6.5,9,11.5,14),text=c("M/M","P/P","A/A","A/M","A/P","J"),line=2)
  mtext(text=c("n =",
               length(MM[MM$vacc_type=="Moderna-Moderna"&MM$sex=="female",AG]),
               length(MM[MM$vacc_type=="Moderna-Moderna"&MM$sex=="male",AG]),
               length(MM[MM$vacc_type=="Pfizer-Pfizer"&MM$sex=="female",AG]),
               length(MM[MM$vacc_type=="Pfizer-Pfizer"&MM$sex=="male",AG]),
               length(MM[MM$vacc_type=="Astra-Astra"&MM$sex=="female",AG]),
               length(MM[MM$vacc_type=="Astra-Astra"&MM$sex=="male",AG]),
               length(MM[MM$vacc_type=="Astra-Moderna"&MM$sex=="female",AG]),
               length(MM[MM$vacc_type=="Astra-Moderna"&MM$sex=="male",AG]),
               length(MM[MM$vacc_type=="Astra-Pfizer"&MM$sex=="female",AG]),
               length(MM[MM$vacc_type=="Astra-Pfizer"&MM$sex=="male",AG]),
               length(MM[MM$vacc_type=="Janssen"&MM$sex=="female",AG]),
               length(MM[MM$vacc_type=="Janssen"&MM$sex=="male",AG])),
        side=3,at=c(0,1,2,3.5,4.5,6,7,8.5,9.5,11,12,13.5,14.5),line=0.5)
  mtext(text=labels[which(AGs==AG)],2,line=2)
  #dev.off()
  
  ##c - split by age group
  #svg(paste("Figure4c.svg",sep=""),16,6)
  par(mfrow=c(1,1),mar=c(4,4,2,2),pty="m")
  BoxMaB(data=list(MM[MM$vacc_type=="Moderna-Moderna"&MM$age_group=="18-25",AG],
                   MM[MM$vacc_type=="Moderna-Moderna"&MM$age_group=="26-45",AG],
                   MM[MM$vacc_type=="Moderna-Moderna"&MM$age_group=="46-65",AG],
                   MM[MM$vacc_type=="Moderna-Moderna"&MM$age_group=="66-79",AG],
                   MM[MM$vacc_type=="Moderna-Moderna"&MM$age_group==">79",AG],
                   MM[MM$vacc_type=="Pfizer-Pfizer"&MM$age_group=="18-25",AG],
                   MM[MM$vacc_type=="Pfizer-Pfizer"&MM$age_group=="26-45",AG],
                   MM[MM$vacc_type=="Pfizer-Pfizer"&MM$age_group=="46-65",AG],
                   MM[MM$vacc_type=="Pfizer-Pfizer"&MM$age_group=="66-79",AG],
                   MM[MM$vacc_type=="Pfizer-Pfizer"&MM$age_group==">79",AG],
                   MM[MM$vacc_type=="Astra-Astra"&MM$age_group=="18-25",AG],
                   MM[MM$vacc_type=="Astra-Astra"&MM$age_group=="26-45",AG],
                   MM[MM$vacc_type=="Astra-Astra"&MM$age_group=="46-65",AG],
                   MM[MM$vacc_type=="Astra-Astra"&MM$age_group=="66-79",AG],
                   MM[MM$vacc_type=="Astra-Astra"&MM$age_group==">79",AG],
                   MM[MM$vacc_type=="Astra-Moderna"&MM$age_group=="18-25",AG],
                   MM[MM$vacc_type=="Astra-Moderna"&MM$age_group=="26-45",AG],
                   MM[MM$vacc_type=="Astra-Moderna"&MM$age_group=="46-65",AG],
                   MM[MM$vacc_type=="Astra-Moderna"&MM$age_group=="66-79",AG],
                   MM[MM$vacc_type=="Astra-Moderna"&MM$age_group==">79",AG],
                   MM[MM$vacc_type=="Astra-Pfizer"&MM$age_group=="18-25",AG],
                   MM[MM$vacc_type=="Astra-Pfizer"&MM$age_group=="26-45",AG],
                   MM[MM$vacc_type=="Astra-Pfizer"&MM$age_group=="46-65",AG],
                   MM[MM$vacc_type=="Astra-Pfizer"&MM$age_group=="66-79",AG],
                   MM[MM$vacc_type=="Astra-Pfizer"&MM$age_group==">79",AG],
                   MM[MM$vacc_type=="Janssen"&MM$age_group=="18-25",AG],
                   MM[MM$vacc_type=="Janssen"&MM$age_group=="26-45",AG],
                   MM[MM$vacc_type=="Janssen"&MM$age_group=="46-65",AG],
                   MM[MM$vacc_type=="Janssen"&MM$age_group=="66-79",AG],
                   MM[MM$vacc_type=="Janssen"&MM$age_group==">79",AG]),
         color=rep(mycol,each=5),
         spacing=c(1:5,6.5:10.5,12:16,17.5:21.5,23:27,28.5:32.5),
         logscale = F,
         ptcex=0.5,
         ptcorral="wrap")
  axis(1,c(1:5,6.5:10.5,12:16,17.5:21.5,23:27,28.5:32.5),labels=F)
  mtext(side=1,line=1.5,at=c(1,3,5,6.5,8.5,10.5,12,14,16,17.5,19.5,21.5,23,25,27,28.5,30.5,32.5),
        text=rep(c("18-25","46-65",">79"),6))
  mtext(side=1,line=0.5,at=c(2,4,7.5,9.5,13,15,18.5,20.5,24,26,29.5,31.5),
        text=rep(c("26-45","66-79"),6))
  mtext(side=1,at=-2,text="age\ngroup",line=1.5)
  mtext(side=1,at=c(3,8.5,14,19.5,25,30.5),text=c("M/M","P/P","A/A","A/M","A/P","J"),line=2.5)
  mtext(text=c("n =",
               length(MM[MM$vacc_type=="Moderna-Moderna"&MM$age_group=="18-25",AG]),
               length(MM[MM$vacc_type=="Moderna-Moderna"&MM$age_group=="26-45",AG]),
               length(MM[MM$vacc_type=="Moderna-Moderna"&MM$age_group=="46-65",AG]),
               length(MM[MM$vacc_type=="Moderna-Moderna"&MM$age_group=="66-79",AG]),
               length(MM[MM$vacc_type=="Moderna-Moderna"&MM$age_group==">79",AG]),
               length(MM[MM$vacc_type=="Pfizer-Pfizer"&MM$age_group=="18-25",AG]),
               length(MM[MM$vacc_type=="Pfizer-Pfizer"&MM$age_group=="26-45",AG]),
               length(MM[MM$vacc_type=="Pfizer-Pfizer"&MM$age_group=="46-65",AG]),
               length(MM[MM$vacc_type=="Pfizer-Pfizer"&MM$age_group=="66-79",AG]),
               length(MM[MM$vacc_type=="Pfizer-Pfizer"&MM$age_group==">79",AG]),
               length(MM[MM$vacc_type=="Astra-Astra"&MM$age_group=="18-25",AG]),
               length(MM[MM$vacc_type=="Astra-Astra"&MM$age_group=="26-45",AG]),
               length(MM[MM$vacc_type=="Astra-Astra"&MM$age_group=="46-65",AG]),
               length(MM[MM$vacc_type=="Astra-Astra"&MM$age_group=="66-79",AG]),
               length(MM[MM$vacc_type=="Astra-Astra"&MM$age_group==">79",AG]),
               length(MM[MM$vacc_type=="Astra-Moderna"&MM$age_group=="18-25",AG]),
               length(MM[MM$vacc_type=="Astra-Moderna"&MM$age_group=="26-45",AG]),
               length(MM[MM$vacc_type=="Astra-Moderna"&MM$age_group=="46-65",AG]),
               length(MM[MM$vacc_type=="Astra-Moderna"&MM$age_group=="66-79",AG]),
               length(MM[MM$vacc_type=="Astra-Moderna"&MM$age_group==">79",AG]),
               length(MM[MM$vacc_type=="Astra-Pfizer"&MM$age_group=="18-25",AG]),
               length(MM[MM$vacc_type=="Astra-Pfizer"&MM$age_group=="26-45",AG]),
               length(MM[MM$vacc_type=="Astra-Pfizer"&MM$age_group=="46-65",AG]),
               length(MM[MM$vacc_type=="Astra-Pfizer"&MM$age_group=="66-79",AG]),
               length(MM[MM$vacc_type=="Astra-Pfizer"&MM$age_group==">79",AG]),
               length(MM[MM$vacc_type=="Janssen"&MM$age_group=="18-25",AG]),
               length(MM[MM$vacc_type=="Janssen"&MM$age_group=="26-45",AG]),
               length(MM[MM$vacc_type=="Janssen"&MM$age_group=="46-65",AG]),
               length(MM[MM$vacc_type=="Janssen"&MM$age_group=="66-79",AG]),
               length(MM[MM$vacc_type=="Janssen"&MM$age_group==">79",AG])),
        side=3,at=c(-0.5,1:5,6.5:10.5,12:16,17.5:21.5,23:27,28.5:32.5),line=0.5)
  mtext(text=labels[which(AGs==AG)],2,line=2)
  #dev.off()
  
}
rm(AGs,AG,labels,day)


#####
##### Figure 5  ####

AGs <- c("WT")
labels <- c(expression("ACE2 binding inhibition (%)"))
for (AG in AGs){
  #svg(paste("Figure5a.svg",sep=""),7.5,5)
  par(mfrow=c(1,1),mar=c(3,4,2,2),pty="m")
  BoxMaB(data=list(MM_inf[MM_inf$vacc_type=="Moderna-Moderna",AG],
                   MM_inf[MM_inf$vacc_type=="Pfizer-Pfizer",AG],
                   MM_inf[MM_inf$vacc_type=="Astra-Astra",AG],
                   MM_inf[MM_inf$vacc_type=="Astra-Moderna",AG],
                   MM_inf[MM_inf$vacc_type=="Astra-Pfizer",AG],
                   MM_inf[MM_inf$vacc_type=="Janssen",AG]),
         color=mycol,
         MWU=F,
         spacing=c(1:6),
         logscale = F,
         ptcex=0.5,
         ptcorral="wrap")
  axis(1,1:6,labels=c("M/M","P/P","A/A","A/M","A/P","J"))
  mtext(text=paste("n=",c(length(MM_inf[MM_inf$vacc_type=="Moderna-Moderna",AG]),
                          length(MM_inf[MM_inf$vacc_type=="Pfizer-Pfizer",AG]),
                          length(MM_inf[MM_inf$vacc_type=="Astra-Astra",AG]),
                          length(MM_inf[MM_inf$vacc_type=="Astra-Moderna",AG]),
                          length(MM_inf[MM_inf$vacc_type=="Astra-Pfizer",AG]),
                          length(MM_inf[MM_inf$vacc_type=="Janssen",AG])),
                   sep=""),
        side=1,at=1:6,line=2)
  mtext(text=c(round(mean(MM_inf[MM_inf$vacc_type=="Moderna-Moderna","time_since_full_vacc"]),1),
               round(mean(MM_inf[MM_inf$vacc_type=="Pfizer-Pfizer","time_since_full_vacc"]),1),
               round(mean(MM_inf[MM_inf$vacc_type=="Astra-Astra","time_since_full_vacc"]),1),
               round(mean(MM_inf[MM_inf$vacc_type=="Astra-Moderna","time_since_full_vacc"]),1),
               round(mean(MM_inf[MM_inf$vacc_type=="Astra-Pfizer","time_since_full_vacc"]),1),
               round(mean(MM_inf[MM_inf$vacc_type=="Janssen","time_since_full_vacc"]),1)),
        side=3,at=1:6,line=1)
  mtext(text=paste("sd=",c(round(sd(MM_inf[MM_inf$vacc_type=="Moderna-Moderna","time_since_full_vacc"]),1),
                           round(sd(MM_inf[MM_inf$vacc_type=="Pfizer-Pfizer","time_since_full_vacc"]),1),
                           round(sd(MM_inf[MM_inf$vacc_type=="Astra-Astra","time_since_full_vacc"]),1),
                           round(sd(MM_inf[MM_inf$vacc_type=="Astra-Moderna","time_since_full_vacc"]),1),
                           round(sd(MM_inf[MM_inf$vacc_type=="Astra-Pfizer","time_since_full_vacc"]),1),
                           round(sd(MM_inf[MM_inf$vacc_first_type=="Janssen","time_since_full_vacc"]),1)),
                   sep=""),
        side=3,at=1:6,line=0)
  mtext(expression("mean"~Delta~"t"[2]),side=3,line=0.5,at=0)
  mtext(text=labels[which(AGs==AG)],2,line=2)
  #dev.off()
  #svg(paste("Figure5b.svg",sep=""),7.5,5)
  par(mfrow=c(1,1),mar=c(3,4,2,2),pty="m")
  BoxMaB(data=list(MM[MM$vacc_type=="Moderna-Moderna",AG],
                   MM[MM$vacc_type=="Pfizer-Pfizer",AG],
                   MM[MM$vacc_type=="Astra-Astra",AG],
                   MM[MM$vacc_type=="Astra-Moderna",AG],
                   MM[MM$vacc_type=="Astra-Pfizer",AG],
                   MM[MM$vacc_type=="Janssen",AG]),
         color=mycol,
         MWU=F,
         spacing=c(1:6),
         logscale = F,
         ptcex=0.5,
         ptcorral="wrap")
  axis(1,1:6,labels=c("M/M","P/P","A/A","A/M","A/P","J"))
  mtext(text=paste("n=",c(length(MM[MM$vacc_type=="Moderna-Moderna",AG]),
                          length(MM[MM$vacc_type=="Pfizer-Pfizer",AG]),
                          length(MM[MM$vacc_type=="Astra-Astra",AG]),
                          length(MM[MM$vacc_type=="Astra-Moderna",AG]),
                          length(MM[MM$vacc_type=="Astra-Pfizer",AG]),
                          length(MM[MM$vacc_type=="Janssen",AG])),
                   sep=""),
        side=1,at=1:6,line=2)
  mtext(text=c(round(mean(MM[MM$vacc_type=="Moderna-Moderna","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Pfizer-Pfizer","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Astra-Astra","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Astra-Moderna","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Astra-Pfizer","time_since_full_vacc"]),1),
               round(mean(MM[MM$vacc_type=="Janssen","time_since_full_vacc"]),1)),
        side=3,at=1:6,line=1)
  mtext(text=paste("sd=",c(round(sd(MM[MM$vacc_type=="Moderna-Moderna","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_type=="Pfizer-Pfizer","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_type=="Astra-Astra","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_type=="Astra-Moderna","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_type=="Astra-Pfizer","time_since_full_vacc"]),1),
                           round(sd(MM[MM$vacc_first_type=="Janssen","time_since_full_vacc"]),1)),
                   sep=""),
        side=3,at=1:6,line=0)
  mtext(expression("mean"~Delta~"t"[2]),side=3,line=0.5,at=0)
  mtext(text=labels[which(AGs==AG)],2,line=2)
  #dev.off()
}
rm(AGs,AG,labels)


#####
##### Figure 6  ####

AGs <- c("RBD","WT")
labels <- c(expression("Normalised RBD"[wt]~"IgG Signal"),
            expression("ACE2 RBD"[wt]~"binding inhibition (%)"))
for (AG in AGs){
  #The following groups exist for delta T post vaccination: 5-12,26-30,54-58,94-103,129-146,176-203
  #svg(paste("Figure6_",AG,".svg",sep=""),7.5,5)
  par(mfrow=c(1,1),mar=c(4,4,2,2),pty="m")
  plot(TP$time_since_full_vacc,TP[,AG],cex=0,log="",xaxt="n",
       xlab=expression(Delta~"t"["2"^"nd"~" vaccination"]~"(days)"),ylab=labels[which(AGs==AG)])
  points(TP[TP$vacc_type=="Pfizer-Pfizer","time_since_full_vacc"],
         TP[TP$vacc_type=="Pfizer-Pfizer",AG],pch=21,col="#00000080",bg=paste(mycol[2],80,sep=""))
  points(TP[TP$vacc_type=="Moderna-Moderna","time_since_full_vacc"],
         TP[TP$vacc_type=="Moderna-Moderna",AG],pch=21,col="#00000080",bg=paste(mycol[1],80,sep=""))
  if(AG=="RBD"){legend("topright",fill=mycol[1:2],legend=c("mRNA-1273","BNT162b2"),bty="n")}
  if(AG=="WT"){legend(160,0.9,fill=mycol[1:2],legend=c("mRNA-1273","BNT162b2"),bty="n")}
  axis(1,at=c(median(TP$time_since_full_vacc[TP$time_since_full_vacc>=4&TP$time_since_full_vacc<=12]),
              median(TP$time_since_full_vacc[TP$time_since_full_vacc>=26&TP$time_since_full_vacc<=30]),
              median(TP$time_since_full_vacc[TP$time_since_full_vacc>=54&TP$time_since_full_vacc<=58]),
              median(TP$time_since_full_vacc[TP$time_since_full_vacc>=94&TP$time_since_full_vacc<=103]),
              median(TP$time_since_full_vacc[TP$time_since_full_vacc>=129&TP$time_since_full_vacc<=146]),
              median(TP$time_since_full_vacc[TP$time_since_full_vacc>=176&TP$time_since_full_vacc<=203])),
       labels=c("5-12","26-30","54-58","94-103","129-146","176-203"))
  #Median line for Pfizer-Pfizer
  lines(c(median(TP$time_since_full_vacc[TP$time_since_full_vacc>=4&TP$time_since_full_vacc<=12]),
          median(TP$time_since_full_vacc[TP$time_since_full_vacc>=26&TP$time_since_full_vacc<=30]),
          median(TP$time_since_full_vacc[TP$time_since_full_vacc>=54&TP$time_since_full_vacc<=58]),
          median(TP$time_since_full_vacc[TP$time_since_full_vacc>=94&TP$time_since_full_vacc<=103]),
          median(TP$time_since_full_vacc[TP$time_since_full_vacc>=129&TP$time_since_full_vacc<=146]),
          median(TP$time_since_full_vacc[TP$time_since_full_vacc>=176&TP$time_since_full_vacc<=203])),
        c(median(TP[TP$time_since_full_vacc>=4&TP$time_since_full_vacc<=12&TP$vacc_type=="Moderna-Moderna",AG]),
          median(TP[TP$time_since_full_vacc>=26&TP$time_since_full_vacc<=30&TP$vacc_type=="Moderna-Moderna",AG]),
          median(TP[TP$time_since_full_vacc>=54&TP$time_since_full_vacc<=58&TP$vacc_type=="Moderna-Moderna",AG]),
          median(TP[TP$time_since_full_vacc>=94&TP$time_since_full_vacc<=103&TP$vacc_type=="Moderna-Moderna",AG]),
          median(TP[TP$time_since_full_vacc>=129&TP$time_since_full_vacc<=146&TP$vacc_type=="Moderna-Moderna",AG]),
          median(TP[TP$time_since_full_vacc>=176&TP$time_since_full_vacc<=203&TP$vacc_type=="Moderna-Moderna",AG])),
        lwd=2,col=mycol[1])
  #Median line for Moderna-Moderna
  lines(c(median(TP$time_since_full_vacc[TP$time_since_full_vacc>=4&TP$time_since_full_vacc<=12]),
          median(TP$time_since_full_vacc[TP$time_since_full_vacc>=26&TP$time_since_full_vacc<=30]),
          median(TP$time_since_full_vacc[TP$time_since_full_vacc>=54&TP$time_since_full_vacc<=58]),
          median(TP$time_since_full_vacc[TP$time_since_full_vacc>=94&TP$time_since_full_vacc<=103]),
          median(TP$time_since_full_vacc[TP$time_since_full_vacc>=129&TP$time_since_full_vacc<=146]),
          median(TP$time_since_full_vacc[TP$time_since_full_vacc>=176&TP$time_since_full_vacc<=203])),
        c(median(TP[TP$time_since_full_vacc>=4&TP$time_since_full_vacc<=12&TP$vacc_type=="Pfizer-Pfizer",AG]),
          median(TP[TP$time_since_full_vacc>=26&TP$time_since_full_vacc<=30&TP$vacc_type=="Pfizer-Pfizer",AG]),
          median(TP[TP$time_since_full_vacc>=54&TP$time_since_full_vacc<=58&TP$vacc_type=="Pfizer-Pfizer",AG]),
          median(TP[TP$time_since_full_vacc>=94&TP$time_since_full_vacc<=103&TP$vacc_type=="Pfizer-Pfizer",AG]),
          median(TP[TP$time_since_full_vacc>=129&TP$time_since_full_vacc<=146&TP$vacc_type=="Pfizer-Pfizer",AG]),
          median(TP[TP$time_since_full_vacc>=176&TP$time_since_full_vacc<=203&TP$vacc_type=="Pfizer-Pfizer",AG])),
        lwd=2,col=mycol[2])
  #dev.off()
}
rm(AG,AGs,labels)





#####
##### Figure 7  ####

##a and c
AGs <- c("RBD","WT")
labels <- c(expression("Normalised RBD"[wt]~"IgG Signal"),
            expression("ACE2 RBD"[wt]~"binding inhibition (%)"))
for (AG in AGs){
  #svg(paste("Figure7_",AG,".svg",sep=""),7.5,5)
  par(mfrow=c(1,1),mar=c(4,4,2,2),pty="m")
  plot(c(LT$time_since_full_vacc_1,LT$time_since_full_vacc_2),c(LT[,paste(AG,"_1",sep="")],LT[,paste(AG,"_2",sep="")]),
       cex=0,log="",xaxt="n",
       xlab=expression(Delta~"t"["2"^"nd"~" vaccination"]~"(days)"),ylab=labels[which(AGs==AG)])
  for (donor in c(1:nrow(LT))[LT$vacc_type=="Pfizer-Pfizer"]){
    lines(c(LT[donor,"time_since_full_vacc_1"],LT[donor,"time_since_full_vacc_2"]),
          c(LT[donor,paste(AG,"_1",sep="")],LT[donor,paste(AG,"_2",sep="")]),col="gray50")
  }
  points(c(LT[LT$vacc_type=="Pfizer-Pfizer","time_since_full_vacc_1"],LT[LT$vacc_type=="Pfizer-Pfizer","time_since_full_vacc_2"]),
         c(LT[LT$vacc_type=="Pfizer-Pfizer",paste(AG,"_1",sep="")],LT[LT$vacc_type=="Pfizer-Pfizer",paste(AG,"_2",sep="")]),
         pch=21,col="#000000FF",bg=paste(mycol[2],"FF",sep=""))
  axis(1,at=seq(0,200,by=20))
  #dev.off()
}
rm(donor,AGs,labels,AG)

##b and d

tmp <- LT[LT$vacc_type=="Pfizer-Pfizer",c("Spike_1","Spike_2","RBD_1","RBD_2","S1_1","S1_2","S2_1","S2_2",
                                          "WT_1","WT_2","ALPHA_1","ALPHA_2","BETA_1","BETA_2",
                                          "GAMMA_1","GAMMA_2","DELTA_1","DELTA_2")]

##log2FCs as Boxplots per Antigen for IgG
#svg(paste("Figure7b.svg",sep=""),6,5)
par(mfrow=c(1,1),mar=c(4,4,2,2),pty="m")
BoxMaB(data=list(log2(tmp[,"Spike_2"]/tmp[,"Spike_1"]),
                 log2(tmp[,"RBD_2"]/tmp[,"RBD_1"]),
                 log2(tmp[,"S1_2"]/tmp[,"S1_1"]),
                 log2(tmp[,"S2_2"]/tmp[,"S2_1"])),
       color=mycol[2],
       MWU=F,
       spacing=c(1:4),
       logscale = F,
       ptcex=0.7,
       ptcorral="wrap")
abline(h=0,lwd=2,col="#4D4D4D80")
mtext(side=2,line=2.5,text="log2FC of antigen IgG signal")
axis(1,1:4,labels=c("Spike","RBD","S1","S2"))
mtext(side=1,at=c(0.25,1:4),line=2,text=paste(c("mean =",
                                                round(mean(log2(tmp[,"Spike_2"]/tmp[,"Spike_1"])),2),
                                                round(mean(log2(tmp[,"RBD_2"]/tmp[,"RBD_1"])),2),
                                                round(mean(log2(tmp[,"S1_2"]/tmp[,"S1_1"])),2),
                                                round(mean(log2(tmp[,"S2_2"]/tmp[,"S2_1"])),2))))
mtext(side=1,at=c(0.25,1:4),line=3,text=paste(c("sd =",
                                                round(sd(log2(tmp[,"Spike_2"]/tmp[,"Spike_1"])),2),
                                                round(sd(log2(tmp[,"RBD_2"]/tmp[,"RBD_1"])),2),
                                                round(sd(log2(tmp[,"S1_2"]/tmp[,"S1_1"])),2),
                                                round(sd(log2(tmp[,"S2_2"]/tmp[,"S2_1"])),2))))
#dev.off()


##Differences as Boxplots per Antigen for ACE2
#svg(paste("Figure7d.svg",sep=""),6,5)
par(mfrow=c(1,1),mar=c(4,4,2,2),pty="m")
BoxMaB(data=list(tmp[,"WT_2"]-tmp[,"WT_1"],
                 tmp[,"ALPHA_2"]-tmp[,"ALPHA_1"],
                 tmp[,"BETA_2"]-tmp[,"BETA_1"],
                 tmp[,"GAMMA_2"]-tmp[,"GAMMA_1"],
                 tmp[,"DELTA_2"]-tmp[,"DELTA_1"]),
       color=mycol[2],
       MWU=F,
       spacing=c(1:5),
       logscale = F,
       ptcex=0.7,
       ptcorral="wrap")
abline(h=0,lwd=2,col="#4D4D4D80")
mtext(side=2,line=2.5,text="Difference in ACE2 RBD binding inhibition (%)")
axis(1,1:5,labels=F)
mtext(side=1,at=c(0.25,1:5),padj=0.5,line=0.5,
      text=c("RBD","wt",expression(alpha),expression(beta),expression(gamma),expression(delta)))
mtext(side=1,at=c(0.25,1:5),line=2,text=paste(c("mean =",
                                                round(mean(tmp[,"WT_2"]-tmp[,"WT_1"]),2),
                                                round(mean(tmp[,"ALPHA_2"]-tmp[,"ALPHA_1"]),2),
                                                round(mean(tmp[,"BETA_2"]-tmp[,"BETA_1"]),2),
                                                round(mean(tmp[,"GAMMA_2"]-tmp[,"GAMMA_1"]),2),
                                                round(mean(tmp[,"DELTA_2"]-tmp[,"DELTA_1"]),2))))
mtext(side=1,at=c(0.25,1:5),line=3,text=paste(c("sd =",
                                                round(sd(tmp[,"WT_2"]-tmp[,"WT_1"]),2),
                                                round(sd(tmp[,"ALPHA_2"]-tmp[,"ALPHA_1"]),2),
                                                round(sd(tmp[,"BETA_2"]-tmp[,"BETA_1"]),2),
                                                round(sd(tmp[,"GAMMA_2"]-tmp[,"GAMMA_1"]),2),
                                                round(sd(tmp[,"DELTA_2"]-tmp[,"DELTA_1"]),2))))
#dev.off()

rm(tmp)

#####
##### Supplementary Figure 1  ####


labels <- c(expression(wt),expression(alpha),expression(beta),
            expression(gamma),expression(delta))
#svg(paste("FigureS1.svg",sep=""),16,6)
par(mfrow=c(1,1),mar=c(4,4,2,2),pty="m")
BoxMaB(data=list(MM[MM$vacc_type=="Moderna-Moderna","WT"],
                 MM[MM$vacc_type=="Moderna-Moderna","ALPHA"],
                 MM[MM$vacc_type=="Moderna-Moderna","BETA"],
                 MM[MM$vacc_type=="Moderna-Moderna","GAMMA"],
                 MM[MM$vacc_type=="Moderna-Moderna","DELTA"],
                 MM[MM$vacc_type=="Pfizer-Pfizer","WT"],
                 MM[MM$vacc_type=="Pfizer-Pfizer","ALPHA"],
                 MM[MM$vacc_type=="Pfizer-Pfizer","BETA"],
                 MM[MM$vacc_type=="Pfizer-Pfizer","GAMMA"],
                 MM[MM$vacc_type=="Pfizer-Pfizer","DELTA"],
                 MM[MM$vacc_type=="Astra-Astra","WT"],
                 MM[MM$vacc_type=="Astra-Astra","ALPHA"],
                 MM[MM$vacc_type=="Astra-Astra","BETA"],
                 MM[MM$vacc_type=="Astra-Astra","GAMMA"],
                 MM[MM$vacc_type=="Astra-Astra","DELTA"],
                 MM[MM$vacc_type=="Astra-Moderna","WT"],
                 MM[MM$vacc_type=="Astra-Moderna","ALPHA"],
                 MM[MM$vacc_type=="Astra-Moderna","BETA"],
                 MM[MM$vacc_type=="Astra-Moderna","GAMMA"],
                 MM[MM$vacc_type=="Astra-Moderna","DELTA"],
                 MM[MM$vacc_type=="Astra-Pfizer","WT"],
                 MM[MM$vacc_type=="Astra-Pfizer","ALPHA"],
                 MM[MM$vacc_type=="Astra-Pfizer","BETA"],
                 MM[MM$vacc_type=="Astra-Pfizer","GAMMA"],
                 MM[MM$vacc_type=="Astra-Pfizer","DELTA"],
                 MM[MM$vacc_type=="Janssen","WT"],
                 MM[MM$vacc_type=="Janssen","ALPHA"],
                 MM[MM$vacc_type=="Janssen","BETA"],
                 MM[MM$vacc_type=="Janssen","GAMMA"],
                 MM[MM$vacc_type=="Janssen","DELTA"]),
       color=paste(rep(mycol,each=5),"80",sep=""),
       MWU=F,
       spacing=c(1:5,6.5:10.5,12:16,17.5:21.5,23:27,28.5:32.5),
       logscale = F,
       points = F)
axis(1,c(1:5,6.5:10.5,12:16,17.5:21.5,23:27,28.5:32.5),labels=F)
mtext(side=1,line=0.5,at=c(-1,1:5,6.5:10.5,12:16,17.5:21.5,23:27,28.5:32.5),
      text=c("RBD",rep(labels,6)),padj=0.5)
mtext(side=1,at=c(3,8.5,14,19.5,25,30.5),text=c("M/M","P/P","A/A","A/M","A/P","J"),line=2)
mtext(text=paste("n=",c(nrow(MM[MM$vacc_type=="Moderna-Moderna",]),
                        nrow(MM[MM$vacc_type=="Pfizer-Pfizer",]),
                        nrow(MM[MM$vacc_type=="Astra-Astra",]),
                        nrow(MM[MM$vacc_type=="Astra-Moderna",]),
                        nrow(MM[MM$vacc_type=="Astra-Pfizer",]),
                        nrow(MM[MM$vacc_type=="Janssen",])),
                 sep=""),
      side=1,at=c(3,8.5,14,19.5,25,30.5),line=3)
mtext(text="ACE2 binding inhibition(%)",2,line=2)
#dev.off()
rm(labels)









#####
##### Supplementary Figure 2  ####

AGs <- c("RBD")
labels <- c(expression("Normalised RBD"[wt]~"IgG Signal"))
for (AG in AGs){
  #svg(paste("FigureS2.svg",sep=""),7.5,5)
  par(mfrow=c(1,1),mar=c(3,4,2,2),pty="m")
  BoxMaB(data=list(MM_inf[MM_inf$vacc_type=="Moderna-Moderna",AG],
                   MM_inf[MM_inf$vacc_type=="Pfizer-Pfizer",AG],
                   MM_inf[MM_inf$vacc_type=="Astra-Astra",AG],
                   MM_inf[MM_inf$vacc_type=="Astra-Moderna",AG],
                   MM_inf[MM_inf$vacc_type=="Astra-Pfizer",AG],
                   MM_inf[MM_inf$vacc_type=="Janssen",AG]),
         color=mycol,
         MWU=F,
         spacing=c(1:6),
         logscale = F,
         ptcex=0.5,
         ptcorral="wrap")
  axis(1,1:6,labels=c("M/M","P/P","A/A","A/M","A/P","J"))
  mtext(text=paste("n=",c(length(MM_inf[MM_inf$vacc_type=="Moderna-Moderna",AG]),
                          length(MM_inf[MM_inf$vacc_type=="Pfizer-Pfizer",AG]),
                          length(MM_inf[MM_inf$vacc_type=="Astra-Astra",AG]),
                          length(MM_inf[MM_inf$vacc_type=="Astra-Moderna",AG]),
                          length(MM_inf[MM_inf$vacc_type=="Astra-Pfizer",AG]),
                          length(MM_inf[MM_inf$vacc_type=="Janssen",AG])),
                   sep=""),
        side=1,at=1:6,line=2)
  mtext(text=c(round(mean(MM_inf[MM_inf$vacc_type=="Moderna-Moderna","time_since_full_vacc"]),1),
               round(mean(MM_inf[MM_inf$vacc_type=="Pfizer-Pfizer","time_since_full_vacc"]),1),
               round(mean(MM_inf[MM_inf$vacc_type=="Astra-Astra","time_since_full_vacc"]),1),
               round(mean(MM_inf[MM_inf$vacc_type=="Astra-Moderna","time_since_full_vacc"]),1),
               round(mean(MM_inf[MM_inf$vacc_type=="Astra-Pfizer","time_since_full_vacc"]),1),
               round(mean(MM_inf[MM_inf$vacc_type=="Janssen","time_since_full_vacc"]),1)),
        side=3,at=1:6,line=1)
  mtext(text=paste("sd=",c(round(sd(MM_inf[MM_inf$vacc_type=="Moderna-Moderna","time_since_full_vacc"]),1),
                           round(sd(MM_inf[MM_inf$vacc_type=="Pfizer-Pfizer","time_since_full_vacc"]),1),
                           round(sd(MM_inf[MM_inf$vacc_type=="Astra-Astra","time_since_full_vacc"]),1),
                           round(sd(MM_inf[MM_inf$vacc_type=="Astra-Moderna","time_since_full_vacc"]),1),
                           round(sd(MM_inf[MM_inf$vacc_type=="Astra-Pfizer","time_since_full_vacc"]),1),
                           round(sd(MM_inf[MM_inf$vacc_first_type=="Janssen","time_since_full_vacc"]),1)),
                   sep=""),
        side=3,at=1:6,line=0)
  mtext(expression("mean"~Delta~"t"[2]),side=3,line=0.5,at=0)
  mtext(text=labels[which(AGs==AG)],2,line=2)
  #dev.off()
}
rm(AGs,AG,labels)






#####
