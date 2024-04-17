setwd("/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/B_Data")

#Load packages
library("lattice")
library("ggplot2")
library("dplyr")
library("readr")
library("rmarkdown")
library("devtools")
library("gghalves")
library("readxl")

# import plotrix package
library("plotrix")

#LMM specific
library("nlme")
library("ggeffects")
library("effects")

ntrls<-7

# Setup working dir and paths
FigOutDir <- ("/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/C_Plots")

# width and height variables for saved plots - leave it as default for now
w = 7
h = 4

#Read data
datfname <- paste("adv_ssst_newmaxf_fixICA_wfids_250_FullPowCondTable_ConsPats_allreps.csv")   
  
LFPdata_long  <- read.csv(datfname, sep=",")
LFPdata_long <- subset(LFPdata_long, double == 0)

LFPdata_long$ID <- as.factor(LFPdata_long$ID)
LFPdata_long$Condition <- as.factor(LFPdata_long$Condition)
LFPdata_long$Diagnosis <- as.factor(LFPdata_long$Diagnosis)

# Ideally it is already in long format

# Get group average across repetitions
gd <- LFPdata_long %>% 
  dplyr::group_by(Diagnosis, Condition) %>% 
  dplyr::summarise(TF_scond = mean(TF_scond))

# Put in CIs..
gd["CI_l"] <- NA
gd["CI_u"] <- NA

ind_c <- 0
for (i in 1:3) {
  for (j in 1:7){
    ind_c <- ind_c + 1
    
    tmp_data <- LFPdata_long %>% filter(Condition == j, Diagnosis == i) 
    
    gd[ind_c, "se"] <- std.error(tmp_data$TF_scond)
    
  }
}

#Pull out trial info
  
#Start jittering
  
set.seed(321)

xj <- jitter(as.numeric(LFPdata_long$Condition))
  
LFPdata_long$xj <- xj

posplotsep <- ggplot(data = gd, aes(x = as.numeric(Condition), y = TF_scond, color = Diagnosis)) +
  
geom_point(size = 5, aes(x = as.numeric(Condition), color = Diagnosis)) + #colour points by group

# geom_errorbar(aes(ymin = TF_scond - se, ymax = TF_scond + se)) +
  
geom_line(data = gd,  alpha = .8, size = 1) +
#geom_line(aes(group = ID, color = Diagnosis), alpha=0.5) +
  
#Note HR will be first (-1)    
scale_colour_manual(values = c("#E69F00", "#00BA38", "#3498db")) +
  
#geom_line(data=x_low, aes(x=Timebin_6, y=fit), size = 1, color='DarkOrange') +
  
#geom_line(data=x_high, aes(x=Timebin_6, y=fit), size = 1, color='DodgerBlue') +
  
ylab("Mean Beta1 Power (% change)") +
  
xlab("Repetition") + 
  
scale_x_continuous(breaks=c(1,2,3,4,5,6,7), labels=c("Dev", "rep1", "rep2",  "rep3", "rep4", "rep5", "rep6")) +
  
theme_classic() + 
theme(axis.title.y = element_text(face="bold", size = 24), axis.title.x = element_text(face="bold", size = 24)) + 
theme(axis.text.x = element_text(face="bold", size = 20, color = "black"), axis.text.y = element_text(face="bold", size = 20, color = "black")) +
  
#Turn caption off
theme(legend.position = "none")

figfname <- paste(FigOutDir, "TF_allreps_3grps_line.tiff", sep="/")

ggsave(figfname, w=12, h=8, dpi=300)

  #Start plot
  
  f1_box <- ggplot(my_sum_box_new, aes(x=Rep, y=TF_scond)) +
    
    geom_point(data = d, aes(x = xj, y = y), color = 'black', size = 2, alpha = .8) +
    geom_line(data = d, aes(x = xj, y = y, group = ID), color = 'black', alpha = 0.2) +
    scale_x_continuous(breaks=c(1,2,3,4,5,6,7), labels=c("Dev", "rep1", "rep2",  "rep3", "rep4", "rep5", "rep6")) +
    
    
    ylab("Mean Beta1 Power (% change)") +
    xlab("Repetition") + 
    theme_classic() +
    theme(axis.title.y = element_text(face="bold", size = 14), axis.title.x = element_text(face="bold", size = 15)) + 
    theme(axis.text.x = element_text(face="bold", size = 14, color = "black"), axis.text.y = element_text(face="bold", size = 14, color = "black")) +
    
    #Turn caption off
    theme(legend.position = "none")
  
  f1_box
  
  figfname <- paste(FigOutDir, '3groups_RepBar.tiff', sep="/")
  
  ggsave(figfname, width = w, height = h, dpi = 300)