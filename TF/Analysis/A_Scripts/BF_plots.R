
# Setup -------------------------------------------------------------------

library("BayesFactor")
library("bayestestR")
library(HDInterval)
library(ggridges)
library(tidyverse)
library(extrafont)
library("see")

# Setup working dir and paths

setwd("/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/B_Data")

FigOutDir <- ("/Users/alistairperry/Documents/Cambridge/Project_2/Fix_T1cor/TF/C_Plots")

LFPdata_wide <- read.csv("adv_ssst_newmaxf_fixICA_wfids_250_noTGBcon_FullPowCondTable_ConsPats.csv", header = TRUE)

CN <- subset(LFPdata_wide, Con_Pat == 1)


## Sample from the corresponding posterior distribution
samples = ttestBF(x = CN$RIFG_Std, y = CN$RIFG_Dev, paired=TRUE, posterior = TRUE, iterations = 10000)

# Plot effect size
# Samples[,3] for effect size
as_tibble(samples[,3]) %>%
  tidyr::pivot_longer(everything()) %>%
  ggplot(aes(x = value)) + 
  ggdist::stat_halfeye(.width = c(0.95),alpha = .8,slab_colour = "black",slab_fill="#3498db", slab_size = .5) + 
  #scale_discrete_manual(values = c("#3498db")) + 

  theme_classic() + 
  theme(axis.text.x = element_text(color = "black", size = 22, face = "bold"), axis.text.y = element_text(color = "black", size = 22, face = "bold"), axis.title.x = element_text(color = "black", size = 26, face = "bold"), axis.title.y = element_text(color = "black", size = 26, face = "bold"))+
  xlab(expression(paste("Condition Effect Size ", delta, " (rep6-dev)"))) +
  ylab("Density") 
  
  # Save out
  figfname <- paste(FigOutDir, "TFcon_conddiff_RIFG_Bayes_Eff.tiff", sep="/")

  ggsave(figfname, dpi=150)

# Again
mydata <- as_tibble(samples[,3]) %>%
  tidyr::pivot_longer(everything())
  
g1<-ggplot(data.frame(x = c(-2, 2)), aes(x)) +
  stat_function(fun = dcauchy, n = 10000, args = list(location = 0, scale = 0.707), size = 1, linetype = "dotted") + 
    ggdist::stat_halfeye(data=mydata, aes(x=value), .width = c(0.95),alpha = .8,slab_colour = "black",slab_fill="#3498db",slab_linewidth = 1, slab_size = .5,  linewidth=5) +
  labs(x = expression(bold(paste("Condition Effect Size ", delta, " (rep6-dev)")))) +
  ylab("Density") +
  theme_classic() + 
  
  theme(axis.text.x = element_text(color = "black", size = 20, face = "bold", family="DejaVu Sans"), axis.text.y = element_text(color = "black", size = 20, face = "bold", family="DejaVu Sans"), axis.title.x = element_text(color = "black", size = 22, face = "bold", family="DejaVu Sans"), axis.title.y = element_text(color = "black", size = 22, face = "bold", family="DejaVu Sans"))

figfname <- paste(FigOutDir, "TFcon_conddiff_RIFG_Bayes_Eff.tiff", sep="/")

ggsave(figfname, g1, w=8, h=6, dpi=100)


# Repeat to get legend
g1 + ggdist::stat_halfeye(data=mydata, aes(x=value), .width = c(0.95),alpha = .8,slab_colour = "black",slab_fill="#3498db",slab_linewidth = 1, slab_size = .5,  linewidth=5, show.legend = TRUE)

figfname <- paste(FigOutDir, "TFcon_conddiff_RIFG_Bayes_Eff_wleg.tiff", sep="/")
ggsave(figfname, w=8, h=6, dpi=150)

# Group differences (CON vs.bvFTD/PSP) ---------------------------------------------

LFPdata_wide <- read.csv("adv_ssst_newmaxf_fixICA_wfids_250_noTGBcon_FullDiffPowTable.csv", header = TRUE)
  
df_nodbl <- subset(LFPdata_wide, double == 0)

bvFTD <- subset(df_nodbl, Diag == 1)
PSP <- subset(df_nodbl, Diag == 2)
CN <- subset(df_nodbl, Diag == 3)

## Sample from the corresponding posterior distribution
samples = ttestBF(x = CN$RIFG, y = bvFTD$RIFG, paired=FALSE,
                    posterior = TRUE, iterations = 10000)
samples_2 = ttestBF(x = CN$RIFG, y = PSP$RIFG, paired=FALSE,
                  posterior = TRUE, iterations = 10000)

mydata <- as_tibble(cbind2(samples[,4], samples_2[,4])) %>%
  tidyr::pivot_longer(everything())

mydata <- data.frame(cbind2(samples[,4], samples_2[,4])) 

colnames(mydata) <- c("bvFTD", "PSP")

g2<-ggplot(data.frame(x = c(-2, 2)), aes(x)) +
  stat_function(fun = dcauchy, n = 10000, args = list(location = 0, scale = 0.707), size = 1, linetype = "dotted") + 
  ggdist::stat_halfeye(data=mydata, aes(x=bvFTD), .width = c(0.95), slab_alpha = .4,slab_colour = "black",slab_fill="#E69F00",interval_color="#E69F00",slab_linewidth = 1, slab_size = .5,  linewidth=20) +
  ggdist::stat_halfeye(data=mydata, aes(x=PSP), .width = c(0.95), slab_alpha = .4,slab_colour = "black",slab_fill="#00BA38", interval_color="#00BA38", slab_linewidth = 1, slab_size = .5,  linewidth=20) +
  labs(x = expression(bold(paste("Group Effect Size ", delta, " (", Delta, "Beta1)")))) +
  ylab("Density") +
  theme_classic() + 
  theme(axis.text.x = element_text(color = "black", size = 25, face = "bold", family = "DejaVu Sans"), axis.text.y = element_text(color = "black", size = 25, face = "bold", family = "DejaVu Sans"), axis.title.x = element_text(color = "black", size = 28, face = "bold", family = "DejaVu Sans"), axis.title.y = element_text(color = "black", size = 28, face = "bold", family = "DejaVu Sans"))

figfname <- paste(FigOutDir, "TFcon_grpeff_RIFG_Bayes_both.tiff", sep="/")

ggsave(figfname, g2, w=9, h=7, dpi=100)

# Repeat with legends

# bvFTD
g2 + ggdist::stat_halfeye(data=mydata, aes(x=bvFTD), .width = c(0.95), slab_alpha = .4,slab_colour = "black",slab_fill="#E69F00",interval_color="#E69F00",slab_linewidth = 1, slab_size = .5,  linewidth=20, show.legend = TRUE)

figfname <- paste(FigOutDir, "TFcon_grpeff_RIFG_Bayes_both_bvFTDleg.tiff", sep="/")

ggsave(figfname, w=8, h=6, dpi=150)

# PSP
g2 + ggdist::stat_halfeye(data=mydata, aes(x=PSP), .width = c(0.95), slab_alpha = .4,slab_colour = "black",slab_fill="#00BA38",interval_color="#00BA38",slab_linewidth = 1, slab_size = .5,  linewidth=20, show.legend = TRUE)

figfname <- paste(FigOutDir, "TFcon_grpeff_RIFG_Bayes_both_PSPleg.tiff", sep="/")

ggsave(figfname, w=8, h=6, dpi=150)



