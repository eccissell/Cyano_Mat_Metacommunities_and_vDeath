#Combined Script for both manuscripts included in this repo
#CODE AUTHOR: ETHAN C. CISSELL
#DATE LAST CODE UPDATED: AUGUST 22, 2025
#####################################################
#SECTION 1
#Code for: Cissell EC, McCoy SJ. In Review @ The American Naturalist. Heterogenous predation, community asynchrony, and metacommunity stability in cyanobacterial mats (Preprint Available Here: https://www.biorxiv.org/content/10.1101/2022.10.07.511315v2)
#####################################################
#Load libraries
require(plyr)
library(ggbreak)
library(ggpubr)
library(readxl)
library(gdata)
library(gtools)
library(gdata)
require(ggplot2)
library(reshape2)
library(colorspace)
library(gtable)
library(labeling)
library(Rcpp)
library(digest)
library(magrittr)
library(stringi)
library(munsell)
library(stringr)
library(scales)
library(grid)
require(ggthemes)
library(lme4)
library(lmerTest)
library(gridExtra)
library(nlme)
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(Rmisc)
library(ggpubr)
library(fishualize)
library(plotly)
library(mgcv)
library(mgcViz)
library(rgl)
library(plot3D)
library(plotrix)
library(dplyr)
library(matrixStats)
library(car)

########
##Metacommunity benthic cover
########
#Read in data
benthic=read.csv('/Users/cissell/Desktop/metacommunity_BCM_benthic_df.csv')
#View
View(benthic)
#Summarize at level of photoquad
benthos = benthic %>%
  group_by(Name, Label, Date, Depth) %>%
  dplyr::summarise(n()/50)
#View(benthos)
colnames(benthos)[5] = "prop"
#Convert Long to Short Format
benthos_wide=spread(benthos,Label,prop)
View(benthos_wide)
summary(benthos_wide)
#Convert NA to 0
benthos_wide[is.na(benthos_wide)]=0
View(benthos_wide)
#Look at all relationships
ggpairs(benthos_wide,column=2:11,progress=FALSE)
##Focus on spread of depth over time
#Calc CI
DEPTH<- summarySE(benthos_wide, groupvars = c("Date"), measurevar = c("Depth"))
View(DEPTH)
#Now Pass that to plot of means
DEPTH %>% ggplot()+
  geom_bar(aes(x = Date, y = Depth), stat= "summary", fun="mean", fill="#393E42") +
  theme_classic() +
  ylab("Depth") +
  scale_fill_brewer(palette="Paired", name = "Substrate") +
  theme(axis.ticks.x=element_blank()) +
  theme(axis.title = element_text(size=18),
        axis.text.x = element_text(size=16, color = "black"),
        axis.text.y = element_text(size=16, color = "black")) +
  theme(legend.text=element_text(size=16), legend.title = element_text(size = 18)) +
  geom_errorbar(aes(x= Date, y= Depth, ymin=Depth-ci, ymax=Depth+ci),
                width=0.09,
                color="#fe5d26",
                size=1)
#Boxplot for quick look at spread
boxplot(Depth~Date,las=1, data=benthos_wide)
dev.print(png,fil="/Users/cissell/Desktop/depth_spread.png",type="quartz",antialias="default",width=5.5,
          height=4.5,units="in",res=1300)
#Model to look at BCM across date and depth
lmbcmmet=lm(BCM~Date*Depth,data=benthos_wide)
summary(lmbcmmet)
#Check Assumptions
library("DHARMa")
check_bcmmet_model1 <- simulateResiduals(fittedModel = lmbcmmet, n = 500)
plot(check_bcmmet_model1)
#Good
#Significance
Anova(lmbcmmet,type="II",test.statistic = "F")
#Run model without interaction because interactions isnt significatn, then compare with LRT
lmbcmmet2=lm(BCM~Date+Depth,data=benthos_wide)
summary(lmbcmmet2)
#Check Assumptions
library("DHARMa")
check_bcmmet_model2 <- simulateResiduals(fittedModel = lmbcmmet2, n = 500)
plot(check_bcmmet_model2)
#Good
#Now test if adding complexity is good using LRT
anova(lmbcmmet2,lmbcmmet,test="LRT") #p=0.269
#Use reduced model
#Significance of reduced model (no interaction)
Anova(lmbcmmet2,type="II",test.statistic = "F")
#Now use Tukey HSD
#Now use pairwise Wilcox with corrections to look at specific differences
pairwise.wilcox.test(benthos_wide$BCM, benthos_wide$Date, p.adjust.method="bonferroni")
#Only 6/11 vs. 6/13, 6/20 and 6/26 differ
#Summarize and Organize Benthic Data by Date
BCM<- summarySE(benthos_wide, groupvars = c("Date"), measurevar = c("BCM"))
colnames(BCM)[3] <- "Mean"
BCM$Group = "BCM"
summary(BCM)

EAM<- summarySE(benthos_wide, groupvars = c("Date"), measurevar = c("EAM"))
colnames(EAM)[3] <- "Mean"
EAM$Group="EAM"

Coral<-summarySE(benthos_wide, groupvars = c("Date"), measurevar = c("HC"))
colnames(Coral)[3] <- "Mean"
Coral$Group="Coral"

Macro<-summarySE(benthos_wide, groupvars = c("Date"), measurevar = c("MAL"))
colnames(Macro)[3] <- "Mean"
Macro$Group="Macroalgae"

Sediment<-summarySE(benthos_wide, groupvars = c("Date"), measurevar = c("Sand"))
colnames(Sediment)[3] <- "Mean"
Sediment$Group="Sediment"

Gorgonian<-summarySE(benthos_wide, groupvars = c("Date"), measurevar = c("SC"))
colnames(Gorgonian)[3] <- "Mean"
Gorgonian$Group="Gorgonian"

Sponge<-summarySE(benthos_wide, groupvars = c("Date"), measurevar = c("SP"))
colnames(Sponge)[3] <- "Mean"
Sponge$Group="Sponge"

Unk<-summarySE(benthos_wide, groupvars = c("Date"), measurevar = c("Unk"))
colnames(Unk)[3] <- "Mean"
Unk$Group="Unknown"

#combine
bs <- rbind(BCM,EAM,Coral,Macro,Sediment,Gorgonian,Sponge,Unk)
View(bs)

#Plot cover of BCM against date
bcm_date=bs %>%
  filter(bs$Group == "BCM") %>%
  ggplot()+
  geom_bar(aes(x = Date, y = Mean), stat= "identity", fill="#393E42",size=1.5,color="black") +
  theme_classic() +
  ylab("Proportional Cover of BCM") +
  scale_fill_brewer(palette="Paired", name = "Substrate") +
  theme(axis.ticks.x=element_blank()) +
  theme(axis.title = element_text(size=16),
        axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        axis.line.y = element_line(size=1.1),
        axis.line.x = element_line(size=1.1)) +
  theme(legend.text=element_text(size=16), legend.title = element_text(size = 18)) +
  geom_errorbar(aes(x= Date, y= Mean, ymin=Mean-ci, ymax=Mean+ci),
                width=0.09,
                color="black",
                size=1) +
  geom_hline(yintercept=0.1956,
             linetype="dashed",
             color="black",
             size=1.07)+
  xlab("")
bcm_date

#Get mean and SD of BCM cover across dates
bs_bcmsub = bs %>%
  subset(Group=="BCM")
mean(bs_bcmsub$Mean) #0.196
sd(bs_bcmsub$Mean) #0.0494

#Summarize and Organize Benthic Data by Depth 
#for Supplementary plot of BCM cover vs depth
BCM2<- summarySE(benthos_wide, groupvars = c("Depth"), measurevar = c("BCM"))
colnames(BCM2)[3] <- "Mean"
BCM2$Group = "BCM"

EAM2<- summarySE(benthos_wide, groupvars = c("Depth"), measurevar = c("EAM"))
colnames(EAM2)[3] <- "Mean"
EAM2$Group="EAM"

Coral2<-summarySE(benthos_wide, groupvars = c("Depth"), measurevar = c("HC"))
colnames(Coral2)[3] <- "Mean"
Coral2$Group="Coral"

Macro2<-summarySE(benthos_wide, groupvars = c("Depth"), measurevar = c("MAL"))
colnames(Macro2)[3] <- "Mean"
Macro2$Group="Macroalgae"

Sediment2<-summarySE(benthos_wide, groupvars = c("Depth"), measurevar = c("Sand"))
colnames(Sediment2)[3] <- "Mean"
Sediment2$Group="Sediment"
View(Sediment2)

Gorgonian2<-summarySE(benthos_wide, groupvars = c("Depth"), measurevar = c("SC"))
colnames(Gorgonian2)[3] <- "Mean"
Gorgonian2$Group="Gorgonian"

Sponge2<-summarySE(benthos_wide, groupvars = c("Depth"), measurevar = c("SP"))
colnames(Sponge2)[3] <- "Mean"
Sponge2$Group="Sponge"

Unk2<-summarySE(benthos_wide, groupvars = c("Depth"), measurevar = c("Unk"))
colnames(Unk2)[3] <- "Mean"
Unk2$Group="Unknown"

#Combine
bs2 <- rbind(BCM2,EAM2,Coral2,Macro2,Sediment2,Gorgonian2,Sponge2,Unk2)
View(bs2)

#Change to factor with levels from Character
bs2$Group=as.factor(bs2$Group)

#Plot BCM cover vs. depth
bs2 %>%
  filter(bs2$Group == "BCM") %>%
  ggplot()+
  geom_smooth(aes(x = Depth, y = Mean),method=lm,se=T,color="black")+
  geom_point(aes(x = Depth, y = Mean), stat= "identity") +
  theme_classic() +
  ylab("Proportional Cover of BCM") +
  xlab("Depth (m)") +
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=13, color = "black"),
        axis.text.y = element_text(size=13, color = "black"),
        axis.line.y = element_line(size=1.1),
        axis.line.x = element_line(size=1.1)) +
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 18))

benthos_wide%>%
  ggplot()+
  geom_smooth(aes(x = Depth, y = BCM),method=lm,se=T,color="black")+
  geom_point(aes(x = Depth, y = BCM), stat= "identity") +
  theme_classic() +
  ylab("Proportional Cover of BCM") +
  xlab("Depth (m)") +
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=13, color = "black"),
        axis.text.y = element_text(size=13, color = "black"),
        axis.line.y = element_line(size=1.1),
        axis.line.x = element_line(size=1.1)) +
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 18))
dev.print(png,fil="/Users/cissell/Desktop/depth_spread.png",type="quartz",antialias="default",width=5.5,
          height=4.5,units="in",res=1300)
lm (BCM~Depth, data = benthos_wide)

########
##Community persistence and coring experiment
########
#Looking at overall mat deaths
df=read.csv('/Users/cissell/Desktop/coring_experiment_df.csv')
View(df)

#Sort by death date
df=df[order(df$Last_Date),]

#Lets plot the mean + SE
p=ggplot(df,aes(x=Treatment,y=Binom)) +
  geom_point(alpha=0.8) +
  geom_jitter(width=0.2,height=0.07) +
  theme_classic()
p= p + stat_summary(fun.data="mean_se",fun.args = list(mult=1),
                    geom="pointrange",color="black",size=1) +
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  scale_x_discrete(labels=c("Wounded","Unmanipulated"))+
  labs(y="Probability of Survival")
p
ggsave("/Users/cissell/Desktop/overall_deaths.svg",width=3.25,height=4,units="in",p)

#Boxplot of deaths by date
p2=ggplot(data=df, aes(x=Last_Date, y=Binom)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))
p2

#Cumulative sum plot of proportion dead by date
#Need to change 0 and 1 assignment for this plot that counts cumulative dead
df_zerooneswap=df %>%
  subset(Treatment=="W")%>%
  mutate(binom_swap = ifelse(Binom==1,0,1))
#Order
#Sort by death date
df_zerooneswap=df_zerooneswap[order(df_zerooneswap$Last_Date),]
#Cumulative sum plot of proportion dead by date
death_cumulsum=ggplot(data=df_zerooneswap, aes(x=Last_Date)) +
  geom_col(aes(y=binom_swap))+
  geom_line(aes(x=Last_Date,y=(cumsum(binom_swap)/26*100),group=1), inherit.aes=FALSE, size=1.1)+
  theme(panel.background = element_blank()) +
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  scale_y_continuous(name="Cumulative % dead")+
  scale_x_discrete(name="Date")
#  geom_vline(xintercept=3,
#             linetype="dashed",
#             color="black",
#             size=1.07,
#             alpha=0.4) +
#  geom_vline(xintercept=4,
#             linetype="dashed",
#             color="black",
#             size=1.07,
#             alpha=0.4)+
#  geom_vline(xintercept=8,
#             linetype="dashed",
#             color="black",
#             size=1.07,
#             alpha=0.4)+
#  geom_vline(xintercept=10,
#             linetype="dashed",
#             color="black",
#             size=1.07,
#             alpha=0.4)+
#  geom_vline(xintercept=14,
#             linetype="dashed",
#             color="black",
#             size=1.07,
#             alpha=0.4)
death_cumulsum

#14 died in unmanipulated treatment across 10 different death dates
#53.8% mortality (14/26)

#Combine this plot with metacommunity persistence plot from above for Fig 1
bcm_ar=ggarrange(bcm_date, NULL, death_cumulsum,
                 ncol=1, nrow=3, align="h",
                 widths=c(1,-0.1,1),
                 heights=c(1,-0.18,1.1),
                 labels=c("a","","b"))
ggsave("/Users/cissell/Desktop/metcom_vs_com.svg",width=6.5,height=8.5,units="in",bcm_ar)

#Use Kruskal-Wallis test to initial test for differences between treatments in start date
kruskal.test(First_Date~Treatment,data=df) #p=0.8821 ; no significant difference
#Use GLM for differences between treatments
trtglm=glm(Binom~Treatment,data=df,family=binomial)
summary(trtglm)
#Check Assumptions
library("DHARMa")
check_trt_model1 <- simulateResiduals(fittedModel = trtglm, n = 500)
plot(check_trt_model1)
#Looks good

#Look for effect of start date
trtglm2=glm(Binom~Treatment*First_Date,data=df,family=binomial)
summary(trtglm2)
#Check Assumptions
library("DHARMa")
check_trt_model2 <- simulateResiduals(fittedModel = trtglm2, n = 500)
plot(check_trt_model2)
#Looks good, and no significatn effects
#Finally test if adding complexity is good using LRT
anova(trtglm,trtglm2,test="LRT") #p=0.02 ; Marginally better use full model

#Now test for significance
Anova(trtglm2,type="II",test.statistic = "F")

#Fit decay curve with glm
exp1=glm(Binom ~ Duration : Treatment-1, data = df, family = binomial(link=log))
#Summarise
summary(exp1)
#Test for significant effect of treatment
anova(exp1, test= "Chisq")
#Extract CI around parameter estimates
confint(exp1)
#Check Assumptions
library("DHARMa")
check_exp1 <- simulateResiduals(fittedModel = exp1, n = 500)
plot(check_exp1)
testDispersion(check_exp1)

#Fails Overdispersion 
#Fit decay curve with quasibinomial to overcome overdispersion
exp2=glm(Binom ~ Duration : Treatment-1, data = df, family = quasibinomial(link=log))
#Summarise
summary(exp2)
#Test for significant effect of treatment
anova(exp2, test= "Chisq")
#Look at CI around parameter estimates
confint(exp2)

#Create predict dataframes and extract predictions for each group
preddfW=data.frame(Treatment=factor("W",levels=c("H","W")), Duration = seq(min(df$Duration),max(df$Duration), length=100))
preddfH=data.frame(Treatment=factor("H",levels=c("H","W")), Duration = seq(min(df$Duration),max(df$Duration), length=100))

#Append the predicts to the dfs
preddfW$Fit=predict(exp2, preddfW, se.fit = TRUE, type = "response")$fit
preddfH$Fit=predict(exp2, preddfH, se.fit = TRUE, type = "response")$fit


#Plot the curves
decay_curve_plot=ggplot(data=df,aes(x=Duration, y=Binom, group=Treatment)) +
  geom_point(aes(color=Treatment),
             size=2,
             position=position_jitter(w=0,h=0.1))+
  scale_color_manual(values=c("#cd9fb2","#455765"),labels=c("W","U"))+
  geom_line(data=preddfW,aes(x=Duration,y=Fit),
            size=1.75,
            color="#455765")+
  geom_line(data=preddfH,aes(x=Duration,y=Fit),
            size=1.75,
            color="#cd9fb2")+
  theme(panel.background = element_blank())+
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line.y = element_line(size=1.1),
        axis.line.x = element_line(size=1.1)) +
  ylab("Probability of survival") +
  xlab("Time (days)")

decay_curve_plot
ggsave("/Users/cissell/Desktop/decay_curve.svg",width=6.5,height=4.5,units="in",decay_curve_plot)

########
##Fish predation pressure
########
#Read in data
matfollow=read.csv('/Users/cissell/Desktop/fish_bites_df.csv')
#Check the data
View(matfollow)
#Summarize to obtain freq by species
library(tidyr)
gathered=gather(matfollow,species,freq,bicolor,yellow_goat,tp_stripe,ip_princess,ip_stripe,queen,surgeon,tang,spotted_goat,ip_redband,rock,french,tp_princess,tp_redband)
View(gathered)
gathered$species=as.factor(gathered$species)
str(gathered)
#Fit model
m1=glmmTMB(data=gathered, freq ~ site+species + (1|site/mat) + offset(log(duration)),family=nbinom2,REML=TRUE)
summary(m1)
Anova(m1,type="II")
#Check model assumptions
simulationOutput <- simulateResiduals(m1)
testDispersion(simulationOutput)
plot(simulationOutput)

ggsave("/Users/cissell/Desktop/site_bites.svg",width=4,height=4.5,units="in",site_bites)

#Obtain SD of total bites among mats to better understand among mat variability suggested in RE
#All mats
sd(matfollow$total.bites) #111.2 bites SD
mean(matfollow$total.bites) #72.8 bites Mean

#Plot of bite count by species for Fig 3
g=gathered %>%
  mutate(species=species%>%
           fct_relevel("bicolor",
                       "tp_stripe",
                       "ip_stripe",
                       "tp_princess",
                       "ip_princess",
                       "tp_redband",
                       "ip_redband",
                       "queen",
                       "surgeon",
                       "tang",
                       "rock",
                       "french",
                       "yellow_goat",
                       "spotted_goat"))
fish_species_bites=ggplot(data=g,aes(x=freq, y=as.factor(species)))+
  labs(y="Species") +
  labs(x="Bite Count") +
  geom_boxplot(notch=FALSE, fill='#A4A4A4',outlier.shape=NA)+
  geom_jitter(data = g, height=0.1, width=0, aes(alpha=8/10),size=3)+
  scale_y_discrete(labels=c("St. partitus", "TP Sc. iseri","IP Sc. iseri","TP Sc. taeniopterus","IP Sc. taeniopterus","TP Sp. aurofrenatum","IP Sp. aurofrenatum","Sc. vetula","A. bahianus","A. coeruleus","H. tricolor","Po. paru","M. martinicus","Ps. maculatus")) +
  theme(panel.grid.major = element_blank (), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.line = element_line(size=1, colour = "black"),
        axis.line.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()) +
  guides (size=FALSE) + guides(alpha=FALSE) +
  theme(axis.text.y=element_text(face="italic")) +
  scale_x_break(c(150,220))
fish_species_bites
ggsave("/Users/cissell/Desktop/species_bites.svg",width=6.5,height=7.5,units="in",fish_species_bites)

#Summarize to obtain rate by species
gathered2=gather(matfollow,rates,rate,r_bi,r_yel,r_tstripe,r_iprince,r_qu,r_surgeon,r_tang,r_spot,r_ired,r_rock,r_french,r_tprince,r_tred)
View(gathered2)
mean(gathered2$rate) #0.2346403
sd(gathered2$rate) #0.922

#Plot of bites across sites
#With damselfish
site_bites_with=ggplot(data=matfollow, aes(x=site, y=total.bites)) +
  geom_boxplot(aes(x=site,y=total.bites, group=site),
               fill="#A4A4A4",
               size=0.75,
               color="black") +
  theme_classic() +
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt")) +
  scale_y_break(c(275,450)) +
  ylab("Total Bites (#)") +
  xlab("Site")
site_bites_with
ggsave("/Users/cissell/Desktop/site_bites.tiff",width=6.5,height=7.5,units="in",site_bites_with)

#Invert axes for consistent plotting in revised figure 4
site_bites_with=ggplot(data=matfollow, aes(x=total.bites, y=site)) +
  geom_boxplot(aes(x=total.bites,y=site, group=site),
               fill="#A4A4A4",
               size=0.75,
               color="black") +
  theme_classic() +
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt")) +
  ylab("Study Site") +
  xlab("Total Bites (#)")
site_bites_with
ggsave("/Users/cissell/Desktop/site_bites.svg",width=6.1,height=3.5,units="in",site_bites_with)

######END OF CODE FOR SECTION 1#########
################################################################
################################################################
################################################################
################################################################
#SECTION 2:
#Code for Cissell EC, McCoy SJ. Viral association with cyanobacterial mat community mortality. Ecology. 2023 Sep;104(9):e4131. doi: 10.1002/ecy.4131. (Link: https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.4131).
################################################################
########
##Viral predation pressure
########
##Load in required packages
library(Nonpareil)

#GENERATE MODEL DAL
npDAL=Nonpareil.curve("/Users/cissell/Desktop/DAL.npo", col = "#f15538")
#OBTAIN PARAMETER ESTIMATES
print(npDAL)
#Check Model Summary
npDAL

#GENERATE MODEL DAB
npDAB=Nonpareil.curve ("/Users/cissell/Desktop/DAB.npo", col = "#fbaf37")
#OBTAIN PARAMETER ESTIMATES
print(npDAB)
#Check Model Summary
npDAB

#GENERATE MODEL DAD
npDAD=Nonpareil.curve ("/Users/cissell/Desktop/DAD.npo", col = "#4fa3a4")
#OBTAIN PARAMETER ESTIMATES
print(npDAD)
#Check Model Summary
npDAD

#Make list for legend
np_list=list(npDAL, npDAB, npDAD)

##GENERATE COMBINED PLOT OF MODEL OUTPUTS
plot(npDAL, col="#f15538",plot.observed=FALSE,plot.model=TRUE,plot.diversity=FALSE,las=1,model.lwd=5,main='',xlim=c(1e7,1e13))
plot(npDAB, col="#fbaf37",add=TRUE, plot.observed=FALSE,plot.model=TRUE,plot.diversity=FALSE,las=1,model.lwd=5, main='')
plot(npDAD, col="#4fa3a4",add=TRUE, plot.observed=FALSE,plot.model=TRUE,plot.diversity=FALSE,las=1,model.lwd=5, main='')
#plot(npDBL, col="#9cbad6",add=TRUE, plot.observed=FALSE,plot.model=TRUE,plot.diversity=FALSE,las=1,model.lwd=5, main='')
Nonpareil.legend(np_list,x=1e11,cex=0.85,bty="n", y=0.6)

#Read in tsvs
setwd("/Users/cissell/Desktop/VIRSORTER2_Death")
deathviral_tsv=dir(patter=".tsv")
n = length(deathviral_tsv)
deathviral_tsv_list=vector("list",n)
for(i in 1:n) deathviral_tsv_list[[i]] = read.table(deathviral_tsv[i],header=T)

#Merge by contig id
DA=merge(deathviral_tsv_list[[1]],deathviral_tsv_list[[2]],by="contig_id",all=T)

#Now apply manual curation
#Keep1:viral_gene>0
DAkeep1=DA%>%
  subset(viral_genes>0)
DAkeep1_remain=DA%>%
  subset(!viral_genes>0)
#Keep2_1:host_gene=0
DAkeep2_1=DAkeep1_remain%>%
  subset(host_genes==0)
DAkeep2_1_remain=DAkeep1_remain%>%
  subset(!host_genes==0)
#Keep2_2:score>=0.95
DAkeep2_2=DAkeep2_1_remain%>%
  subset(max_score>=0.95)
DAkeep2_2_remain=DAkeep2_1_remain%>%
  subset(!max_score>=0.95)
#Keep2_3:hallmark>2
DAkeep2_3=DAkeep2_2_remain%>%
  subset(hallmark>2)
DAkeep2_3_remain=DAkeep2_2_remain%>%
  subset(!hallmark>2)
#Manual:Host<=1 & length >=10000
DAkeepman=DAkeep2_3_remain%>%
  subset(host_genes<=1 & length>=10000)
#Man merge
DAmanmerge=rbind(DAkeepman,DAkeep2_3)
#Now screen those contigs in manmerge manually
#remove DAL_41||full
##
DAmanmerge=DAmanmerge%>%
  filter(!row_number() %in% c(2))

#Now merge all the kept phages (yay!)
deathviral_curated=bind_rows(DAkeep1,DAkeep2_1,DAkeep2_2,DAmanmerge)

#Ensure deduped
deathviral_dedup=deathviral_curated[!duplicated(deathviral_curated),] #Fine

#Export combined as a tsv
write.table(deathviral_curated[,1], file="/Users/cissell/Desktop/VIRSORTER2_Death/deathcurated_ids.tsv", quote=F, sep='\t', col.names=NA)

##Now we have depreplicated at 95% within mat
deathdedup_ids=read.table("/Users/cissell/Desktop/VIRSORTER2_Death/death_dedup_headers.tsv",header=T)

#Remove IDS after || in curated dataframe for merging
#Just do it in Excel becuase it is easier.
write.table(deathviral_curated, file="/Users/cissell/Desktop/VIRSORTER2_Death/curat.tsv", quote=F, sep='\t', col.names=NA)
#Now read back in
deathviral_curated_simp=read.table("/Users/cissell/Desktop/VIRSORTER2_Death/curat.tsv",header=T)

#Now merge to subset by clustered at 95%
deathviral_clust_curate=merge(deathviral_curated_simp,deathdedup_ids,by="contig_id")

#Let's clean up lengths for pro vs non pro regions
#Start by splitting into 2 dataframes, pro vs non-pro
deathviral_clust_curate_pro=deathviral_clust_curate%>%
  subset(provirus=="Yes")
deathviral_clust_curate_non=deathviral_clust_curate%>%
  subset(provirus=="No")

#clean up columns from non
deathviral_clust_curate_non_clean=deathviral_clust_curate_non%>%
  select(-c(proviral_length,
            host_length,
            region_types,
            region_lengths,
            region_coords_bp,
            region_coords_genes,
            region_viral_genes,
            region_host_genes,
            length))
#clean up columns in pro and rename length column
deathviral_clust_curate_pro_clean=deathviral_clust_curate_pro%>%
  select(-c(host_length,
            region_types,
            region_lengths,
            region_coords_bp,
            region_coords_genes,
            region_viral_genes,
            region_host_genes,
            length,
            contig_length)) %>%
  rename(contig_length=proviral_length)
#recombine
deathviral_clust_curate_clean=rbind(deathviral_clust_curate_non_clean,deathviral_clust_curate_pro_clean)

#pull mean length
mean(deathviral_clust_curate_clean$contig_length) #6647.2
median(deathviral_clust_curate_clean$contig_length) #3638


#####Get viral abundances from DNA mapping and make violins
#We will do this separately when looking at cyanophage hosts so that total host + phage adds to 1e6 and here total phage adds to 1e6 so not directly comparable
#Read in abundances
death_viral_abunds=read.table("/Users/cissell/Desktop/VIRSORTER2_Death/DA_phage_counts.tsv",header=T)
#Make new df
death_po=death_viral_abunds
str(death_po)
#Put length in KB
death_po$Plength=death_po$Plength/1000
#Create RPK
death_po$PDAB=death_po$PDAB/death_po$Plength
death_po$PDAD=death_po$PDAD/death_po$Plength
death_po$PDAL=death_po$PDAL/death_po$Plength
#create scaling factors
sumDABp=(sum(death_po$PDAB)/1000000)
sumDADp=(sum(death_po$PDAD)/1000000)
sumDALp=(sum(death_po$PDAL)/1000000)

#Create TPM
death_po$PDAB=death_po$PDAB/sumDABp
death_po$PDAD=death_po$PDAD/sumDADp
death_po$PDAL=death_po$PDAL/sumDALp

#Bring it back to a sensical name now
Deathcount_tpm_phageonly=death_po
#Check for sense
sum(Deathcount_tpm_phageonly$PDAB)
sum(Deathcount_tpm_phageonly$PDAD)
sum(Deathcount_tpm_phageonly$PDAL)

#Pivot longer for plotting and do log
Deathcount_tpm_phageonly_long=Deathcount_tpm_phageonly%>%
  select(-Plength)%>%
  pivot_longer(
    cols=c(PDAB,PDAD,PDAL),
    names_to="Sample",
    values_to="TPM")%>%
  mutate(logTPM=log10(TPM+1))
Deathcount_tpm_phageonly_long$Sample=as.factor(Deathcount_tpm_phageonly_long$Sample)
#Redo the order of Sample for good-good plotting
Deathcount_tpm_phageonly_long=Deathcount_tpm_phageonly_long%>%
  mutate(Sample=fct_relevel(Sample,
                            "PDAL",
                            "PDAB",
                            "PDAD"))
#Find median numbers
aggregate(logTPM~Sample,data=Deathcount_tpm_phageonly_long,summary)
#Median DAL = 2.27
#Median DAB = 1.64
#Median DAD = 2.41

#Now violin
#with zeros
death_virusviolin=ggplot(Deathcount_tpm_phageonly_long, aes(y=logTPM,x=Sample))+
  geom_violin(aes(fill=Sample),
              outlier.shape = NA,
              size=1,
              alpha=0.9,
              draw_quantiles=c(0.5))+
  scale_fill_manual(values=c("#f15538","#fbaf37","#4fa3a4"),name="Time")+
  geom_point(aes(fill=Sample),position=position_jitterdodge(0.2),colour="black",pch=21,
             alpha=0.7,
             size=1.5)+
  geom_violin(aes(fill=Sample),
              outlier.shape = NA,
              size=1,
              alpha=0,
              draw_quantiles=c(0.25,0.5,0.75))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Viral CPM ("~log[10]~')'),
       x="Sample")+
  theme(legend.position="none")
death_virusviolin
ggsave("/Users/cissell/Desktop/viral_props_death.svg",death_virusviolin,width=6,height=3)

#Lets also see if we remove those that ever are 0
Deathcount_tpm_phageonly_longnozero=Deathcount_tpm_phageonly_long%>%
  filter(logTPM>0)
#Summary
aggregate(logTPM~Sample,data=Deathcount_tpm_phageonly_longnozero,summary)

#Violin without zeros
death_virusviolinnozero=ggplot(Deathcount_tpm_phageonly_longnozero, aes(y=logTPM,x=Sample))+
  geom_violin(aes(fill=Sample),
              outlier.shape = NA,
              size=1,
              alpha=0.9,
              draw_quantiles=c(0.5))+
  scale_fill_manual(values=c("#f15538","#fbaf37","#4fa3a4"),name="Time")+
  geom_point(aes(fill=Sample),position=position_jitterdodge(0.2),colour="black",pch=21,
             alpha=0.7,
             size=1.5)+
  geom_violin(aes(fill=Sample),
              outlier.shape = NA,
              size=1,
              alpha=0,
              draw_quantiles=c(0.25,0.5,0.75))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Viral CPM ("~log[10]~')'),
       x="Sample")+
  theme(legend.position="none")
death_virusviolinnozero
ggsave("/Users/cissell/Desktop/viral_props_death_nozero.png",death_virusviolinnozero,width=6.5,height=4)

#Create 3 separate rank abundance curve plots to look at patterns in eveness
#Sort by abundance from each sample in separate df
death_virus_DAL_rank_sort = Deathcount_tpm_phageonly %>%
  select(contig_id,PDAL) %>%
  arrange(desc(PDAL))
death_virus_DAL_rank_sort$contig_id=factor(death_virus_DAL_rank_sort$contig_id, levels = death_virus_DAL_rank_sort$contig_id)

death_virus_DAB_rank_sort = Deathcount_tpm_phageonly %>%
  select(contig_id,PDAB) %>%
  arrange(desc(PDAB))
death_virus_DAB_rank_sort$contig_id=factor(death_virus_DAB_rank_sort$contig_id, levels = death_virus_DAB_rank_sort$contig_id)

death_virus_DAD_rank_sort = Deathcount_tpm_phageonly %>%
  select(contig_id,PDAD) %>%
  arrange(desc(PDAD))
death_virus_DAD_rank_sort$contig_id=factor(death_virus_DAD_rank_sort$contig_id, levels = death_virus_DAD_rank_sort$contig_id)

#Now plot the abundance 'curves'
rank_DAL=ggplot(death_virus_DAL_rank_sort, aes(y=log10(1+PDAL),x=contig_id))+
  geom_bar(stat="identity", 
           color="#f15538",
           fill = "#f15538")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank())+
  labs(y=expression("Viral CPM ("~log[10]~')'),
       x="Viral rank")+
  theme(legend.position="none")
rank_DAL
ggsave("/Users/cissell/Desktop/RAC_DAL.svg",rank_DAL,width=6.5,height=3)

rank_DAB=ggplot(death_virus_DAB_rank_sort, aes(y=log10(1+PDAB),x=contig_id))+
  geom_bar(stat="identity", 
           color="#fbaf37",
           fill = "#fbaf37")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank())+
  labs(y=expression("Viral CPM ("~log[10]~')'),
       x="Viral rank")+
  theme(legend.position="none")
rank_DAB
ggsave("/Users/cissell/Desktop/RAC_DAB.svg",rank_DAB,width=6.5,height=3)

rank_DAD=ggplot(death_virus_DAD_rank_sort, aes(y=log10(1+PDAD),x=contig_id))+
  geom_bar(stat="identity", 
           color="#4fa3a4",
           fill = "#4fa3a4")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank())+
  labs(y=expression("Viral CPM ("~log[10]~')'),
       x="Viral rank")+
  theme(legend.position="none")
rank_DAD
ggsave("/Users/cissell/Desktop/RAC_DAD.svg",rank_DAD,width=6.5,height=3)


#Predicting phage-host linkages by combining multiple in-silico approaches
######
#Read in files
######
#Read in files
Death_crispr=read.table("/Users/cissell/Desktop/VIRSORTER2_Death/DA_blast_spacer.txt",header=T)
Death_nt=read.table("/Users/cissell/Desktop/VIRSORTER2_Death/DA_blast_nt_id.txt",header=T)
Death_dstar=read.csv("/Users/cissell/Desktop/VIRSORTER2_Death/DA_predictions.csv")

#Filter spacer by mismatch ≤1 and e-value ≤1, then best bitscore
Death_crispr_filter=Death_crispr %>%
  filter(mismatch<=1 & e_value <=1) %>%
  group_by(Phage) %>%
  filter(bitscore==max(bitscore)) %>%
  ungroup() %>%
  mutate(Source="crispr")

#Order by phage ID
Death_crispr_filter=Death_crispr_filter[order(Death_crispr_filter$Phage),]

#Filter nt by ≥50 bitscore and ≤0.001 e value
Death_nt_filter=Death_nt%>%
  filter(bitscore>=50 & e_value <=0.001)

#Extract only highest bitscore for each unique phage ID
Death_nt_high=Death_nt_filter%>%
  group_by(Phage) %>%
  filter(bitscore==max(bitscore)) %>%
  ungroup() %>%
  mutate(Source="nt")

#Order by phage ID
Death_nt_high=Death_nt_high[order(Death_nt_high$Phage),]

#Remove blanks from Kmer dataframe
Death_dstar=na.omit(Death_dstar)

#Find median k-mer matches number to filter by
median(Death_dstar$common.kmers)

#Filter d2star by common.kmer>305.5
Death_dstar_filter=Death_dstar%>%
  filter(common.kmers>305.5) %>%
  mutate(Source="kmer")

#Now create smaller dataframes from each for merging
Death_crisp_merger=Death_crispr_filter[,c(1,2,13)]
Death_nt_merger=Death_nt_high[,c(1,2,13)]
Death_dstar_merger=Death_dstar_filter[,c(1,2,6)]

#Now merge
Death_ph_merge1=merge(Death_crisp_merger,Death_nt_merger, by="Phage",all=T)
Death_ph_merge_2=merge(Death_ph_merge1,Death_dstar_merger,by="Phage",all=T)

#Now create new column with single best prediction
Death_ph=Death_ph_merge_2
Death_ph=dplyr::rename(Death_ph,MAG_cr=MAG.x)
Death_ph=dplyr::rename(Death_ph,MAG_nt=MAG.y)
Death_ph$MAG_cr=as.character(Death_ph$MAG_cr)
Death_ph$MAG_nt=as.character(Death_ph$MAG_nt)
Death_ph$MAG=as.character(Death_ph$MAG)
Death_ph=Death_ph%>%
  mutate(best=ifelse(is.na(MAG_cr) & is.na(MAG_nt),MAG[],
                     ifelse(is.na(MAG_cr),MAG_nt[],MAG_cr[])))

#Get length of unique phages to make sure all persist to end
length(unique(Death_ph$Phage))#352

##Eliminate redundancy based on conditions
#Create dataframe without duplicates for later combining
Death_ph_nd=Death_ph[!c(duplicated(Death_ph$Phage) | duplicated(Death_ph$Phage, fromLast=T)),]
#Check
duplicated(Death_ph_nd$Phage)
length(unique(Death_ph_nd$Phage))#302
#Create dataframe with all duplicates
Death_ph_d=Death_ph[duplicated(Death_ph$Phage) | duplicated(Death_ph$Phage, fromLast=T),]
#Check
length(unique(Death_ph_d$Phage))#50, adds to 352 so still good
#Focus dataframe on just meaningful columns
Death_ph_d=Death_ph_d[,-c(3,5:8)]
#Create dummy column and fill with NA
Death_ph_d$Con=NA
#Fill dummy column with Yes if crispr and nt match
Death_ph_d$Con[Death_ph_d$MAG_cr==Death_ph_d$MAG_nt] = "Yes"
#Fill dummy column with no if crisp and nt do not match
Death_ph_d$Con[Death_ph_d$MAG_cr!=Death_ph_d$MAG_nt] = "No"

#Extract Concencus dataframe for later combining
ConcensusDeath=Death_ph_d%>%
  filter(Con=="Yes")%>%
  mutate(best=MAG_cr)%>%
  select(Phage,best)%>%
  filter(!duplicated(Phage))

#5 of 50 recovered
#These  9 cases had a duplicate 'No' that we need to remove from the frame for downstream
Death_ph_d_non=Death_ph_d%>%
  group_by(Phage)%>%
  filter(!Phage%in%ConcensusDeath$Phage)


#Pipeline to remove those that had concencus and then if not bring down to a single column with redundant hits
Death_ph_dup=Death_ph_d_non%>%
  group_by(Phage)%>%
  filter(is.na(Con))%>%
  mutate(Con=ifelse(is.na(MAG_nt),MAG_cr[],MAG_nt[]))%>%
  select(-MAG_cr,-MAG_nt) %>%
  mutate(ID=row_number(),
         V="v")%>%
  unite(vid,c(V,ID),sep="")
length(unique(Death_ph_dup$Phage)) #40 of  50 recovered, that means No should have 1

#Pipeline to extract duplicates supported with multiple identical assignments and make nonredundant for later combining and to create vector for subsetting other dataframe
Death_ph_d2=Death_ph_dup%>%
  group_by(Phage)%>%
  filter(c(duplicated(Con)) | duplicated(Con, fromLast=T))%>%
  filter(duplicated(Phage))%>%
  mutate(best=Con)%>%
  select(-Con,-vid)%>%
  filter(!duplicated(Phage))
length(unique(Death_ph_d2$Phage)) #27

#Create vector for subsetting
rdeath=as.vector(Death_ph_d2$Phage)

#Pipeline to remove those already in above dataframe and then pick the first of remaining duplicates
Death_ph_nd_2=Death_ph_dup%>%
  group_by(Phage)%>%
  filter(!Phage%in%rdeath)%>%
  filter(!duplicated(Phage))%>%
  mutate(best=Con)%>%
  select(-Con,-vid)
length(unique(Death_ph_nd_2$Phage)) #13 recovered, Death_ph_nd_2+Death_ph_d2= recovered

#Pipeline to extract thos assigned 'No', use CRISPR as host, and remove redundancy
#First extract those supported multiple times by Crispr
NoncensusDeath=Death_ph_d_non%>%
  filter(Con=="No")%>%
  mutate(best=MAG_cr)%>%
  select(Phage,best) %>%
  mutate(ID=row_number(),
         V="v")%>%
  unite(vid,c(V,ID),sep="")
length(unique(NoncensusDeath$Phage))

NoncensusDeath_2=NoncensusDeath%>%
  group_by(Phage)%>%
  filter(c(duplicated(best)) | duplicated(best, fromLast=T)) %>%
  filter(duplicated(Phage))%>%
  mutate(best=best)%>%
  select(-vid)%>%
  filter(!duplicated(Phage))
length(unique(NoncensusDeath_2$Phage)) #5

#Create another vector for subsetting
rnoncenDeath=as.vector(NoncensusDeath_2$Phage)

#Pipeline to remove those already in above dataframe and then pick the first of remaining duplicates
NoncensusDeath_3=NoncensusDeath%>%
  group_by(Phage)%>%
  filter(!Phage%in%rnoncenDeath)%>%
  filter(!duplicated(Phage))%>%
  select(-vid)
length(unique(NoncensusDeath_3$Phage)) #0

#Make original dataframe fit the columnns
Death_ph_nd_o=Death_ph_nd%>%
  select(Phage,best)
#Create final dataframe
Death_phage_host_final=bind_rows(Death_ph_nd_o,
                                 Death_ph_d2,
                                 Death_ph_nd_2,
                                 ConcensusDeath,
                                 NoncensusDeath_2,
                                 NoncensusDeath_3)
#Make hosts factor
Death_phage_host_final$best=as.factor(Death_phage_host_final$best)
#Rename best to Host
Death_phage_host_final=dplyr::rename(Death_phage_host_final,Host=best)

#Ensure there are no duplicates
length(unique(Death_phage_host_final$Phage)) #352
duplicated(Death_phage_host_final$Phage)

########################
#Time to get the host and phage abundances going for some plots
########################
##############################################
#Read in host counts
Death_host_counts=read.xls("/Users/cissell/Desktop/VIRSORTER2_Death/genome_tpm.xlsx")
#Phage abundances were read in earlier
#death_viral_abunds
###NOW LETS BUILD THE PER MAT FRAMES WITH PAIRS,TAX,COUNTS,and LENGTHS and calculate TPM
death_viral_abunds2=death_viral_abunds%>%
  dplyr::rename(Phage=contig_id)
Death_tax_counts1=merge(death_viral_abunds2,Death_phage_host_final,by="Phage")
Death_tax_counts2=merge(Death_tax_counts1,Death_host_counts,by="Host",all.x=T)
#Make new df
Death_mag=Death_tax_counts2
str(Death_mag)
#Put length in KB
Death_mag$Plength=Death_mag$Plength/1000
Death_mag$Hlength=Death_mag$Hlength/1000
#Create RPK for Mat 1
Death_mag$PDAB=Death_mag$PDAB/Death_mag$Plength
Death_mag$HDAB=Death_mag$HDAB/Death_mag$Hlength
Death_mag$PDAD=Death_mag$PDAD/Death_mag$Plength
Death_mag$HDAD=Death_mag$HDAD/Death_mag$Hlength
Death_mag$PDAL=Death_mag$PDAL/Death_mag$Plength
Death_mag$HDAL=Death_mag$HDAL/Death_mag$Hlength
#Need to extract a new dataframe of only unique so each MAG only occurs once
Deathhpk1=Death_mag[!duplicated(Death_mag$Host),]
#create scaling factors
#Mat 1
Deathmsum11=((sum(Death_mag$PDAB)+sum(Deathhpk1$HDAB))/1000000)
Deathmsum13=((sum(Death_mag$PDAD)+sum(Deathhpk1$HDAD))/1000000)
Deathmsum15=((sum(Death_mag$PDAL)+sum(Deathhpk1$HDAL))/1000000)
#Create TPM for Mat 1
Death_mag$PDAB=Death_mag$PDAB/Deathmsum11
Death_mag$HDAB=Death_mag$HDAB/Deathmsum11
Death_mag$PDAD=Death_mag$PDAD/Deathmsum13
Death_mag$HDAD=Death_mag$HDAD/Deathmsum13
Death_mag$PDAL=Death_mag$PDAL/Deathmsum15
Death_mag$HDAL=Death_mag$HDAL/Deathmsum15
#Wont make sense just yet becuase of duplicates, but hold on a second there sport

#Bring it back to a sensical name now
Death_magcount_tpm=Death_mag

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Create dataframe to merge back after summarizing
Death_host_lost=Death_magcount_tpm%>%
  select(c(Host,HDAB,HDAD,HDAL))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
Death_hcount_tpm2=Death_magcount_tpm%>%
  select(-Phage,Plength,Hlength)%>%
  group_by(Host)%>%
  summarise_at(c("PDAB","PDAD","PDAL"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
Death_count_merger=merge(Death_hcount_tpm2,Death_host_lost,by="Host",all=T)
str(Death_count_merger)
#Check for sense in TPM just to be sure
sum(Death_count_merger$PDAB)+sum(Death_count_merger$HDAB)
sum(Death_count_merger$PDAD)+sum(Death_count_merger$HDAD)
sum(Death_count_merger$PDAL)+sum(Death_count_merger$HDAL)
#Nice

#Create each individual VMR and logVMR for each pair
Deathcount_vmrs=Death_count_merger%>%
  mutate(VMRB=PDAB/HDAB)%>%
  mutate(logVMRB=log10((PDAB+1)/HDAB))%>%
  mutate(VMRD=PDAD/HDAD)%>%
  mutate(logVMRD=log10((PDAD+1)/HDAD))%>%
  mutate(VMRL=PDAL/HDAL)%>%
  mutate(logVMRL=log10((PDAL+1)/HDAL))%>%
  dplyr::rename(P_DAB=PDAB,
                P_DAD=PDAD,
                P_DAL=PDAL,
                H_DAB=HDAB,
                H_DAD=HDAD,
                H_DAL=HDAL,
                VMR_DAB=VMRB,
                VMR_DAD=VMRD,
                VMR_DAL=VMRL,
                logVMR_DAB=logVMRB,
                logVMR_DAD=logVMRD,
                logVMR_DAL=logVMRL)

#Pivot longer for later merging and plotting
Death_count_vmr_long=Deathcount_vmrs%>%
  pivot_longer(
    cols=-c(Host),
    names_to=c(".value","Sample"),
    names_sep="_")
Death_count_vmr_long=as.data.frame(Death_count_vmr_long)
Death_count_vmr_long$Sample=as.factor(Death_count_vmr_long$Sample)
str(Death_count_vmr_long)
Death_count_vmr_long2=Death_count_vmr_long%>%
  pivot_longer(cols=c(P,H),
               names_to="Tax",
               values_to="TPM")%>%
  mutate(logTPM=log10(TPM))
#Nice

#Redo the order of Sample for good-good plotting
Death_count_vmr_long2=Death_count_vmr_long2%>%
  mutate(Sample=fct_relevel(Sample,
                            "DAL",
                            "DAB",
                            "DAD"))

#Now we plot
#VMR
death_virus_mag_vmr=ggplot(Death_count_vmr_long2, aes(y=logVMR,x=Sample))+
  geom_point(aes(color=Host),
             size=5)+
  geom_line(aes(group=Host,color=Host),
            size=2)+
  scale_color_manual(values=c("#588b8b","#c8553d"))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("VMR ("~log[10]~')'),
       x="Sample")+
  theme(legend.position = "none")
death_virus_mag_vmr
ggsave("/Users/cissell/Desktop/deathvmr.svg",death_virus_mag_vmr,width=3.25,height=4)

#TPM
death_virus_mag_tpm=ggplot(Death_count_vmr_long2, aes(y=logTPM,x=Sample))+
  geom_point(aes(color=Host,shape=Tax),
             size=5)+
  geom_line(aes(group=interaction(Host,Tax),color=Host),
            size=2)+
  scale_color_manual(values=c("#588b8b","#c8553d"))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("CPM ("~log[10]~')'),
       x="Sample")+
  theme(legend.position = "none")+
  scale_y_continuous(breaks=seq(3.5,5.5,0.5))

death_virus_mag_tpm
ggsave("/Users/cissell/Desktop/deathtpm2.svg",death_virus_mag_tpm,width=3.25,height=4)



##########################
##RNA Data
##########################
#First let's look at expression of interesting physiological pathways in dominant cMAGs
#We will have to do some merging to get the feature count matrix
#Read in all the counts
rna_death=read.xls('/Users/cissell/Desktop/VIRSORTER2_Death/RNA/DeathMAGs_RNA_counts.xlsx')
#Read in the gene id
rna_deathkegg=read.xls('/Users/cissell/Desktop/VIRSORTER2_Death/RNA/Death_annotations.xlsx')

#First calculate TPM in RNA including all alignments
#Make new df
rpkdeath=rna_death
str(rpkdeath)
#Put length in KB
rpkdeath$Length=rpkdeath$Length/1000
#Create RPK
rpkdeath$RAB=rpkdeath$RAB/rpkdeath$Length
rpkdeath$RAD=rpkdeath$RAD/rpkdeath$Length
#Need to extract a new dataframe of only unique so each MAG only occurs once
deathhpkmar=rpkdeath[!duplicated(rpkdeath$Host),]
#create scaling factors
#Mat 1
deathsumRAB=(sum(rpkdeath$RAB)/1000000)
deathsumRAD=(sum(rpkdeath$RAD)/1000000)
#Create TPM for Mat 1
rpkdeath$RAB=rpkdeath$RAB/deathsumRAB                     
rpkdeath$RAD=rpkdeath$RAD/deathsumRAD

#Check for sense
sum(rpkdeath$RAB) #Good
sum(rpkdeath$RAD) #Good

#Merge counts with gene ids
death_RNA_merged=merge(rpkdeath,rna_deathkegg, by="contig",all.x=T)

#Remove rows without KEGG ID
deathsmerged=death_RNA_merged %>%
  na_if("") %>%
  na.omit
#Drop contig
deathsmerged=subset(deathsmerged,select=-c(contig,Length))
#Summarize within KEGG ID within MAG
deathsummerged=deathsmerged1 %>%
  group_by(kegg_id,Host) %>%
  dplyr::summarise(across(everything(), sum))

#Apply an abundance normalization (from DNA) to TPM data
#First grab abundance data from above DNA
Death_abund_normal_matrix = Death_count_vmr_long2 %>%
  select(c(Host,Sample,Tax,TPM)) %>%
  subset(Tax=="H") %>%
  select(-Tax) %>%
  subset(Sample=="DAB" | Sample=="DAD")
#DAB_MAG_2 Normalization abundances
#DAB - 25933.953
#DAD - 14153.746
#DAL_MAG_7 Normalization abundances
#DAB - 3643.213
#DAD - 2492.928
#Now apply normalization
Death_abund_normalizedDAB = deathsummerged %>%
  subset(Host=="DAB_2") %>%
  mutate(RAB = RAB/25933.953) %>%
  mutate(RAD = RAD/14153.746)

Death_abund_normalizedDAL = deathsummerged %>%
  subset(Host=="DAL_7") %>%
  mutate(RAB = RAB/3643.213) %>%
  mutate(RAD = RAD/2492.928)

Death_abund_normalized=rbind(Death_abund_normalizedDAB,Death_abund_normalizedDAL)
#Manually subset to KEGG regions I am interested in (pull from Cisell & McCoy In Review Molecular Ecology)
Death_abund_normal_subset=Death_abund_normalized %>%
  subset(kegg_id=="K02703" |
           kegg_id=="K02706" |
           kegg_id=="K02705" |
           kegg_id=="K02689" |
           kegg_id=="K02690" |
           kegg_id=="K02691" |
           kegg_id=="K02692")

#Now smaller subsets for plotting
#Photosynthesis
Death_abund_normal_photosubset=Death_abund_normal_subset%>%
  pivot_longer(cols=c(RAD,RAB),names_to="sample",values_to="Abund")

#plot
cmag_photosynth_tpm=ggplot(Death_abund_normal_photosubset, aes(y=log10(Abund),x=sample,fill=Host))+
  geom_boxplot(size=1.5)+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Abundance normalized CPM ("~log[10]~')'),
       x="Sample") +
  scale_fill_manual(values=c("#588b8b","#c8553d"))
#  theme(legend.position = "none")=
cmag_photosynth_tpm
ggsave("/Users/cissell/Desktop/deathphotosynthtpm.svg",cmag_photosynth_tpm,width=5.5,height=5)

###############
####Now do global phage RNA abundance across transect and VMAR analysis for cMAG:Phage pairs
#Read in Phage RNA counts
death_RNA_phages=read.xls("/Users/cissell/Desktop/VIRSORTER2_Death/RNA/death_virus_counts.xlsx")

#First we need to summarize ORF length and mapping to get single length and RNA abund per vOTU
death_RNA_phages_summarized=death_RNA_phages %>%
  select(-contig_id) %>%
  group_by(Phage) %>%
  mutate(Length=sum(Length)) %>%
  mutate(RAB=sum(RAB)) %>%
  mutate(RAD=sum(RAD)) %>%
  ungroup() %>%
  distinct(Phage, .keep_all=T)
#Make new df
rnadeath_phage=death_RNA_phages_summarized
str(rnadeath_phage)
#Put length in KB
rnadeath_phage$Length=rnadeath_phage$Length/1000
#Create RPK
rnadeath_phage$RAB=rnadeath_phage$RAB/rnadeath_phage$Length
rnadeath_phage$RAD=rnadeath_phage$RAD/rnadeath_phage$Length
#create scaling factors
sumRABpo=(sum(rnadeath_phage$RAB)/1000000)
sumRADpo=(sum(rnadeath_phage$RAD)/1000000)
#Create TPM
rnadeath_phage$RAB=rnadeath_phage$RAB/sumRABpo
rnadeath_phage$RAD=rnadeath_phage$RAD/sumRADpo

#Check for sense
sum(rnadeath_phage$RAB) #Good
sum(rnadeath_phage$RAD) #Good


#Pivot longer for plotting and do log
Death_Rtpm_phageonly_long=rnadeath_phage%>%
  select(-Length)%>%
  pivot_longer(
    cols=c(RAB,RAD),
    names_to="Sample",
    values_to="TPM")%>%
  mutate(logTPM=log10(TPM+1))
Death_Rtpm_phageonly_long$Sample=as.factor(Death_Rtpm_phageonly_long$Sample)

#Add RAL as a blank level for good spacing on the plot for combining with DNA
Death_Rtpm_phageonly_long$Sample = factor(Death_Rtpm_phageonly_long$Sample,
                                          levels = c(levels(Death_Rtpm_phageonly_long$Sample),"RAL"))

#Redo the order of Sample for good-good plotting
Death_Rtpm_phageonly_long=Death_Rtpm_phageonly_long%>%
  mutate(Sample=fct_relevel(Sample,
                            "RAL",
                            "RAB",
                            "RAD"))
#Find median numbers
aggregate(logTPM~Sample,data=Death_Rtpm_phageonly_long,summary)
#Median DAB = 1.70
#Median DAD = 1.59

#Now violin
#with zeros
RNAdeath_virusviolin=ggplot(Death_Rtpm_phageonly_long, aes(y=logTPM,x=Sample))+
  geom_violin(aes(fill=Sample),
              outlier.shape = NA,
              size=1,
              alpha=0.9,
              draw_quantiles=c(0.5))+
  scale_fill_manual(values=c("black","#fbaf37","#4fa3a4"),name="Time",
                    drop=FALSE)+
  geom_point(aes(fill=Sample),position=position_jitterdodge(0.2),colour="black",pch=21,
             alpha=0.7,
             size=1.5)+
  geom_violin(aes(fill=Sample),
              outlier.shape = NA,
              size=1,
              alpha=0,
              draw_quantiles=c(0.25,0.5,0.75))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Viral TPM ("~log[10]~')'),
       x="Sample")+
  theme(legend.position="none") +
  scale_x_discrete(drop=FALSE)
RNAdeath_virusviolin
ggsave("/Users/cissell/Desktop/viral_props_deathRNA.svg",RNAdeath_virusviolin,width=6,height=3)

#Now do RNA with lineage resolved pairs
#death_viral_abunds
###NOW LETS BUILD THE PER MAG FRAMES WITH PAIRS,TAX,COUNTS,and LENGTHS and calculate TPM
Death_tax_countsRNA1=merge(death_RNA_phages_summarized,Death_phage_host_final,by="Phage")
#Add distinction for Phage (P)
Death_tax_countsRNA1=Death_tax_countsRNA1%>%
  dplyr::rename(PRAB=RAB) %>%
  dplyr::rename(PRAD=RAD) %>%
  dplyr::rename(PLength=Length)
#Summarize host and add distinction for Host (H)
deathhost_rna_counts = rna_death %>%
  select(-contig) %>%
  group_by(Host)%>%
  mutate(RAB=sum(RAB)) %>%
  mutate(RAD=sum(RAD)) %>%
  mutate(Length=sum(Length)) %>%
  ungroup()%>%
  dplyr::rename(HRAB=RAB) %>%
  dplyr::rename(HRAD=RAD) %>%
  dplyr::rename(HLength=Length) %>%
  distinct(Host, .keep_all=T)
levels (deathhost_rna_counts$Host) = c("DAB_MAG_2","DAL_MAG_7")
#Merge
Death_tax_countsRNA2=merge(Death_tax_countsRNA1,deathhost_rna_counts,by="Host",all.x=T)
#Make new df
Death_rna_mag=Death_tax_countsRNA2
str(Death_rna_mag)
#Put length in KB
Death_rna_mag$PLength=Death_rna_mag$PLength/1000
Death_rna_mag$HLength=Death_rna_mag$HLength/1000
#Create RPK for Mat 1
Death_rna_mag$PRAB=Death_rna_mag$PRAB/Death_rna_mag$PLength
Death_rna_mag$HRAB=Death_rna_mag$HRAB/Death_rna_mag$HLength
Death_rna_mag$PRAD=Death_rna_mag$PRAD/Death_rna_mag$PLength
Death_rna_mag$HRAD=Death_rna_mag$HRAD/Death_rna_mag$HLength

#Need to extract a new dataframe of only unique so each MAG only occurs once
Deathhpk1rna=Death_rna_mag[!duplicated(Death_rna_mag$Host),]
#create scaling factors
#Mat 1
Deathmsum11rnaph=((sum(Death_rna_mag$PRAB)+sum(Deathhpk1rna$HRAB))/1000000)
Deathmsum13rnaph=((sum(Death_rna_mag$PRAD)+sum(Deathhpk1rna$HRAD))/1000000)

#Create TPM
Death_rna_mag$PRAB=Death_rna_mag$PRAB/Deathmsum11rnaph
Death_rna_mag$HRAB=Death_rna_mag$HRAB/Deathmsum11rnaph
Death_rna_mag$PRAD=Death_rna_mag$PRAD/Deathmsum13rnaph
Death_rna_mag$HRAD=Death_rna_mag$HRAD/Deathmsum13rnaph

#Wont make sense just yet becuase of duplicates, but hold on a second there sport

#Bring it back to a sensical name now
Death_magcount_tpm_rnaph=Death_rna_mag

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Create dataframe to merge back after summarizing
Death_host_lost_rna=Death_magcount_tpm_rnaph%>%
  select(c(Host,HRAB,HRAD))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
Death_magcount_tpm_rnaph2=Death_magcount_tpm_rnaph%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("PRAB","PRAD"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
Death_count_merger_rnaph=merge(Death_magcount_tpm_rnaph2,Death_host_lost_rna,by="Host",all=T)
str(Death_count_merger_rnaph)
#Check for sense in TPM just to be sure
sum(Death_count_merger_rnaph$PRAB)+sum(Death_count_merger_rnaph$HRAB)
sum(Death_count_merger_rnaph$PRAD)+sum(Death_count_merger_rnaph$HRAD)
#Nice

#Create each individual VMR and logVMR for each pair
Deathcount_vmrs_rnaph=Death_count_merger_rnaph%>%
  mutate(VMRB=PRAB/HRAB)%>%
  mutate(logVMRB=log10((PRAB+1)/HRAB))%>%
  mutate(VMRD=PRAD/HRAD)%>%
  mutate(logVMRD=log10((PRAD+1)/HRAD))%>%
  dplyr::rename(P_RAB=PRAB,
                P_RAD=PRAD,
                H_RAB=HRAB,
                H_RAD=HRAD,
                VMR_RAB=VMRB,
                VMR_RAD=VMRD,
                logVMR_RAB=logVMRB,
                logVMR_RAD=logVMRD)

#Pivot longer for later merging and plotting
Death_count_vmr_long_rnaph=Deathcount_vmrs_rnaph%>%
  pivot_longer(
    cols=-c(Host),
    names_to=c(".value","Sample"),
    names_sep="_")
Death_count_vmr_long_rnaph=as.data.frame(Death_count_vmr_long_rnaph)
Death_count_vmr_long_rnaph$Sample=as.factor(Death_count_vmr_long_rnaph$Sample)
str(Death_count_vmr_long_rnaph)
Death_count_vmr_long_rnaph2=Death_count_vmr_long_rnaph%>%
  pivot_longer(cols=c(P,H),
               names_to="Tax",
               values_to="TPM")%>%
  mutate(logTPM=log10(TPM))
#Nice

#Redo the order of Sample for good-good plotting
Death_count_vmr_long_rnaph2=Death_count_vmr_long_rnaph2%>%
  mutate(Sample=fct_relevel(Sample,
                            "RAB",
                            "RAD"))

#Now we plot
#VMR
death_virus_mag_vmr_rnaph=ggplot(Death_count_vmr_long_rnaph2, aes(y=logVMR,x=Sample))+
  geom_point(aes(color=Host),
             size=5)+
  geom_line(aes(group=Host,color=Host),
            size=2)+
  scale_color_manual(values=c("#588b8b","#c8553d"))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("VMAR ("~log[10]~')'),
       x="Sample")+
  theme(legend.position = "none")
death_virus_mag_vmr_rnaph
ggsave("/Users/cissell/Desktop/deathvmr_rnaph.svg",death_virus_mag_vmr_rnaph,width=3.25,height=4)

#TPM
death_virus_mag_tpm_rnaph=ggplot(Death_count_vmr_long_rnaph2, aes(y=logTPM,x=Sample))+
  geom_point(aes(color=Host,shape=Tax),
             size=5)+
  geom_line(aes(group=interaction(Host,Tax),color=Host),
            size=2)+
  scale_color_manual(values=c("#588b8b","#c8553d"))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("TPM ("~log[10]~')'),
       x="Sample")+
  theme(legend.position = "none")

death_virus_mag_tpm_rnaph
ggsave("/Users/cissell/Desktop/deathtpm2_rnaph.svg",death_virus_mag_tpm_rnaph,width=3.25,height=4)
