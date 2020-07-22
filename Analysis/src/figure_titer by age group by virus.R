library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(Hmisc)
library(ggsignif)
library(lmtest)
library(multiwayvcov)

#################################################

grouping = 1  

age_group_1 = c(0, 2, 10, 20, 30, 40, 50, 60, 70, 90)


breaks1 = seq(0, 10, 2)
labels1 = 2^(breaks1)*10

breaks2 = seq(0, 10, 2)
labels2 = 2^breaks2

breaks_elisa = seq(0, 5)
labels_elisa = 2^(breaks_elisa)

breaks_frnt = seq(0, 4)
labels_frnt = 2^(breaks_frnt)*10

breaks_ella = seq(0, 5)
labels_ella = 2^(breaks_ella)*10


#############################################################################
# Load data
#############################################################################

sera = read.table("../../Data/titer/titer.csv", sep=",", head=T, stringsAsFactors = F)
frnt = sera[,c(1,2,4,5)]
elisa = sera[,c(1,2,6,7)]
ella = sera[,c(1,2,8,9)]

colnames(frnt) = gsub("FRNT_", "", colnames(frnt))
frnt = melt(frnt, id.vars = c("Sample_ID", "Age"))

colnames(elisa) = gsub("ELISA_HA_", "", colnames(elisa))
elisa = melt(elisa, id.vars = c("Sample_ID", "Age"))

colnames(ella) = gsub("ELLA_", "", colnames(ella))
ella = melt(ella, id.vars = c("Sample_ID", "Age"))

msera = merge(frnt, elisa, by=c("Sample_ID", "Age", "variable"))
msera = merge(msera, ella, by=c("Sample_ID", "Age", "variable"))

colnames(msera) = c("Sample_ID", "Age", "Virus", "FRNT", "ELISA", "ELLA")
msera$l2FRNT = log2(msera$FRNT/10)
msera$l2ELISA = log2(msera$ELISA)
msera$l2ELLA = log2(msera$ELLA/10)

msera = transform(msera, Age_group = cut(Age, age_group_1))
msera$Age_group = gsub("\\(", "", msera$Age_group)
msera$Age_group = gsub("\\]", "", msera$Age_group)

cc = strsplit(msera$Age_group, ",")
part1 = unlist(cc)[2*(1:nrow(msera))-1]
part2 = unlist(cc)[2*(1:nrow(msera)) ]
part1 = as.numeric(part1)+1
msera$Age_group = paste0(part1, "-", part2)
msera$Age_group = factor(msera$Age_group, 
                        levels = c("1-2", "3-10", "11-20", "21-30", "31-40", 
                                   "41-50", "51-60", "61-70", "71-90"))

################################################
#Fig 1

element_text_size = 13
xaxis_size = 12

virus = "3c2.A"
msera_sub = msera[msera$Virus == virus,]


fig1a = ggplot(msera_sub, aes(x=Age_group, y=l2ELISA)) +
  geom_hline(yintercept=log2(0.7), linetype="dashed") +
  geom_jitter(width = 0.2, height=0, alpha=0.2) +
  stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width=0.5) +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",
               size=1, width=0.4) +
  scale_y_continuous(breaks=breaks2, labels = labels2) +
  xlab("Age group (years)") +
  ylab("Relative IgG concentration") +
  coord_cartesian(ylim=c(-1, 11)) +
  theme_classic() +
  theme(text = element_text(size=element_text_size),
        axis.text.x = element_text(size=xaxis_size, angle=90),
        axis.text.y = element_text(size = xaxis_size)) 
  

fig1b = ggplot(msera_sub, aes(x=Age_group, y=l2FRNT)) +
  geom_hline(yintercept=log2(20/10), linetype="dashed") +
  geom_jitter(width = 0.2, height=0, alpha=0.2) +
  stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width=0.5) +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",
               size=1, width=0.4) +
  scale_y_continuous(breaks=breaks1, labels = labels1) +
  xlab("Age group (years)") +
  ylab("FRNT") +
  coord_cartesian(ylim=c(0, 8)) +
  theme_classic() +
  theme(text = element_text(size=element_text_size),
        axis.text.x = element_text(size=xaxis_size, angle=90),
        axis.text.y = element_text(size = xaxis_size)) 
  


fig1c = ggplot(msera_sub, aes(x=Age_group, y=l2ELLA)) +
  geom_hline(yintercept=log2(20/10), linetype="dashed") +
  geom_jitter(width = 0.2, height=0, alpha=0.2) +
  stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width=0.5) +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",
               size=1, width=0.4) +
  scale_y_continuous(breaks=breaks1, labels = labels1) +
  xlab("Age group (years)") +
  ylab("ELLA") +
  coord_cartesian(ylim=c(0, 8))+
  theme_classic() +
  theme(text = element_text(size=element_text_size),
        axis.text.x = element_text(size=xaxis_size, angle=90),
        axis.text.y = element_text(size = xaxis_size)) 
  

################################################################
virus = "A2"
msera_sub = msera[msera$Virus == virus,]

fig1d = ggplot(msera_sub, aes(x=Age_group, y=l2ELISA)) +
  geom_hline(yintercept=log2(0.7), linetype="dashed") +
  geom_jitter(width = 0.2, height=0, alpha=0.2) +
  stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width=0.5) +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",
               size=1, width=0.4) +
  scale_y_continuous(breaks=breaks2, labels = labels2) +
  xlab("Age group (years)") +
  ylab("Relative IgG concentration") +
  coord_cartesian(ylim=c(-1, 11)) +
  theme_classic()+
  theme(text = element_text(size=element_text_size),
        axis.text.x = element_text(size=xaxis_size, angle=90),
        axis.text.y = element_text(size = xaxis_size)) 

fig1e = ggplot(msera_sub, aes(x=Age_group, y=l2FRNT)) +
  geom_hline(yintercept=log2(20/10), linetype="dashed") +
  geom_jitter(width = 0.2, height=0, alpha=0.2) +
  stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width=0.5) +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",
               size=1, width=0.4) +
  scale_y_continuous(breaks=breaks1, labels = labels1) +
  xlab("Age group (years)") +
  ylab("FRNT") +
  coord_cartesian(ylim=c(0, 8)) +
  theme_classic() +
  theme(text = element_text(size=element_text_size),
        axis.text.x = element_text(size=xaxis_size, angle=90),
        axis.text.y = element_text(size = xaxis_size)) 


fig1f = ggplot(msera_sub, aes(x=Age_group, y=l2ELLA)) +
  geom_hline(yintercept=log2(20/10), linetype="dashed") +
  geom_jitter(width = 0.2, height=0, alpha=0.2) +
  stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width=0.5) +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar",
               size=1, width=0.4) +
  scale_y_continuous(breaks=breaks1, labels = labels1) +
  xlab("Age group (years)") +
  ylab("ELLA") +
  coord_cartesian(ylim=c(0, 8))+
  theme_classic()+
  theme(text = element_text(size=element_text_size),
        axis.text.x = element_text(size=xaxis_size, angle=90),
        axis.text.y = element_text(size = xaxis_size)) 


##############################################################
fig1g = ggplot(msera) +
  geom_smooth(aes(x=Age, y=l2ELISA, col=Virus)) +
  scale_y_continuous(breaks=breaks_elisa, labels=labels_elisa) +
  coord_cartesian(ylim=c(0, 5)) +
  ylab("Relative IgG concentration") +
  xlab("Age (years)") +
  scale_color_manual(guide = F, values = c("#006ddb", "#920000"))+
  theme_classic() +
  theme(text = element_text(size=element_text_size),
        axis.text.x = element_text(size=xaxis_size),
        legend.position = "non") 
  

fig1h = ggplot(msera) +
  geom_smooth(aes(x=Age, y=l2FRNT, col=Virus)) +
  scale_y_continuous(breaks=breaks_frnt, labels = labels_frnt) +
  coord_cartesian(ylim=c(0, 4)) +
  scale_color_manual(guide = F, values = c("#006ddb", "#920000")) +
  ylab("FRNT") +
  xlab("Age (years)") +
  theme_classic() +
  theme(text = element_text(size=element_text_size),
        axis.text.x = element_text(size=xaxis_size),
        legend.position = "non") 
   


fig1i = ggplot(msera) +
  geom_smooth(aes(x=Age, y=l2ELLA, col=Virus)) +
  scale_y_continuous(breaks=breaks_ella, labels=labels_ella) +
  coord_cartesian(ylim=c(0, 5)) + 
  scale_color_manual(name = "", values = c("#006ddb", "#920000"),
                     labels = c("3c2.A", "3c2.A2")) +
  ylab("ELLA") +
  xlab("Age (years)") +
  theme_classic() +
  theme(text = element_text(size=element_text_size),
        axis.text.x = element_text(size=xaxis_size),
        legend.position = "non") +
  
  theme(legend.position = c(0.75, 0.9),
        legend.background = element_rect(fill=NA))

  
  
return (ggarrange(fig1a, fig1b, fig1c, 
                  fig1d, fig1e, fig1f,
                  fig1g, fig1h, fig1i,
                  ncol = 3, nrow = 3, 
                  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I")))
  

ggsave(paste("../figure/Figure1_A-I.png"), width=10, height=10)
ggsave(paste("../figure/Figure1_A-I.pdf"), width=10, height=10)


