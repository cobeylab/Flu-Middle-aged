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

fig_dir = "../figure/"
tb_dir = "../table/"
age_group_1 = c(0, 2, 10, 20, 30, 40, 50, 60, 70, 90)


breaks1 = seq(0, 10, 2)
labels1 = 2^(breaks1)*10

breaks2 = seq(0, 10, 2)
labels2 = 2^breaks2

#############################################################################
# Load data
#############################################################################

# sera
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



# covariates
xvarf = "../../Data/covariate/imprinting_covariates_FRNT_HA_viruses.csv"
imp = read.table(xvarf, head=T, sep=",", stringsAsFactors = F )

imp$Age = 2017 - imp$y
#imp$Virus = imp$virus
imp$Virus = ifelse(imp$Virus == "WT", "3c2.A", imp$Virus)
imp$Virus = ifelse(imp$Virus == "T131K_R142K_R261Q", "A2", imp$Virus)
#unique(imp$Virus)


#Merge input data
df = merge(msera, imp, by = c("Virus", "Age"))
df = transform(df, Age_group = cut(Age, breaks = c(0, 4, 17, 44, 64, 90)))
df = transform(df, Age_group2 = cut(Age, breaks = c(0, 64, 90)))
df$std_Age = (df$Age - mean(df$Age))/sd(df$Age)
df$detectable = ifelse(df$l2FRNT == 0, 0, 1)
df$similarity = df$episim

#CI
cif = "../../Data/covariate/subtype_imprinting_probabilities_CI.csv"
ci = read.table(cif, head=T, sep=",", stringsAsFactors = F)
head(ci)
ci$Age = 2017-ci$yob
ci_tv1 = ci[,c("Age", "similarity_tv1", "similarity_tv1_ci1", "similarity_tv1_ci2")]
ci_tv1$Virus = "3c2.A"
colnames(ci_tv1) = c("Age", "Similarity_rmvn", "Similarity_ci1", "Similarity_ci2", "Virus")
ci_tv2 = ci[,c("Age", "similarity_tv2", "similarity_tv2_ci1", "similarity_tv2_ci2")]
ci_tv2$Virus = "A2"
colnames(ci_tv2) = c("Age", "Similarity_rmvn", "Similarity_ci1", "Similarity_ci2", "Virus")

ci_similarity = rbind(ci_tv1, ci_tv2)

ci_imp = ci[,c("Age", "imp_H3N2_ci1", "imp_H3N2_ci2")]

df = merge(df, ci_similarity, by=c("Age", "Virus"))
df = merge(df, ci_imp, by="Age")

p1 = ggplot(df) +
  geom_smooth(aes(x=Age, y=l2FRNT, col=Virus), span=0.6) +
  scale_y_continuous(breaks = c(0,1,2,3), labels = 2^c(0,1,2,3)*10, name = "FRNT") +
  scale_color_manual(labels = c("3c2.A",
                                "3c2.A2"),
                     values = c("#006ddb", "#920000"),
                     name = "") +
  theme_classic() +
  theme(legend.position = c(0.7, 0.9),legend.background = element_rect(fill=NA)) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 2.0), "lines")) +
  theme(axis.title.x = element_blank())


p2 = ggplot(df) +
  geom_line(aes(x=Age, y=imp_h3n2)) +
  geom_ribbon(aes(x=Age, ymin = imp_H3N2_ci1, ymax = imp_H3N2_ci2), alpha=0.55) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  theme_classic() +
  ylab("Imprinting probability to H3") +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 1.6), "lines")) +
  theme(axis.title.x = element_blank())
  
p12 = ggplot(df) +
  geom_smooth(aes(x=Age, y=l2FRNT, col=Virus), span=0.6) +
  geom_line(aes(x=Age, y=imp_h3n2*3, col="Imprinting")) +
  geom_ribbon(aes(x=Age, ymin = imp_H3N2_ci1*3, ymax = imp_H3N2_ci2*3), alpha=0.55) +
  scale_y_continuous(breaks = c(0,1,2,3), labels = 2^c(0,1,2,3)*10, name = "FRNT",
                     sec.axis = sec_axis(~./3, 
                                         name = "Imprinting probability to H3", 
                                         breaks= c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
                                         labels = c(0, 0.2, 0.4, 0.6, 0.8, 1.0) )) +
  scale_color_manual(labels = c("3c2.A",
                                "3c2.A2",
                                "Imprinting to H3"),
                     values = c("#006ddb", "#920000", "#000000"),
                     name = "") +
  theme_classic() +
  theme(plot.margin = unit(c(0.2, 0.9, 0.2, 1.3), "lines")) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = c(0.7, 0.9),legend.background = element_rect(fill=NA)) 
  

p3 = ggplot(df) +
  geom_line(aes(x=Age, y=Similarity_rmvn, col=Virus)) +
  geom_ribbon(aes(x=Age, ymin=Similarity_ci1, ymax=Similarity_ci2, fill=Virus), alpha=0.5 ) +
  
  scale_color_manual(labels = c("Similarity between imprinting H3 and 3c2.A HA",
                                  "Similarity between imprinting H3 and 3c2.A2 HA"),
                     values = c("#006ddb", "#920000"),
                      name = "") +
  scale_fill_manual(labels = c("Similarity between imprinting H3 and 3c2.A HA",
                                "Similarity between imprinting H3 and 3c2.A2 HA"),
                     values = c("#006ddb", "#920000"),
                     name = "") +
  theme_classic() +
  theme(legend.position = c(0.65, 0.9),legend.background = element_rect(fill=NA)) +
  theme(plot.margin = unit(c(0.2, 3.1, 0.2, 0.95), "lines")) +
  ylab("Amino-acid similarity")
  

library(ggpubr)

ggarrange(p12,  p3,
          nrow=2,
          heights = c(3,2))


ggsave("../figure/Fig3B_Imprinting to H3.png", width=6, height=5.5)
ggsave("../figure/Fig3B_Imprinting to H3.pdf", width=6, height=5.5)


##########################################################################
df$vac_cov = 0.01*df$vac_cov
df$children = ifelse(df$Age < 18, 1, 0)

#simple check
fit = glm( df$l2FRNT ~  
                   df$imp_h3n2 + 
                   df$Similarity_rmvn + 
                   df$Virus +
                   df$vac_cov )


df$res = fit$residuals
df$fitted = fit$fitted.values

ggplot(df) +
  geom_point(aes(x=Age, y=res)) +
  geom_smooth(aes(x=Age, y=res), span=0.3) +
  geom_smooth(aes(x=Age, y=fitted), span=0.3, col="red")

#####################################
# regression using felm to cluster error by individual


#Imprinting and similarity

fit_all_simil= glm( df$l2FRNT ~  
      df$imp_h3n2 + 
        df$similarity + 
      df$Virus +
      df$vac_cov)

# Imprinting, similarity, and age

fit_all_simil_a = glm(df$l2FRNT ~ 
                  df$imp_h3n2 + 
                    df$similarity + 
                  df$Virus +
                  df$vac_cov + 
                  df$Age)


fit_all_simil_ch = glm(df$l2FRNT ~ 
                  df$imp_h3n2 + 
                    df$similarity + 
                  df$Virus +
                  df$vac_cov + 
                  df$children)



summary(fit_all_simil)
summary(fit_all_simil_a)
summary(fit_all_simil_ch)


