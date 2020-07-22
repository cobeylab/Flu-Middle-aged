library(reshape2)
library(ggplot2)
library(dplyr)
library(Hmisc)
library(ggsignif)
library(lmtest)
library(multiwayvcov)

source("utils.R")
source("boot_meanTiter.R")

#################################################

assay = "HA FRNT"

width_o=0.2
height_o=0.2

width_v = 0.2
height_v = 0.2

width_a = 0.2
heith_a = 0.2

width_va = width_a
height_va = heith_a

width_s = width_v
height_s = height_v


grouping = 1 

if(grouping == 1){
  fig_dir = "../extended_figure/frnt/"
  tb_dir = "../extended_table/frnt/"
  age_group_1 = c(0,2,10,20,30,40,50,60,70,90)
  age_group_2 = c(0, 17, 40, 60, 90)
}

#############################################################################
# Load data
#############################################################################

sera = read.table("../../Data/titer/titer.csv", sep=",", head=T, stringsAsFactors = F)

m_sera = melt(sera, id.vars = c("Sample_ID", "Age"))
colnames(m_sera)[grep("variable", colnames(m_sera))] = c("Test")
colnames(m_sera)[grep("value", colnames(m_sera))] = c("Titer")

# FRNT: 3c2a v. a2 titer
frnt = m_sera[grep("FRNT", m_sera$Test),]
frnt$Test = gsub("FRNT_", "", frnt$Test)
frnt$Titer = log2(frnt$Titer/10)
frnt = frnt[frnt$Test != "3c1",]

breaks = seq(0, 10, 2)
labels = 2^(breaks)*10

#Plot overview
ggplot(frnt, aes(x=Age, y=Titer)) +
  geom_smooth(aes(col=Test)) +
  geom_jitter(aes(col=Test), width=0, height=0.2) +
  theme(text = element_text(size=18)) +
  scale_color_discrete(name="") +
  scale_y_continuous(breaks=breaks, labels = labels) +
  ggtitle(assay)
#ggsave(paste(fig_dir, "overview.png"), width=6, height=4)

#assign age group
frnt = transform( frnt, Age_group = cut( Age, age_group_1 ) )
frnt$Age_group = gsub("\\(", "", frnt$Age_group)
frnt$Age_group = gsub("\\]", "", frnt$Age_group)

cc = strsplit(frnt$Age_group, ",")
part1 = unlist(cc)[2*(1:nrow(frnt))-1]
part2 = unlist(cc)[2*(1:nrow(frnt)) ]
part1 = as.numeric(part1)+1
frnt$Age_group = paste0(part1, "-", part2)

if("3-10" %in% frnt$Age_group){
  frnt$Age_group = factor(frnt$Age_group, 
                          levels = c("1-2", "3-10", "11-20", "21-30", "31-40", 
                                     "41-50", "51-60", "61-70", "71-90"))
}

###########################################################################
# Resample 
###########################################################################

#make wide data frame
frnt_w = reshape(frnt,
                 timevar = "Test",
                 idvar = c("Sample_ID", "Age", "Age_group"),
                 direction = "wide")

#make dataframe with new age group assignment

frnt2 = frnt[frnt$Sex != "",]

frnt2 = transform( frnt2, Age_group = cut( Age, age_group_2) )
frnt2$Age_group = gsub("\\(", "", frnt2$Age_group)
frnt2$Age_group = gsub("\\]", "", frnt2$Age_group)
cc = strsplit(frnt2$Age_group, ",")
part1 = unlist(cc)[2*(1:nrow(frnt2))-1]
part2 = unlist(cc)[2*(1:nrow(frnt2)) ]
part1 = as.numeric(part1)+1
frnt2$Age_group = paste0(part1, "-", part2)

frnt2_w = reshape(frnt2,
                  timevar = "Test",
                  idvar = c("Sample_ID", "Age", "Age_group"),
                  direction = "wide")

##############################################################

data = frnt
data2 = frnt2
data_w = frnt_w
data2_w = frnt2_w

source("analysis_for_assay.R")

