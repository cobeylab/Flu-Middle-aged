
library(reshape2)
library(ggplot2)
library(dplyr)
library(Hmisc)
library(ggsignif)
library(lmtest)
library(multiwayvcov)

source("utils.R")
source("boot_meanTiter.R")

assay = "HA ELISA"



width_o=0
height_o=0

width_v = 0.2
height_v = 0

width_a = 0.2
heith_a = 0.2

width_va = width_a
height_va = heith_a

width_s = width_v
height_s = height_v

grouping = 1

if(grouping == 1){
  fig_dir = "../figure/elisa/"
  tb_dir = "../extended_table/elisa/"
  age_group_1 = c(0,2,10,20,30,40,50,60,70,90)
  age_group_2 = c(0, 17, 40, 60, 90)
}
###############################################################################################

sera = read.table("../../Data/titer/titer.csv", sep=",", head=T, stringsAsFactors = F)

# elisa: 3c2a v. a2 titer
elisa = data.frame(sera[, c("Sample_ID", "Age")], sera[, grep("ELISA", colnames(sera))])
colnames(elisa) = gsub("ELISA_HA_", "", colnames(elisa))


m_elisa = melt(elisa, id.vars = c("Sample_ID", "Age"))
colnames(m_elisa)[grep("variable", colnames(m_elisa))] = c("Test")
colnames(m_elisa)[grep("value", colnames(m_elisa))] = c("Titer")
m_elisa$Titer = log2(m_elisa$Titer)

breaks = seq(0, 10, 2)
labels = (2^(breaks))

#overview
ggplot(m_elisa, aes(x=Age, y=Titer)) +
  geom_jitter(aes(col=Test), width=0, height=0) +
  geom_smooth(aes(col=Test)) +
  theme(text = element_text(size=18)) +
  scale_color_discrete(name="") +
  scale_y_continuous(breaks=breaks, labels = labels) +
  ggtitle(assay)
ggsave(paste0(fig_dir, "overview.png"), width=6, height=4)


#assign age group
elisa = m_elisa

elisa = transform( elisa, Age_group = cut( Age, age_group_1 ) )
elisa$Age_group = gsub("\\(", "", elisa$Age_group)
elisa$Age_group = gsub("\\]", "", elisa$Age_group)
cc = strsplit(elisa$Age_group, ",")
part1 = unlist(cc)[2*(1:nrow(elisa))-1]
part2 = unlist(cc)[2*(1:nrow(elisa)) ]
part1 = as.numeric(part1)+1
elisa$Age_group = paste0(part1, "-", part2)


if("3-10" %in% elisa$Age_group){
  elisa$Age_group = factor(elisa$Age_group, 
                          levels = c("1-2", "3-10", "11-20", "21-30", "31-40", 
                                     "41-50", "51-60", "61-70", "71-90"))
}
####################################################################
# Resample
#######################################################################

elisa_w = reshape(elisa,
                  timevar = "Test",
                  idvar = c("Sample_ID", "Age", "Age_group"),
                  direction = "wide")


#make dataframe with new age group assignment

elisa2 = elisa[elisa$Sex != "",]

elisa2 = transform( elisa2, Age_group = cut( Age, age_group_2) )
elisa2$Age_group = gsub("\\(", "", elisa2$Age_group)
elisa2$Age_group = gsub("\\]", "", elisa2$Age_group)
cc = strsplit(elisa2$Age_group, ",")
part1 = unlist(cc)[2*(1:nrow(elisa2))-1]
part2 = unlist(cc)[2*(1:nrow(elisa2)) ]
part1 = as.numeric(part1)+1
elisa2$Age_group = paste0(part1, "-", part2)

elisa2_w = reshape(elisa2,
                  timevar = "Test",
                  idvar = c("Sample_ID", "Age", "Age_group"),
                  direction = "wide")

#############################################################

data = elisa
data2 = elisa2
data_w = elisa_w
data2_w = elisa2_w

source("analysis_for_assay.R")
