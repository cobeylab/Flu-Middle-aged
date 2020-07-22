library(reshape2)
library(ggplot2)
library(dplyr)
library(Hmisc)
library(ggsignif)
library(lmtest)
library(multiwayvcov)

source("utils.R")
source("boot_meanTiter.R")

assay = "NA ELLA"


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
  fig_dir = "../extended_figure/ella/"
  tb_dir = "../extended_table/ella/"
  age_group_1 = c(0,2,10,20,30,40,50,60,70,90)
  age_group_2 = c(0, 17, 40, 60, 90)
}
#########################################################################

sera = read.table("../../Data/titer/titer.csv", sep=",", head=T, stringsAsFactors = F)
ella = data.frame(sera[, c("Sample_ID", "Age")], sera[, grep("ELLA", colnames(sera))])
ella = na.omit(ella)
colnames(ella) = gsub("ELLA_", "", colnames(ella))
m_ella = melt(ella, id.vars = c("Sample_ID", "Age"))
colnames(m_ella)[grep("variable", colnames(m_ella))] = c("Test")
colnames(m_ella)[grep("value", colnames(m_ella))] = c("Titer")
m_ella$Titer = log2(m_ella$Titer/10)

breaks = c(seq(0, 10, 2))
labels = (2^(breaks))*10


ggplot(m_ella, aes(x=Age, y=Titer)) +
  geom_jitter(aes(col=Test), width=width_o, height=height_o) +
  geom_smooth(aes(col=Test)) +
  theme(text = element_text(size=18)) +
  scale_color_discrete(name="") +
  scale_y_continuous(breaks=breaks, labels = labels) +
  ggtitle(assay)
#ggsave(paste0(fig_dir, "overview.png"), width=6, height=4)


ella = m_ella

ella = transform (ella, Age_group = cut(Age, age_group_1))
ella$Age_group = gsub("\\(", "", ella$Age_group)
ella$Age_group = gsub("\\]", "", ella$Age_group)
cc = strsplit(ella$Age_group, ",")
part1 = unlist(cc)[2*(1:nrow(ella))-1]
part2 = unlist(cc)[2*(1:nrow(ella)) ]
part1 = as.numeric(part1)+1
ella$Age_group = paste0(part1, "-", part2)

if("3-10" %in% ella$Age_group){
  ella$Age_group = factor(ella$Age_group, 
                          levels = c("1-2", "3-10", "11-20", "21-30", "31-40", 
                                     "41-50", "51-60", "61-70", "71-90"))
}
####################################################################
# Resample
#######################################################################

ella_w = reshape(ella,
                  timevar = "Test",
                  idvar = c("Sample_ID", "Age", "Age_group"),
                  direction = "wide")


#make dataframe with new age group assignment

ella2 = ella[ella$Sex != "",]

ella2 = transform( ella2, Age_group = cut( Age, age_group_2) )
ella2$Age_group = gsub("\\(", "", ella2$Age_group)
ella2$Age_group = gsub("\\]", "", ella2$Age_group)
cc = strsplit(ella2$Age_group, ",")
part1 = unlist(cc)[2*(1:nrow(ella2))-1]
part2 = unlist(cc)[2*(1:nrow(ella2)) ]
part1 = as.numeric(part1)+1
ella2$Age_group = paste0(part1, "-", part2)

ella2_w = reshape(ella2,
                   timevar = "Test",
                   idvar = c("Sample_ID", "Age", "Age_group"),
                   direction = "wide")

######################################################################3

data = ella
data2 = ella2
data_w = ella_w
data2_w = ella2_w

source("analysis_for_assay.R")
