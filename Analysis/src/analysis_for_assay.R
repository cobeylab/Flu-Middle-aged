
#############################################
# Resample and calculate mean by age group
assay_w = data_w

titer_tc1 = "Titer.3c2.A"
titer_tc2 = "Titer.A2"
titer_tc1 = as.symbol(titer_tc1)
titer_tc2 = as.symbol(titer_tc2)

mean_by_AG = 
  assay_w %>% 
  group_by(Age_group) %>%
  summarise(mean1 = mean(eval(titer_tc1)),
            mean2 = mean(eval(titer_tc2)) )
mean_by_AG = data.frame(mean_by_AG)

mean_by_group = mean_by_AG
numBoot = 20000

filename = paste0(tb_dir, "resampled_mean12.rds")
if(file.exists(filename)) {
  resampled_mean12 = readRDS(filename)
}else {
  resampled_mean12 = resample_AG_getMean_by_virus(assay_w, titer_tc1, titer_tc2, numBoot)
  saveRDS(resampled_mean12, filename)
}


#############################################################################

###########################################################################
# Titer by virus
###########################################################################

mxtiter = max(data$Titer)
anal = "titer_by_virus"

#Regression
fit = lm(Titer ~ Test, data=data)
cl = cluster.vcov(fit, data$Sample_ID)
coef = coeftest(fit, cl)
pval = coef[2, 4]
star = get_stars(pval)
write.csv(coef, paste0(tb_dir, anal, ".csv"))

#Bootstrap test
mean12 = colMeans(mean_by_AG[,c("mean1", "mean2")])  
theta_hat = mean12[1] - mean12[2]

theta_stars = rowMeans(resampled_mean12[[1]]) - rowMeans(resampled_mean12[[2]])
theta_stars = theta_stars - theta_hat
theta_stars = sort(theta_stars)

if (theta_hat >= 0 ){
  pval_boot = bootstrap_test_greater_pval(theta_hat, theta_stars, numBoot)
}else if(theta_hat < 0){
  pval_boot = bootstrap_test_less_pval(theta_hat, theta_stars, numBoot)
}
star_boot = get_stars(pval_boot)

#check if two tests agree
stopifnot(star_boot == star)

#write mean difference and p value
pval_boot = ifelse(pval_boot < 0.001, "<0.001", pval_boot)
write.csv( c(theta_hat, pval_boot), paste0(tb_dir, anal, "_bootstrap.csv") )


#Figure

mxtiter = max(data$Titer)
plot = plot_titers_by_age_and_significance(data, star_boot, start = "3c2.A", end = "A2", 
                                           x_by = as.symbol("Test"), y_pos = mxtiter-1, 
                                           j_w=width_v, j_h=height_v, facet_age = 0)
plot = plot + ggtitle(assay) 
plot
#ggsave(paste0(fig_dir, anal, ".png"), width=5, height=4)  

##########################################################################
# Between each age group
##########################################################################
anal = "titer_by_age_"

#3c2.A

tv = "3c2.A"
data_sub = data[data$Test == tv,]

p_and_diff = boot_test_every_pair(resampled_mean12[[1]], resampled_mean12[[1]], "mean1", "mean1", mean_by_group)
dt = get_stat_table(p_and_diff)
write.csv(dt, paste0(tb_dir, tv, "_dt.csv"))
plot = plot_test_by_ageGroup(data_sub, fig_dir, anal, tv, p_and_diff[[1]])
plot
#ggsave(paste0(fig_dir, anal, tv, ".png"), width=6, height=4)
#saveRDS(plot, paste0(fig_dir, anal, tv, ".rds"))

ptitle = paste0("Significantly higher column compared to row",
                "\n(", tv, ")")
plot = plot_significance_every_pair(p_and_diff[[1]], ptitle, p_and_diff[[2]])
plot
#ggsave(paste0(fig_dir, "significance_", anal, tv, ".png"), width=5, height=4)


#A2

tv = "A2"
data_sub = data[data$Test == tv,]

p_and_diff = boot_test_every_pair(resampled_mean12[[2]], resampled_mean12[[2]], "mean2", "mean2", mean_by_group)
dt = get_stat_table(p_and_diff)
write.csv(dt, paste0(tb_dir, tv, "_dt.csv"))
plot = plot_test_by_ageGroup(data_sub, fig_dir, anal, tv, p_and_diff[[1]])
plot
#ggsave(paste0(fig_dir, anal, tv, ".png"), width=6, height=4)
#saveRDS(plot, paste0(fig_dir, anal, tv, ".rds"))

ptitle = paste0("Significantly higher column compared to row",
                "\n(", tv, ")")
plot = plot_significance_every_pair(p_and_diff[[1]], ptitle, p_and_diff[[2]])
plot
#ggsave(paste0(fig_dir, "significance_", anal, tv, ".png"), width=5, height=4)

#################################################################################
# By virus for each age group
####################################################################################
# anal = "titer_by_virus_for_each_AG"
# 
# #Regression
# 
# AG = sort(unique(data$Age_group))
# pvals = c()
# for(i in 1:length(AG)) {
#   df = data[data$Age_group == AG[i],]
#   fit = lm(Titer ~ Test, df)
#   vcov = cluster.vcov(fit, df$Sample_ID)
#   coef = coeftest(fit, vcov)
#   pvals = c(pvals, coef[2,4])
# }
# stars = get_stars(pvals)
# 
# #Bootstrap
# p_and_diff = boot_test_each_group(resampled_mean12[[1]], resampled_mean12[[2]], "mean1", "mean2", mean_by_group)
# stars_boot = get_stars(p_and_diff[[1]])
# 
# x_by = as.symbol("Test")
# plot = plot_titers_by_age_and_significance(data, stars_boot, start = "3c2.A", end = "A2", 
#                                            x_by = x_by, y_pos = mxtiter-1, j_w = width_va, j_h = height_va,
#                                            facet_age = 1, col_by = as.symbol("Test"))
# plot = plot + ggtitle(assay)
# plot
# ggsave(paste0(fig_dir, anal, ".png"), width=9, height=5)
# 
# 
# #mean difference and p value
# #sig_idx = grep("\\*", stars_boot)
# sig_idx = grep("*", stars_boot)
# 
# ag = sort(unique(data$Age_group))[sig_idx]
# mean_diffs = signif(p_and_diff[[2]][sig_idx], 3)
# pvals = ifelse(p_and_diff[[1]][sig_idx] < 0.001, "<0.001",  
#                round(p_and_diff[[1]][sig_idx], 3) )
# df_each_ag = data.frame(ag, mean_diffs, pvals)
# colnames(df_each_ag) = c("Age group", "Mean difference", "p-value")
# write.csv(df_each_ag, paste0(tb_dir, anal, ".csv") )
# 
# 
