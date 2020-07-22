

bootstrap_test_greater_pval = function(theta_hat, theta_stars, numBoot){
  
  #if (theta_hat == 0){
  #  return (1)
  #}
  
  if(theta_hat < theta_stars[1]){
    return (1)
  }
  for(i in 1:(length(theta_stars)-1)){
    if (theta_hat > theta_stars[i] & theta_hat <= theta_stars[i+1]){
      p = (length(theta_stars)-i)/length(theta_stars)
      return (p)
    }
  }
  if(theta_hat > theta_stars[length(theta_stars)]){
    return (0)
  }
  print (theta_hat)
}

bootstrap_test_less_pval = function(theta_hat, theta_stars, numBoot){
  if (theta_hat >0 | theta_hat == 0){
    return (1)
  }
  
  if(theta_hat > theta_stars[length(theta_stars)]){
    return (1)
  }
  
  for(i in 1:(length(theta_stars)-1)){
    if (theta_hat >= theta_stars[i] & theta_hat < theta_stars[i+1]){
      p = i/length(theta_stars)
      return (p)
    }
  }
  if(theta_hat < theta_stars[length(theta_stars)]){
    return (0)
  }
}

boot_test_every_pair = function(resampled_mean1s, resampled_mean2s, mean1, mean2, mean_by_group){
  #When test is done for every pair of groups, one-sided test("greater") is performed
  
  
  diffs = data.frame(matrix(ncol=0, nrow=nrow(mean_by_group)))
  ci_1s = data.frame(matrix(ncol=0, nrow=nrow(mean_by_group)))
  ci_2s = data.frame(matrix(ncol=0, nrow=nrow(mean_by_group)))
  p_results = data.frame(matrix(ncol=0, nrow=nrow(mean_by_group)))
  
  idx_ci1 = numBoot*0.025
  idx_ci2 = numBoot*0.975
  
  for(i in 1:ncol(resampled_mean1s)){
    mean1s = resampled_mean1s[,i]
    diff_i = c()
    ci_i1 = c()
    ci_i2 = c()
    p_i = c()
    
    for(j in 1:ncol(resampled_mean2s)){
      mean2s = resampled_mean2s[,j]
      
      #theta hat (effect size)
      theta_hat = mean_by_group[i, mean1] - mean_by_group[j, mean2]
      diff_i = c(diff_i, theta_hat)
      
      #theta stars
      theta_stars = mean1s-mean2s
      theta_stars = sort(theta_stars)
      theta_stars = theta_stars - theta_hat
      
      #CI
      CI1 = theta_stars[idx_ci1] + theta_hat
      CI2 = theta_stars[idx_ci2] + theta_hat
      ci_i1 = c(ci_i1, CI1)
      ci_i2 = c(ci_i2, CI2)
      
      #p value
      pval = bootstrap_test_greater_pval(theta_hat, theta_stars, numBoot)
      p_i = c(p_i, pval) 
    }
    diffs = data.frame(diffs, diff_i)
    ci_1s = data.frame(ci_1s, ci_i1)
    ci_2s = data.frame(ci_2s, ci_i2)
    p_results = data.frame(p_results, p_i)
    
  }
  colnames(p_results) = mean_by_group$Age_group
  colnames(diffs) = mean_by_group$Age_group
  p_results$row = mean_by_group$Age_group
  diffs$row = mean_by_group$Age_group
  
  list(p_results, diffs, ci_1s, ci_2s)
}

boot_test_every_pair_less = function(resampled_mean1s, resampled_mean2s, mean1, mean2, mean_by_group){
  #When test is done for every pair of groups, one-sided test("greater") is performed
  
  p_results = data.frame(matrix(ncol=0, nrow=nrow(mean_by_group)))
  
  for(i in 1:ncol(resampled_mean1s)){
    mean1s = resampled_mean1s[,i]
    p_i = c()
    for(j in 1:ncol(resampled_mean2s)){
      mean2s = resampled_mean2s[,j]
      
      theta_hat = mean_by_group[i, mean1] - mean_by_group[j, mean2]
      
      
      theta_stars = mean1s-mean2s
      theta_stars = sort(theta_stars)
      
      theta_stars = theta_stars - theta_hat
      
      pval = bootstrap_test_less_pval(theta_hat, theta_stars, numBoot)
      p_i = c(p_i, pval) 
    }
    p_results = data.frame(p_results, p_i)
  }
  colnames(p_results) = mean_by_group$Age_group
  p_results$row = mean_by_group$Age_group
  
  p_results
}

boot_test_each_group= function(resampled_mean1s, resampled_mean2s, mean1, mean2, mean_by_group){
  #When test is doen for each group, two-side test is performed
  p_results = c()
  diffs = c()
  for(i in 1:ncol(resampled_mean1s)){
    mean1s = resampled_mean1s[,i]
    mean2s = resampled_mean2s[,i]
    
    theta_hat = mean_by_group[i, mean1] - mean_by_group[i, mean2]
    diffs = c(diffs, theta_hat)
    
    theta_stars = mean1s-mean2s
    theta_stars = sort(theta_stars)
    theta_stars = theta_stars - theta_hat
    
    if (theta_hat >= 0 ){
      pval = bootstrap_test_greater_pval(theta_hat, theta_stars, numBoot)
    }else if(theta_hat < 0){
      pval = bootstrap_test_less_pval(theta_hat, theta_stars, numBoot)
    }
    
    p_results = c(p_results, pval)
  }
  
  list(p_results, diffs)
}

resample_AG_getFD_by_virus = function (assay_w, titer_tc1, titer_tc2, numBoot){
  resampled_fds = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))

  for(i in 1:numBoot){
    #resample individuals
    resampled_fd = 
      assay_w %>% 
      group_by(Age_group) %>%
      do(sample_n(., size=length(.$Age_group), replace=T)) %>%
      summarise(fd = mean(eval(titer_tc1) - eval(titer_tc2)) )
    
    resampled_fds = rbind(resampled_fds, t(resampled_fd$fd))
  }
  return ( resampled_fds )
}

resample_AG_getMean_by_virus = function (assay_w, titer_tc1, titer_tc2, numBoot){
  resampled_mean1s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  resampled_mean2s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  
  for(i in 1:numBoot){
    #resample individuals
    resampled_mean = 
      assay_w %>% 
      group_by(Age_group) %>%
      do(sample_n(., size=length(.$Age_group), replace=T)) %>%
      summarise(mean1 = mean(eval(titer_tc1)),
                mean2 = mean(eval(titer_tc2)) )
    
    resampled_mean1s = rbind(resampled_mean1s, t(resampled_mean$mean1))
    resampled_mean2s = rbind(resampled_mean2s, t(resampled_mean$mean2))
  }
  return ( list(resampled_mean1s, resampled_mean2s) )
}


resample_AG_get_condFrac_by_virus = function (assay_w, numBoot){
  resampled_frac1s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  resampled_frac2s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  
  
  for(i in 1:numBoot){
    #resample individuals
    resampled_frac = 
      assay_w %>% 
      group_by(Age_group) %>%
      do(sample_n(., size=length(.$Age_group), replace=T)) %>%
      summarise(frac1 = sum( (Titer.3c2.A <1) & (Titer.3c1>=1) )/sum(Titer.3c2.A <1) ,
                frac2 = sum((Titer.A2 <1) & (Titer.3c1>=1))/sum(Titer.A2 <1) 
      )
    resampled_frac1s = rbind(resampled_frac1s, t(resampled_frac$frac1))
    resampled_frac2s = rbind(resampled_frac2s, t(resampled_frac$frac2))
  }
  return ( list(resampled_frac1s, resampled_frac2s) )
}


resample_AG_get_Frac_by_virus = function (assay_w, numBoot){
  resampled_frac1s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  resampled_frac2s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  
  
  for(i in 1:numBoot){
    #resample individuals
    resampled_frac = 
      assay_w %>% 
      group_by(Age_group) %>%
      do(sample_n(., size=length(.$Age_group), replace=T)) %>%
      summarise(frac1 = sum( (Titer.3c2.A <1) & (Titer.3c1>=1) )/length(Age_group) ,
                frac2 = sum((Titer.A2 <1) & (Titer.3c1>=1)) /length(Age_group)
      )
    resampled_frac1s = rbind(resampled_frac1s, t(resampled_frac$frac1))
    resampled_frac2s = rbind(resampled_frac2s, t(resampled_frac$frac2))
  }
  return ( list(resampled_frac1s, resampled_frac2s) )
}



resample_Gender_within_AG_getMean_by_gender = function (assay_w, titer_tc, numBoot){
  resampled_mean1s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  resampled_mean2s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  
  for(i in 1:numBoot){
    resampled_mean = 
      assay_w %>% 
      group_by(.dots=c("Age_group", "Sex")) %>%
      do(sample_n(., size=length(.$Sex), replace=T)) %>%
      summarise(mean = mean(eval(titer_tc)))
    
    resampled_mean = reshape(data.frame(resampled_mean), 
                             timevar = "Sex",
                             idvar = c("Age_group"),
                             direction = "wide")
    
    resampled_mean1s = rbind(resampled_mean1s, t(resampled_mean$mean.FEMALE))
    resampled_mean2s = rbind(resampled_mean2s, t(resampled_mean$mean.MALE))
  }
  return ( list(resampled_mean1s, resampled_mean2s) )
}

plot_titers_by_age_and_significance = function(titer, stars, start, end, x_by, y_pos,
                                               j_w=0.2, j_h=0.2, facet_age = 1, col_by =0, is_frnt = 1){

  if(col_by == 0){
    plot =  ggplot(titer, aes(x=eval(x_by), y=Titer)) +
      geom_jitter(width=j_w, height=j_h, alpha=0.4) +
      #geom_boxplot(alpha=0) +
      stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean,
                   geom = "crossbar", width=0.7, color="black") +
      stat_summary(fun.data=mean_cl_boot, geom="errorbar",
                   size=1, width=0.5, color = "black") +
      scale_y_continuous(breaks=breaks, labels=labels) +
      theme(text= element_text(size=15)) +
      xlab("")
  }else{
    plot =  ggplot(titer, aes(x=eval(x_by), y=Titer)) +
      geom_jitter(aes(col=eval(col_by)), width=j_w, height=j_h, alpha = 0.6) +
      #geom_boxplot(aes(col=eval(col_by)), alpha=0) +
      stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean,
                   geom = "crossbar", width=0.7, aes(col=eval(col_by))) +
      stat_summary(fun.data=mean_cl_boot, geom="errorbar",
                   size=1, width=0.5, aes(col=eval(col_by))) +
      scale_y_continuous(breaks=breaks, labels=labels) +
      theme(text= element_text(size=15)) +
      xlab("")
  }
  
  if (facet_age == 1){
    plot = plot + facet_wrap(~Age_group, nrow=1, strip.position="bottom") +
      xlab("Age group") +
      theme(panel.spacing=unit(0.1, "lines"),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      scale_color_discrete(name="")
  }
    
  if(length(grep("\\*", stars)) != 0){
    
    if (facet_age == 1) {
      annotation_df = data.frame(Age_group=sort(unique(titer$Age_group))[grep("\\*", stars)], 
                                 start=start, 
                                 end=end,
                                 y=y_pos,
                                 label=stars[grep("\\*", stars)])  
    }else{
      annotation_df = data.frame(
                                 start=start, 
                                 end=end,
                                 y=y_pos,
                                 label=stars[grep("\\*", stars)])  
    }
    
    plot = plot +
      geom_signif(data=annotation_df,
                  aes(xmin=start, xmax=end, annotations=label, y_position=y),
                  manual=TRUE, size=1.2, textsize=8, vjust = 0.5) 
  }

  plot
}

get_ci_from_bootstrap = function(resampled){
  ci = data.frame(lower=numeric(), upper=numeric())
  for(c in 1:ncol(resampled)){
    rs = resampled[,c]
    rs = sort(rs)
    nr = length(rs)
    lower = rs[0.025*nr]
    upper = rs[0.975*nr]
    ci = rbind(ci, c(lower, upper))
  }
  colnames(ci) = c("Lower", "Upper")
  ci
}

resample_AG_getMean_by_ratio = function (assay_w, r1, r2, r3, r4, r5, r6, numBoot){
  resampled_mean1s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  resampled_mean2s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  resampled_mean3s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  resampled_mean4s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  resampled_mean5s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  resampled_mean6s = data.frame(matrix(ncol=length(unique(assay_w$Age_group)), nrow=0))
  
  for(i in 1:numBoot){
    #resample individuals
    resampled_mean = 
      assay_w %>% 
      group_by(Age_group) %>%
      do(sample_n(., size=length(.$Age_group), replace=T)) %>%
      summarise(mean1 = mean(eval(r1)),
                mean2 = mean(eval(r2)),
                mean3 = mean(eval(r3)),
                mean4 = mean(eval(r4)),
                mean5 = mean(eval(r5)),
                mean6 = mean(eval(r6)))
    
    resampled_mean1s = rbind(resampled_mean1s, t(resampled_mean$mean1))
    resampled_mean2s = rbind(resampled_mean2s, t(resampled_mean$mean2))
    resampled_mean3s = rbind(resampled_mean3s, t(resampled_mean$mean3))
    resampled_mean4s = rbind(resampled_mean4s, t(resampled_mean$mean4))
    resampled_mean5s = rbind(resampled_mean5s, t(resampled_mean$mean5))
    resampled_mean6s = rbind(resampled_mean6s, t(resampled_mean$mean6))
  }
  return ( list(resampled_mean1s, resampled_mean2s, resampled_mean3s, resampled_mean4s, resampled_mean5s, resampled_mean6s) )
}