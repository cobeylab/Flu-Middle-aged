get_stat_table = function (p_and_diff) {
  
  dt = data.frame(matrix("0", nrow=9, ncol=9), stringsAsFactors = F)
  
  for (i in 1:nrow(p_and_diff[[1]])) {
    onerow = c()
    for (j in 1: (ncol(p_and_diff[[2]]) - 1) ) {
      
      diff = format( round(p_and_diff[[2]][i,j], 1), nsmall=1, scientific=F )
      
      ci1 = format( round(p_and_diff[[3]][i,j], 1), nsmall=1, scientific=F )
      ci2 = format( round(p_and_diff[[4]][i,j], 1), nsmall=1, scientific=F )
      
      ci = paste0("[", ci1, ",", ci2, "]")
      
      p = p_and_diff[[1]][i,j]
      if (p < 0.0007) {
        p = "< 0.0007"
      }else {
        p = format( round(p, 3), nsmall=3, scientific=F )
      }
      p = paste("(", p, ")")
      
      dt[i,j] = paste(diff, ci, p )
      if (i == j){
        dt[i,j] = ""
      }
    }
    
  }
  return (dt)
  
}

get_stars = function(pvals, sbtr = 0){
  stars = c()
  for(i in 1:length(pvals)) {
    if(pvals[i] < 0.05/(length(pvals)-sbtr) ) {
      stars = c(stars, "*")
    }else{
      stars = c(stars, "NS")
    }
  }
  stars 
}


get_p_from_bootstrap_fraction_greater = function(x1, x2, n1, n2, numBoot){
  theta_hat = x1/n1 - x2/n2
  if (theta_hat < 0 | ((x1==0) & (x2==0))){
    return (1)
  }
  
  theta_stars = rbinom(numBoot, size=n1, prob=(x1/n1))/n1 -
    rbinom(numBoot, size=n2, prob=(x2/n2))/n2
  theta_stars = sort(theta_stars - theta_hat)
  
  if(theta_hat < theta_stars[1]){
    return (1)
  }
  for(i in 1:(length(theta_stars)-1)){
    if (theta_hat >= theta_stars[i] & theta_hat < theta_stars[i+1]){
      p = (numBoot-i)/numBoot
      return (p)
    }
  }
  if(theta_hat >= theta_stars[numBoot]){
    return (0)
  }
}


bootstrap_test_result_for_every_pair = function(frac_df, num, leng, numBoot){
  pt_results = data.frame(matrix(ncol=0, nrow=nrow(frac_df)))
  for(i in 1:nrow(frac_df)){
    pt_i = c()
    x1 = frac_df[i, num] 
    n1 = frac_df[i, leng]
    for(j in 1:nrow(frac_df)){
      x2 = frac_df[j, num]
      n2 = frac_df[j, leng]
      p = get_p_from_bootstrap_fraction_greater(x1, x2, n1, n2, numBoot)
      pt_i = c(pt_i, p) 
      print (p)
    }
    pt_results = data.frame(pt_results, pt_i)

  } 
  colnames(pt_results) = frac_df$Age_group
  pt_results$row = frac_df$Age_group
  pt_results
}

plot_test_by_ageGroup = function(data_sub, fig_dir, anal, tv, p){
  
  #max titer
  mxtiter = max(data_sub$Titer) + 2
  
  #extract only significant results for geom_signif 
  p_m = melt(p, id.vars = c("row"))
  colnames(p_m) = c("xmax","xmin", "p")
  stars = get_stars(p_m$p, sbtr =nrow(p))
  p_m$stars = stars
  p_m = p_m[grep("\\*", stars),]
  
  print (p_m)
  if(nrow(p_m)==0){
    p1 = ggplot(data_sub, aes(x=Age_group, y=Titer)) +
      geom_jitter(width=0.2, height=0.2, alpha=0.4) +
      #geom_boxplot(alpha=0) +
      stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean,
                   geom = "crossbar", width=0.7, color="black") +
      stat_summary(fun.data=mean_cl_boot, geom="errorbar",
                   size=1, width=0.5, color = "black") +
      #ggtitle( paste(tv, assay) ) +
      theme_classic() +
      theme(text = element_text(size=13)) +
      scale_y_continuous(breaks=breaks, labels = labels) + 
      xlab("Age group")
  }else {
    #gap = (mxtiter/2)/(nrow(p_m))
    gap = 0.5
    
    p1 = ggplot(data_sub, aes(x=Age_group, y=Titer)) +
      geom_jitter(width=0.2, height=0.2, alpha=0.2) +
      #geom_boxplot(alpha=0) +
      stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean,
                   geom = "crossbar", width=0.7, color="black") +
      stat_summary(fun.data=mean_cl_boot, geom="errorbar",
                   size=1, width=0.5, color = "black") +
      #ggtitle( paste(tv, assay) ) +
      theme_classic() +
      coord_fixed(ratio = 2/4, ylim=c(0, 12)) +
      theme(text = element_text(size=13)) +
      scale_y_continuous(breaks=breaks, labels = labels) + 
      
      xlab("Age group") +
      geom_signif(annotation = p_m$stars, textsize= 5,
                  xmin = p_m$xmin, xmax= p_m$xmax, vjust = 0.6,
                  y_pos = seq(mxtiter- (nrow(p_m)-1)*gap + 1, mxtiter + 1.1, gap) ) 
      
    
    
  }
  
  
  p1
}


plot_test_by_ageGroup_fraction = function(data_sub, fig_dir, anal, tv, p){
  
  #max titer
  mxtiter = max(data_sub$Fraction) 
  
  #extract only significant results for geom_signif 
  p_m = melt(p, id.vars = c("row"))
  colnames(p_m) = c("xmax","xmin", "p")
  stars = get_stars(p_m$p, sbtr =nrow(p))
  p_m$stars = stars
  p_m = p_m[grep("\\*", stars),]
  
  p1 = ggplot(data_sub, aes(x=Age_group, y=Fraction)) +
    geom_point(size=2) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper), size=1.2) +
    ylab("Fraction FRNT titer >=20") +
    xlab("Age group") +
    theme(text = element_text(size=13)) +
    geom_signif(annotation = p_m$stars,
                xmin = p_m$xmin, xmax= p_m$xmax, 
                y_pos = seq(1.1, 1.1+ (nrow(p_m)-1)*0.05, 0.05),
                vjust = 0.5)
  
  p1
}

plot_significance_every_pair = function(test_results, plot_title, diffs){
  
  #melt p values
  m_test_results = melt(test_results, id.vars = c("row"))
  colnames(m_test_results)[c(2,3)] = c("column", "p")
  
  #color == 1 if significant after multiple test correction
  alpha = 0.05/(nrow(test_results)*(nrow(test_results)-1))
  m_test_results$col = ifelse(m_test_results$p < alpha , 1, 0)
  m_test_results$col = as.factor(m_test_results$col)
  
  #difference between means or proportions
  m_diffs = melt(diffs, id.vars=c("row"))
  m_test_results$diff = m_diffs$value

  #labels to report p values and difference  

  m_test_results$p = ifelse( m_test_results$p < 0.0001, "<0.0001", 
                            formatC( round(m_test_results$p, 4), format="f") )
  m_test_results$diff = round(m_test_results$diff, 3)
  m_test_results$label = paste0(m_test_results$p, "\n", "(", m_test_results$diff, ")") 
  #m_test_results$label = ifelse( m_test_results$row == m_test_results$column, "", m_test_results$label )
  m_test_results$label = ifelse( m_test_results$col == 0, "", m_test_results$label )
  
                                  
  p_sig = ggplot(m_test_results, aes(x=column, y=row)) +
    geom_tile(aes(fill=col), col="white") +
    geom_text(aes(label=label), size = 2.4) +
    scale_fill_manual(values = c("white", "red"), labels = c("", "Significant"),
                      name= "") +
    ylab("") +xlab("") +
    ggtitle(plot_title) +
    theme(text=element_text(size=13)) +
    theme(legend.position="bottom")
}

plot_numbers_and_significance = function(titer_bin, stars, start, end, titer){
  if (length(grep("\\*", stars)) == 0){
    plot = ggplot(titer_bin[titer_bin$Test == start | titer_bin$Test == end,], aes(x=Test)) +
      geom_bar(aes(fill=YN)) +
      facet_wrap(~Age_group) +
      scale_fill_discrete(labels=c("Titer <20", "Titer >=20"), name = "")
  }else{
    annotation_df = data.frame(Age_group=unique(titer_bin$Age_group)[grep("\\*", stars)], 
                               start=start, 
                               end=end,
                               y=120,
                               label=stars[grep("\\*", stars)])    
    plot = ggplot(titer_bin[titer_bin$Test == start | titer_bin$Test == end,], aes(x=Test)) +
      geom_bar(aes(fill=YN)) +
      facet_wrap(~Age_group) +
      geom_signif(data=annotation_df,
                  aes(xmin=start, xmax=end, annotations=label, y_position=y),
                  manual=TRUE) +
      scale_fill_discrete(labels=c("Titer <20", "Titer >=20"), name = "")
  }
  plot
}

 


plot_fraction_and_significanc = function(m_frac_no_titer, stars, start, end, exclude, clade1, clade2){
  annotation_df = data.frame(Age_group=unique(m_frac_no_titer$Age_group)[grep("\\*", stars)], 
                             start=start, 
                             end=end,
                             y=0.9,
                             label=stars[grep("\\*", stars)])
  plot = ggplot(m_frac_no_titer[m_frac_no_titer$Clade != exclude,], aes(x=Clade, y=Fraction)) +
    geom_point()+
    geom_errorbar(aes(ymin=Lower, ymax=Upper)) +
    facet_wrap(~Age_group) +
    scale_x_discrete(labels=c(clade1, clade2)) +
    geom_signif(data=annotation_df,
                aes(xmin=start, xmax=end, annotations=label, y_position=y),
                manual=TRUE) +
    ylim(c(0,1))
  plot
}

show_bootstrap_every_pair = function(stat_df, m_stat_df, testClade, 
                                    tcFraction, tcNumNoTiter, leng,
                                    pTitle, testPTitle, numBoot, dir, anal){
  plot1 = ggplot(m_stat_df[m_stat_df$Clade == tcFraction,], aes(x=Age_group, y=Fraction)) +
    geom_point(size=2) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper), size=1) +
    ggtitle(pTitle) +
    theme(text = element_text(size=15)) +
    xlab("Age group")
  
  pt_results = bootstrap_test_result_for_every_pair(stat_df, tcNumNoTiter, leng, numBoot)
  print (pt_results)
  plot2 = plot_significance(pt_results, testPTitle)
  
  print (plot1)
  ggsave(paste0(dir, "frac_", testClade, "_", anal, ".png"), width=5, height=4)
  print (plot2)
  ggsave(paste0(dir, "significance_frac_", testClade, "_", anal, ".png"), width=5, height=4)
}

get_m_stat_df = function(stat_df, nums, fracs, lengs){

  cis = c()
  for (i in 1:length(nums)){
    ci = binconf(stat_df[,nums[i]], stat_df[,lengs[i]])
    cis = rbind(cis, ci)
  }

  
  m_stat_df = melt(stat_df[,c("Age_group", fracs, lengs)], id.vars = c("Age_group", lengs),
                   value.name = "Fraction", variable.name = "Clade")
  m_stat_df = data.frame(m_stat_df, cis)
  
  m_stat_df
}