### 1. explore the impact on substructuring in DRC for COI
### 2. Confirm for most this is negligable

# working through sensitivity analysis to produce assured data set
summary <- readRDS("analysis/data/derived/cat_rmcl.rds")
summary$region_denom[nchar(summary$region_denom)>1] <- 1
filter_for_regional_discrepancy <- function(summary, gt=0.05){

  # subset the data for the genotype threshold
  sum_gt <- summary[summary$gt==gt,]

  # spread the regions outwide and decide which have discordance between region
  reg_wide <- spread(sum_gt[,c(1,2,6)],region,median)
  regional_diffs <- which(apply(reg_wide[,-1],1,
                                function(x) length(na.omit(unique(x))))>1)

  # subset to the disagreement
  regional_diff_df <- sum_gt[which(sum_gt$name %in% reg_wide$name[regional_diffs]),]
  regional_diff_df[order(regional_diff_df$name,regional_diff_df$region_denom),]

  # subeset to those that aren't 1
  reg_diff_1 <- regional_diff_df[which(regional_diff_df$name %in%
                                         c(unique(regional_diff_df$name[regional_diff_df$median==1]))),]

  # plot these out
  differences_by_region <- group_by(reg_diff_1,name,region_denom) %>%
    summarise(med = median(median),
              den = as.numeric(region_denom)) %>%
    ggplot(aes(x=den,y=med,group=name)) + geom_line() + facet_wrap(~name,scale="free_y")

  message(paste(length(unique(reg_diff_1$name)),": ",
                paste0(unique(reg_diff_1$name),collapse=" ")))
  invisible(list("names"=unique(reg_diff_1$name),
                 "plot"=differences_by_region,
            "all"=unique(regional_diff_df$name)))

}

concerns <- vector("list", 3); names(concerns) <- c("gt0.05", "gt0.1", "gt0.15")
for(g in c(0.05, 0.1, 0.15)){
  concerns[[paste0("gt",g)]] <- filter_for_regional_discrepancy(summary,g)
}
check <- lapply(concerns,function(x) x$names) %>% unlist
check <- as.character(check[duplicated(check)]) %>% unique
concerns[[1]]$plot <- concerns[[1]]$plot + ggtitle("Hetozygote MAF Threshold = 0.05") +
  xlab("DRC Clusters") + ylab("COI")
concerns[[2]]$plot <- concerns[[2]]$plot + ggtitle("Hetozygote MAF Threshold = 0.1") +
  xlab("DRC Clusters") + ylab("COI")
concerns[[3]]$plot <- concerns[[3]]$plot + ggtitle("Hetozygote MAF Threshold = 0.15") +
  xlab("DRC Clusters") + ylab("COI")

cowplot::plot_grid(cowplot::plot_grid(plotlist = lapply(concerns[1:2],function(x) x$plot),ncol = 1),concerns[[3]]$plot,ncol=2)


discrepancies_across_clusters <- ggplot(
  p$samples[p$samples$ID %in% check,],aes_string(x="long",y="lat")) +
  borders(fill="gray22") + theme_bw() +
  geom_point(col="white") +
  coord_fixed(xlim = c(10,max(p$samples$long,na.rm=TRUE)),
              ylim = c(min(p$samples$lat,na.rm=TRUE),max(p$samples$lat,na.rm=TRUE))) +
  scale_color_discrete(name="Region")


## now we are happy that the number of discrepancies for median = 1 is small across
## substructure, then we can take regions as a whole
summary$country <- gsub("([A-z]*)([0-9]\\.[0-9])","\\1",summary$region)
reg_summarised <- group_by(summary, name,gt,country, region_denom) %>%
  summarise(mean=mean(mean),
            median=median(median),
            sd=mean(sd),
            quantile0.025=min(quantile0.025),
            quantile0.975=max(quantile0.975))


monoclonal_rating <- function(quantile0.975,median,gt) {
  if(sum(quantile0.975>1)==0){
    return(1)
  } else if(sum(median>1)==0){
    return(2)
  } else if(sum(median[gt>0.05]>1)==0){
    return(3)
  } else if(sum(median[gt>0.1]>1)==0){
    return(4)
  } else {
    return(5)
  }
}

monoclonals <- group_by(reg_summarised,name,country) %>%
  summarise(mono_rating = monoclonal_rating(quantile0.975,median,gt),
            coi = median(median),
            coi_low=median[gt==0.15],
            coi_middle=median[gt==0.10],
            coi_high=median[gt==0.05],
            coi_0.025=median(quantile0.025),
            coi_0.975=median(quantile0.975))

monoclonals$regional_diffs <- FALSE
monoclonals$regional_diffs[match(check,monoclonals$name)] <- TRUE

p$samples$coi <- monoclonals$coi
p$samples$coi_0.025 <- monoclonals$coi_0.025
p$samples$coi_0.975 <- monoclonals$coi_0.975
p$samples$monoclonal_rating <- haven::labelled(
  monoclonals$mono_rating,
  c("Definite monoclonal: COI 95%CI == 1 for all het threshold"=1,
    "Likely monoclonal: COI median == 1 for all het threshold"=2,
    "Probably monoclonal: COI median == 1 when het threshold >= 0.1"=3,
    "Weak monoclonal: COI median == 1 when het threshold >= 0.15"=4,
    "Not monoclonal"=5)
)
p$samples$regional_discrepancy <- monoclonals$regional_diffs
p$samples$coi_low <- monoclonals$coi_low
p$samples$coi_middle <- monoclonals$coi_middle
p$samples$coi_high <- monoclonals$coi_high


# save outputs
saveRDS(p$samples, "analysis/data/derived/monoclonal_samples.rds")
saveRDS(check_plots, "analysis/plots/monoclonal_changes.rds")
saveRDS(discrepancies_across_clusters, "analysis/plots/monoclonal_discrpancy_map.rds")

mono_summary <- data.frame("Monoclonality"=names(attr(p$samples$monoclonal_rating,"label")),
                           "Frequency"=as.numeric(p$samples$monoclonal_rating %>% table))
saveRDS(mono_summary,"analysis/plots/monoclonal_summaries.rds")
