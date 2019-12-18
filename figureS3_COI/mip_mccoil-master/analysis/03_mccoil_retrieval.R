dir <- "/home/oj/net/Malaria/OJ/McCOILR_Results/MIPS"
sub_dirs <- list.dirs(dir)[-1]
files <- lapply(sub_dirs, list.files, full.names = TRUE)
names(files) <- gsub(".*/MIPS/","",sub_dirs)

# get categorical summaries
get_cat_sens <- function(fils){

  sums <- lapply(fils[grep("summary",fils)],data.table::fread)
  sums <- lapply(sums, function(x){
    x$group <- gsub("(.*)/(.*)","\\2",gsub("(.*)(_rep_[0-9]\\..*)","\\1",x$file))
    return(x)
  })
  sums <- data.table::rbindlist(sums)

  sums$gt <- gsub("gtcall_([0-9]\\.[0-9]{1,2})_.*","\\1",sums$group)
  sums$region <- gsub("gtcall_([0-9]\\.[0-9]{1,2})_(.*)","\\2",sums$group)
  sums$region_denom <- gsub("(.*)\\.([0-9])","\\2",sums$region)

  group <- dplyr::group_by(sums[sums$CorP=="C",],name,region,region_denom,gt) %>%
    summarise(mean = mean(mean),
              median = median(median),
              sd = mean(sd),
              quantile0.025 = median(quantile0.025),
              quantile0.975 = median(quantile0.975),
              n = length(region))

  return(group)

}

summary <- get_cat_sens(files$cat_output)
saveRDS(summary,"analysis/data/derived/cat_rmcl.rds")

drc_summary <- get_cat_sens(grep("DRC_",files$cat_output,value=TRUE))
saveRDS(drc_summary,"analysis/data/derived/drc_cat_rmcl.rds")

