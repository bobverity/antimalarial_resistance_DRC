load("analysis/data/derived/mccoil_inputs.RData")


### cluster set up
workdir <- "/home/oj/net/Malaria/OJ/McCOILR_Results"
new_dir <- "MIPS"
dir.create(paste0(workdir,"/",new_dir))
dir.create(paste0(workdir,"/",new_dir,"/cat_output"))
didehpc::didehpc_config_global(workdir=workdir,
                               credentials="/home/oj/.smbcredentials",
                               temp=didehpc::path_mapping("tmp",
                                                          "/home/oj/net/temp",
                                                          "//fi--didef3.dide.ic.ac.uk/tmp",
                                                          "T:"),
                               home=didehpc::path_mapping("OJ",
                                                          "/home/oj/net/Malaria",
                                                          "//fi--didef3.dide.ic.ac.uk/Malaria",
                                                          "M:"),
                               cluster="small")
didehpc::web_login()
root <- file.path(workdir, "contexts")
packages.vector <- c("Rcpp","McCOILR")
context::context_log_start()

## set up context
ctx <- context::context_save(root,
                             packages = packages.vector,
                             package_sources = provisionr::package_sources(
                               github=c("OJWatson/McCOILR")
                             )
)

config <- didehpc::didehpc_config(use_workers = TRUE)
obj <- didehpc::queue_didehpc(ctx, config = config)

# send some workers
id <- obj$rrq$message_send("TIMEOUT_SET",1*24*60*60)

## submission
paramList <- list()
n_chains <- 5
length(paramList) <- n_chains


try_fail_catch <- function(expr, attempts = 3){
  r <- NULL
  attempt <- 1
  while( is.null(r) && attempt <= 3 ) {
    attempt <- attempt + 1
    try(
      r <- eval(expr)
    )
  }

}

# loop over genotype thresholds for calling heterozygote and region allocation
gt_call <- c(0.05,0.1,0.15)
p$samples$KMC1 <- p$samples$Country

# ALL ----

for(g in gt_call){

  gt_data <- assignGTforREALMcCOIL(wsaf,err=g)
  rownames(gt_data) <- strsplit(rownames(gt_data),"-") %>%
    lapply(function(x) paste0(head(x,-2),collapse="-")) %>%
    unlist

  for(k in c(1:5)) {


    if(k >= 2){
      regs <- grep("DRC",unique(p$samples[[paste0("KMC",k)]]),value=TRUE)
    } else {
      regs <- unique(p$samples[[paste0("KMC",k)]])
    }


    for(j in regs){

      for(i in 1:n_chains){

        paramList[[i]] <- (list(data = gt_data[which(p$samples[[paste0("KMC",k)]] == j),],
                                maxCOI=25, totalrun=100000, burnin=1000,
                                M0=15, err_method=3,
                                path="M:/OJ/McCOILR_Results",
                                thin=0.1,
                                output=paste0("MIPS/cat_output/gtcall_",g,"_",j,"_rep_",i,".txt")))

      }

      try_fail_catch(grp2 <- queuer::qlapply(X = paramList,obj = obj, timeout=0,
                                             FUN = function(x){
                                               return(McCOIL_categorical(data = x$data,
                                                                         maxCOI = x$maxCOI,
                                                                         totalrun = x$totalrun,
                                                                         burnin = x$burnin,
                                                                         M0 = x$M0,
                                                                         threshold_ind = 0, threshold_site = 0,
                                                                         thin = x$thin,
                                                                         err_method = x$err_method,
                                                                         path = x$path,
                                                                         output = x$output))
                                             },
                                             name = paste0("MIPS_cat_final_gtcall_",g,"_",j),overwrite = TRUE))


    }
  }
}

