# trends in COI and missingness
library(tidyverse)
library(rdhs)

figure_dir <- "analysis/manuscript_figures/"

# load our samples
samples <- readRDS("analysis/data/derived/monoclonal_samples.rds")
samples$symptomatic <- TRUE
samples$symptomatic[samples$Country %in% c("DRC","Zambia")] <- FALSE

source("analysis/facet_zoom.R")

# COI disstributions first
coi_dist <- ggplot(samples, aes(x = Country, y = coi)) + geom_violin(trim=FALSE) +
  stat_summary(fun.data = "mean_cl_boot", fun.args = list(B=10000), geom = "pointrange", color = "red") +
  ylab("COI") + scale_y_continuous(breaks = 0:8) + theme(panel.grid.major.y = element_line(color="grey"))

cowplot::save_plot(paste0(figure_dir,"/sup_fig3_v1.png"),
                   coi_dist, base_height = 6, base_width = 9, dpi=300)
cowplot::save_plot(paste0(figure_dir,"/sup_fig3_v1.svg"),
                   coi_dist, base_height = 6, base_width = 9, dpi=300)

cowplot::save_plot(paste0(figure_dir,"/sup_fig3_v2.png"),
                   coi_dist + facet_zoom2(ylim=c(1.25,3),zoom.size = 0.5),
                   base_height = 6, base_width = 12, dpi=300)
cowplot::save_plot(paste0(figure_dir,"/sup_fig3_v2.svg"),
                   coi_dist + facet_zoom2(ylim=c(1.25,3),zoom.size = 0.5),
                   base_height = 6, base_width = 12, dpi=300)

# and values of boostrap
boots <- group_by(samples,Country) %>% summarise(b = list(Hmisc::smean.cl.boot(coi,B=10000)),
                                                 mean = b[[1]][1],
                                                 low = b[[1]][2],
                                                 high = b[[1]][3]) %>%
  select(Country,mean,low,high)
write.table(cbind(boots[,1],round(boots[,2:4],3)),paste0(figure_dir,"/country_bootstrap_ci.csv"),)

boot_symp <- group_by(samples,symptomatic) %>% summarise(b = list(Hmisc::smean.cl.boot(coi,B=10000)),
                                                         mean = b[[1]][1],
                                                         low = b[[1]][2],
                                                         high = b[[1]][3]) %>%
  select(symptomatic,mean,low,high)


# alternative way of looking at significant differences between counties by coi ranks
coi_df <- data.frame("coi"=samples$coi, "country"=samples$Country)
coi_df_wide <- mltools::one_hot(data.table::as.data.table(coi_df))
mod <- lm(rank(coi) ~ 0 + country_DRC + country_Ghana + country_Tanzania + country_Uganda + country_Zambia, data = coi_df_wide)
kr <- PMCMR::posthoc.kruskal.nemenyi.test(samples$coi, factor(samples$Country), dist = "Tukey")


## COMPARE TO DHS COVARIATES OF PREAVALENCE
# ------------------------------------------------------------------------------
# load our samples and subset to DRC

samples <- samples[samples$Country=="DRC",]

match_clean <- function(a,b){
  a <- gsub("[[:punct:][:space:]]","",tolower(stringi::stri_trans_general(a, "latin-ascii")))
  b <- gsub("[[:punct:][:space:]]","",tolower(stringi::stri_trans_general(b, "latin-ascii")))
  ret <- match(a,b)
  if(sum(is.na(ret)>0)){
    dists <- stringdist::seq_distmatrix(lapply(a,utf8ToInt),lapply(b,utf8ToInt))
    ret[is.na(ret)] <- apply(dists[which(is.na(ret)),,drop=FALSE],1,which.min)
    print(unique(cbind(a,b[ret])))
  }
  return(ret)
}

# grab dhs relevant surveys and extract relevant information
drcs <- dhs_datasets(countryIds = "CD",fileFormat = "FL", fileType = c("GE","PR"), surveyYear = 2013)
dats <- get_datasets(drcs)
vars <- get_variable_labels(dats$CDPR61FL)
barcs <- search_variable_labels(names(dats),"age|prov|smear|date|weight|residence|Cluster number|cta|anemia|net|malaria measurement|rapid test|blood sample|Case Identification|^Line number$")
extr  <- extract_dhs(barcs,add_geo = TRUE)[[1]]
vars <- get_variable_labels(extr)

dhs_admin_rdt <- rdhs::dhs_data(countryIds = "CD",surveyYear = "2013",indicatorIds = c("ML_PMAL_C_RDT"),breakdown = "subnational")
dhs_admin_micro <- rdhs::dhs_data(countryIds = "CD",surveyYear = "2013",indicatorIds = c("ML_PMAL_C_MSY"),breakdown = "subnational")
extr$rdt_dhs <- dhs_admin_rdt$Value[match_clean(names(attr(extr$hv024,"labels"))[(extr$hv024)],dhs_admin_rdt$CharacteristicLabel)]
extr$micro_dhs <- dhs_admin_micro$Value[match_clean(names(attr(extr$hv024,"labels"))[(extr$hv024)],dhs_admin_rdt$CharacteristicLabel)]


# individual data
samples$age <- extr$hml16a[match(paste0(samples$hhid,samples$hvidx),paste0(extr$hhid,extr$hvidx))]
samples$age[samples$MissingChildren] <- sample(61:72,size = sum(samples$MissingChildren),replace=TRUE)
samples$weight <- extr$hv005[match(samples$hv001,extr$hv001)]
samples$strata <- paste0(samples$ADM1NAME,samples$hv001)
options("survey.lonely.psu"="certainty")

# add seasonality
ses <- sapply(names(attr(extr$hv024,"labels"))[as.numeric(unique(extr$hv024))],
              function(x) {
                if(x=="orientale") x <- "Province orientale"
                magenta:::seasonal_profile(x,"Democratic Republic of the Congo")
              })
seasonal_add <- function(months, region, ses){
  if(is.na(months)){
    return(NA)
  }
  months <- as.numeric(months)
  region <- names(attr(extr$hv024,"labels"))[as.numeric(region)]
  ring <- c(12,1:11)
  mts <- ring[months]
  days <- lubridate::days_in_month(1:12) %>% cumsum
  days <- sapply(1:12,function(x) c(0,days[-12])[x]:days[x])
  return(mean(ses[unlist(days[mts]),region]))
}

extr$seasonal_adm_month <- apply(extr,1,function(x){seasonal_add(x[which(names(extr)=="hc18")],x[which(names(extr)=="hv024")],ses)})

# turn into complex survey design
extr$strata <- paste0(extr$ADM1NAME,extr$hv001)
extr$hv005 <- extr$hv005  /1e6
svy <- extr[extr$hml35 %in% 0:1 & extr$hml32 %in% 0:1,] %>%
  srvyr::as_survey_design(hv001, weights=hv005)
svy_samps <- samples[!is.na(samples$weight),] %>%
  srvyr::as_survey_design(strata=strata, weights=weight)


meta_cluster <- left_join(
  svy %>% group_by(hv001) %>%
    summarise(micro = srvyr::survey_mean(hml35),
              rdt = srvyr::survey_mean(hml32),
              n = srvyr::survey_total(hml35),
              age = srvyr::survey_mean(hml16a)),
  svy$variables %>% group_by(hv001) %>%
    summarise(months = paste(unique(na.omit(hc18)),collapse=" "),
              region = names(attr(extr$hv024,"labels"))[unique(hv024)],
              seasonal = mean(seasonal_adm_month,na.rm=TRUE),
              admin = unique(ADM1NAME)),
  by="hv001") %>%
  left_join(
    svy_samps %>% group_by(hv001) %>%
      summarise(coi = srvyr::survey_mean(coi),
                coi_low = srvyr::survey_mean(coi_low),
                coi_middle = srvyr::survey_mean(coi_middle),
                coi_high = srvyr::survey_mean(coi_high)),
    by="hv001")

meta_admin <- left_join(
  svy %>% group_by(ADM1NAME) %>%
    summarise(micro = srvyr::survey_mean(hml35),
              rdt = srvyr::survey_mean(hml32),
              n = srvyr::survey_total(hml35),
              age = srvyr::survey_mean(hml16a)),
  svy$variables %>% group_by(ADM1NAME) %>%
    summarise(months = paste(unique(na.omit(hc18)),collapse=" "),
              region = names(attr(extr$hv024,"labels"))[unique(hv024)],
              seasonal = mean(seasonal_adm_month,na.rm=TRUE),
              admin = unique(ADM1NAME)),
  by="ADM1NAME") %>%
  left_join(
    svy_samps %>% group_by(ADM1NAME) %>%
      summarise(coi = srvyr::survey_mean(coi),
                coi_low = srvyr::survey_mean(coi_low),
                coi_middle = srvyr::survey_mean(coi_middle),
                coi_high = srvyr::survey_mean(coi_high)),
    by="ADM1NAME")


## ADMIN LEVEL

admin_plot1 <-
  ggplot(meta_admin,
         aes(micro,coi,
             xmin=micro-1.96*micro_se,xmax=micro+1.96*micro_se,
             ymin=coi-1.96*coi_se,ymax=coi+1.96*coi_se)) +
  geom_point(aes(size=n),alpha=0.6) +
  geom_linerange(alpha=0.5) + geom_errorbarh(alpha=0.5) +
  geom_smooth(aes(weight=n),span=3) +
  theme_bw() +
  xlab("Admin Microscopy Prevalence") +
  ylab("COI") + scale_size_continuous(name="Sample Size")

admin_plot2 <-
  ggplot(meta_admin,
         aes(micro,coi)) +
  geom_point(aes(size=n),alpha=0.6) +
  geom_smooth(aes(weight=n),span=3) +
  theme_bw() +
  xlab("Admin Microscopy Prevalence") +
  ylab("COI") + scale_size_continuous(name="Sample Size")


saveRDS(admin_plot,"analysis/plots/admin_coi.rds")

## CLUSTER

# remove clusters we did not have coi estimates for
meta_cluster <- meta_cluster[!is.na(meta_cluster$coi) & meta_cluster$micro>0,]
meta_cluster$prev_bin <- cut(meta_cluster$micro,seq(0,1,0.1))
prev_bin <- ggplot(meta_cluster,aes(x=prev_bin,y=coi)) + geom_boxplot() + theme_bw() +
  xlab("Prevalence Range") + ylab("COI")
saveRDS(prev_bin,"analysis/plots/prev_bin.rds")

# CLuster plot
g <- ggplot(meta_cluster, aes(x=micro,y=coi,weight=weights)) +
  geom_point(aes(size=n),alpha=0.6) +
  geom_smooth(span=10,se=TRUE) +
  theme_bw() + xlab("Cluster Microscopy Prevalence") +
  ylab("COI") + scale_size_continuous(name="Sample Size")
saveRDS(g, "analysis/plots/cluster_coi.rds")

cowplot::save_plot(paste0(figure_dir,"/sup_fig4_v1.png"),
                   cowplot::plot_grid(admin_plot1,g,labels = "auto"),
                   base_height = 6, base_width = 14, dpi=300)
cowplot::save_plot(paste0(figure_dir,"/sup_fig4_v1.svg"),
                   cowplot::plot_grid(admin_plot1,g,labels = "auto"),
                   base_height = 6, base_width = 14, dpi=300)
cowplot::save_plot(paste0(figure_dir,"/sup_fig4_v2.png"),
                   cowplot::plot_grid(admin_plot2, g,labels = "auto"),
                   base_height = 6, base_width = 14, dpi=300)
cowplot::save_plot(paste0(figure_dir,"/sup_fig4_v2.svg"),
                   cowplot::plot_grid(admin_plot2, g,labels = "auto"),
                   base_height = 6, base_width = 14, dpi=300)

## statistical relationship against prevalence
glm_clust <- lme4::glmer(coi ~ micro + (1|region), weights = weights+2, data = meta_cluster, family = gaussian(link="log"),nAGQ = 5)
glm_admin <- lme4::glmer(coi ~ micro + (1|ADM1NAME), weights = weights+2, data = meta_admin, family = gaussian(link="log"),nAGQ = 5)
summary(glm_clust)
summary(glm_admin)

## EXTRA ANALYSIS COMPARING TO IMPERIAL TRANSMISSION MODEL
## -----------------------------------------------------------------------------

# design the relationship based on imperial model predicted EIR and prev relationship
df <- data.frame("EIR"=magenta:::age_brackets(200,100))
df$prev <- vapply(magenta:::age_brackets(200,100), function(x){
  eq <- magenta:::equilibrium_init_create(magenta:::age_brackets(),het_brackets = 5,EIR = x,model_param_list = magenta:::model_param_list_create(num_int = 1),ft = 0.4)
  return(sum(eq$init_D[5:11,,1] + eq$init_T[5:11,,1] + eq$init_A[5:11,,1]*eq$p_det_eq[5:11,])/sum(eq$den[5:11]))
},numeric(1))

# functional form for imperial model relationship
rhs <- function(x,a,b){a*x^b}
m <- nls(EIR ~ I(1+exp(1)^(a +b*prev)), data = df, start = list(a=1,b=1),trace = T)
m_sum <- summary(m)
m2 <- nls(EIR ~ rhs(prev,a,b), data = df, start = list(a=100,b=0.1),trace = T)
m2_sum <- summary(m2)
df$exp <- predict(m,type="response")
df$power <- predict(m2,type="response")
imp_model <- ggplot(df[seq(1,100,length.out = 50),],aes(prev,EIR)) + geom_point() +
  geom_line(mapping=aes(prev,exp,color="green")) +
  geom_line(mapping=aes(prev,power,color="red")) +
  scale_colour_manual(name="Model Fit",
                      values = c('green'='green','red' = 'red'),
                      labels = c("Expoenetial", "Power Model")) +
  theme_bw() + xlab("Microscopy Prevalence") + ylab("EIR")
saveRDS(imp_model,"analysis/plots/imper_model.rds")

# use functional form to predict the relationship and compare it the loess
glm <- glm(coi~I(exp(1)^(m_sum$coefficients[1] + m_sum$coefficients[2]*micro)),weights=n,data=meta_cluster)
glm2 <- glm(coi~I((m2_sum$coefficients[1])*(micro^(m2_sum$coefficients[2]))),weights = n,data=meta_cluster)
glmm <- lme4::lmer(coi~I(exp(1)^(m_sum$coefficients[1] + m_sum$coefficients[2]*micro))+(1|region),data=meta_cluster,weights = n)
glmm2 <- lme4::lmer(coi~I((m2_sum$coefficients[1])*(micro^m2_sum$coefficients[2]))+(1|region),data=meta_cluster,weights = n)
meta_cluster$coi_predict_power <- meta_cluster$coi_predict_exp <- NA
meta_cluster$coi_predict_exp <- predict(glm,type="response")
meta_cluster$coi_predict_power <- predict(glm2,type="response")
meta_cluster$coi_predict_glmm_exp <- predict(glmm,type="response")
meta_cluster$coi_predict_glmm_power <- predict(glmm2,type="response")


g1 <- ggplot(meta_cluster, aes(x=micro,y=coi,weight=n)) +
  geom_point(aes(size=n),alpha=0.6) +
  geom_smooth(span=1,aes(col="blue"),se=FALSE) +
  theme_bw() + xlab("Cluster Microscopy Prevalence") +
  ylab("COI") + scale_size_continuous(name="Cluster Size") +
  #geom_smooth(mapping=aes(y=coi_predict_glmm_power,col="green"),se=F)+
  geom_smooth(mapping=aes(y=coi_predict_glmm_exp,col="red"),se=F) +
  scale_colour_manual(name="Model Fit",
                      values = c('blue'='blue','red' = 'red','green' = 'green'),
                      labels = expression(LOESS,a.prev^b,e^(a+b.prev)))


# compare fits
coi_model_fits <- data.frame(
  "Model"=c("Loess","Power","Exponential","GLMM Power","GLMM Exponential"),
  "RSS"=c(sum(residuals(loess(coi~micro, meta_cluster))^2),
          sum((meta_cluster$coi_predict_power-meta_cluster$coi)^2),
          sum((meta_cluster$coi_predict_exp-meta_cluster$coi)^2),
          sum((meta_cluster$coi_predict_glmm_power-meta_cluster$coi)^2),
          sum((meta_cluster$coi_predict_glmm_exp-meta_cluster$coi)^2)
  )
)

saveRDS(g1,"analysis/plots/coi_prev_plot.rds")
saveRDS(coi_model_fits[c(1,4,5),],"analysis/plots/coi_prev_models.rds")
