# few needed packages
#devtools::install_github("mrc-ide/MIPanalyzer")
#devtools::install_github("andrewparkermorgan/rplasmodium", ref="16b1776")

# load them for ease
library(tidyverse)
library(MIPanalyzer)

# get the data in and the meta data for DRC extracted from the DHS
p <- readRDS("analysis/data/raw/biallelic_distances.rds")


# make the drug sites we need
drug_sites <- rplasmodium::pf_3d7_PutDrugRxSites
names(drug_sites) <- drug_sites[1,]
drug_sites <- drug_sites[-1,]
drug_sites <- rbind(drug_sites, c("Pf3D7_06_v3",1213948,1216005,"-","aat1","PF3D7_0629500"))
drug_sites$chrnum <- paste0("chr",as.numeric(substr(drug_sites$chr,7,8)))
drug_sites$start <- as.numeric(drug_sites$start)
drug_sites$end <- as.numeric(drug_sites$end)
rownames(drug_sites) <- NULL

# add the drugx loci info to the loci meta
dist <- c("10kb"=1e4, "50kb"=5e4, "500kb"=5e5)

# assign for each kb region
for (j in seq_len(length(dist))) {

  p$loci[[paste0("drugrx",names(dist[j]))]] <- NA

  for (i in 1:nrow(drug_sites)) {

    ds <- which( (p$loci$CHROM == drug_sites$chrnum[i]) &
                   ((p$loci$POS > (drug_sites$start[i] - dist[j])) &
                      (p$loci$POS < (drug_sites$end[i] + dist[j]))))

    if (length(ds) > 0) {
    p$loci[[paste0("drugrx",names(dist[j]))]][ds] <- drug_sites$name[i]
    }
  }

}



## split up DRC into n clusters based on geog

# check we are not missing any geo data for DRC
sum(is.na(p$samples$lat[p$samples$Country=="DRC"])) # 0
sum(is.na(p$samples$long[p$samples$Country=="DRC"])) # 0

## assign clustering based on Country and geographical k means for DRC
for (k in c(2:5)){
  p$samples[[paste0("KMC",k)]] <- p$samples$Country
  p$samples[[paste0("KMC",k)]][p$samples[[paste0("KMC",k)]]=="DRC"] <-
    paste0("DRC",
           kmeans(p$samples[p$samples$Country=="DRC",c("lat","long")],
                  centers = k,nstart = 600,iter.max = 10000)$cluster,
           ".",k)
}

# quick plot to check it looks sensible
map <- list()
for(k in seq_len(4)){
col <- paste0("KMC",k+1)
p$samples[[col]] <- as.factor(p$samples[[col]])
map[[k]] <- ggplot(p$samples,aes_string(x="long",y="lat",col=col)) +
  borders(fill="gray22") + theme_bw() +
  geom_point(alpha=0.6) +
  coord_fixed(xlim = c(min(p$samples$long,na.rm=TRUE),max(p$samples$long,na.rm=TRUE)),
              ylim = c(min(p$samples$lat,na.rm=TRUE),max(p$samples$lat,na.rm=TRUE))) +
  scale_color_discrete(name="Region") + theme(legend.position = "non")
}
regions <- cowplot::plot_grid(plotlist = map)
cowplot::save_plot(plot = regions,filename = "plots/kmc_regions.png",dpi=300,base_width=6,base_height=6)

si_fig <- ggplot(p$samples,aes_string(x="long",y="lat",col=paste0("KMC",2))) +
  borders(fill="#F0F0F0") + theme_bw() +
  geom_point(alpha=0.6) +
  coord_fixed(xlim = c(min(p$samples$long,na.rm=TRUE),max(p$samples$long,na.rm=TRUE)),
              ylim = c(min(p$samples$lat,na.rm=TRUE),max(p$samples$lat,na.rm=TRUE))) +
  scale_color_discrete(name="Region",labels = c("DRC East","DRC West","Ghana","Tanzania","Uganda","Zambia")) +
  theme(axis.line = element_blank(),axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.background = element_rect(fill="#9FCAE1"),
        panel.border = element_blank())

cowplot::save_plot(plot = si_fig,filename = "analysis/plots/si_kmc_regions.png",dpi=300,base_width=6,base_height = 6)

si_fig <- ggplot(p$samples[p$samples$Country=="DRC",],aes_string(x="long",y="lat",col=paste0("KMC",2))) +
  borders(fill="#F0F0F0") + theme_bw() +
  geom_point(alpha=0.6) +
  coord_fixed(xlim = c(min(p$samples$long[p$samples$Country=="DRC"],na.rm=TRUE),max(p$samples$long[p$samples$Country=="DRC"],na.rm=TRUE)),
              ylim = c(min(p$samples$lat[p$samples$Country=="DRC"],na.rm=TRUE),max(p$samples$lat[p$samples$Country=="DRC"],na.rm=TRUE))) +
  scale_color_discrete(name="Region",labels = c("DRC East","DRC West","Ghana","Tanzania","Uganda","Zambia")) +
  theme(axis.line = element_blank(),axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.background = element_rect(fill="#9FCAE1"),
        panel.border = element_blank())

cowplot::save_plot(plot = si_fig,filename = "analysis/plots/si_drc_kmc_regions.png",dpi=300,base_width=6,base_height = 6)

## create our mccoil inputs
A1 <- p$counts
A1[is.na(A1)] <- -1
A2 <- p$coverage-p$counts
A2[is.na(A2)] <- -1
wsaf <- MIPanalyzer::get_wsaf(p, impute=FALSE)

save(A1,A2,wsaf,p, file="analysis/data/derived/mccoil_inputs.RData")
save(map,file = "analysis/plots/regions.rds")
