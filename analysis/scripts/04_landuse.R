library(odem.data)

setwd('~/Documents/DSI/odem.data/')


library(parallel)
library(MASS)


library(ggplot2)
library(RColorBrewer)
library(patchwork)

lstm <- read.csv('lstm/all_observations_results.csv')
lstm.sites <- read.csv('lstm/per_site_results.csv')

lake.link <- read.csv('../landuse/lake_link.csv')
lake.information <- read.csv('../landuse/lake_information.csv')
lake.landuse <- read.csv('../landuse/zone_landuse1.csv')

# landuse_lakeatlas <- read.csv('../landuse/lts_atlas_glcp_waste_bathy_glev.csv')
# landuse_nlcd <- read.csv('../landuse/lts_nadp_nlcd_glcp_waste_bathy_glev.csv')

landuse <- data.frame('water' = lake.landuse$nlcd_openwater11_pct + lake.landuse$nlcd_icesnow12_pct,
                      'developed' =  lake.landuse$nlcd_devopen21_pct + lake.landuse$nlcd_devlow22_pct +
                        lake.landuse$nlcd_devmed23_pct + lake.landuse$nlcd_devhi24_pct,
                      'barren' = lake.landuse$nlcd_barren31_pct,
                      'forest' = lake.landuse$nlcd_fordec41_pct + lake.landuse$nlcd_forcon42_pct +
                        lake.landuse$nlcd_formix43_pct,
                      'shrubland' = lake.landuse$nlcd_shrub52_pct,
                      'herbaceous' = lake.landuse$nlcd_grass71_pct ,
                      'cultivated' = lake.landuse$nlcd_past81_pct + lake.landuse$nlcd_cultcrop82_pct,
                      'wetlands' = lake.landuse$nlcd_wetwood90_pct + lake.landuse$nlcd_wetemerg95_pct)

id.of.interest <- gsub(".*_","",lstm.sites$site_id)
id.model.information <- match(id.of.interest, lake.information$lake_nhdid)

id.of.landuse <- match(lake.information$ws_zoneid[id.model.information], lake.landuse$zoneid)



landuse.of.interest <- landuse[id.of.landuse,]

landuse.df <- data.frame('nhdr_id' = lstm.sites$site_id)
landuse.df <- cbind(landuse.df, landuse.of.interest)
landuse.df$LANDUSE <- colnames(landuse.df[2:9])[apply(landuse.df[,2:9],1,which.max)]

write.csv(landuse.df, file = 'lstm/landuse.csv', quote = F, row.names = F)

landuse.df[match('nhdhr_143249470', landuse.df$nhdr_id),]

mendota <-  "143249470"
lake.information[ match(mendota, lake.information$lake_nhdid), ]

trout <- "69886228"
lake.information[ match(trout, lake.information$lake_nhdid), ]
