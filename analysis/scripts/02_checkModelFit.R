library(odem.data)

setwd('~/Documents/DSI/odem.data/')

all.nml <- list.files('inst/extdata/pball_nml/pball_nml_files/')
all.nml_short <- sub(".*pball_", "", all.nml)
all.nml_short2 <- sub(".nml.*", "", all.nml_short)

library(parallel)
library(MASS)

if (!exists("password")){
  password <- as.character(read.delim('analysis/scripts/password.txt', header = FALSE, stringsAsFactor = FALSE))
}

all.dne <- list.files('analysis/')
all.dne_all <- all.dne[grepl('nhdhr', all.dne)]

info.df <- c()
for (idx in all.dne_all){
  if (file.exists(paste0('analysis/',idx,'/lakeinfo.txt'))){
    info.df <- rbind(info.df,read.csv(paste0('analysis/',idx,'/lakeinfo.txt')))
  }
}

 empt.lst <- print(all.dne_all[!(all.dne_all %in%  info.df$lake_id)])
 print(empt.lst)
# for (inls in empt.lst){
#   print(list.files(paste0('analysis/',inls)))
#   if (length(list.files(paste0('analysis/',inls)))==1 && list.files(paste0('analysis/',inls)) == 'input.txt'){
#     unlink(paste0('analysis/',inls), recursive = T, force = T)
#   }
# }


library(ggplot2)
library(RColorBrewer)
library(patchwork)

g1 <- ggplot(info.df) +
  geom_point(aes(fit_train, log10(area), col = depth, size = n_obs)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylab('log10 Surface Area') + xlab("Training RMSE (70%)") +
  ggtitle(paste0('Mean RMSE: ', round(mean(info.df$fit_train),1), ' g/m3')) +
  theme_bw(); g1

g2 <- ggplot(info.df) +
  geom_point(aes(fit_test, log10(area), col = depth, size = n_obs)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylab('log10 Surface Area') + xlab("Testing RMSE (30%)") +
  ggtitle(paste0('Mean RMSE: ', round(mean(info.df$fit_test, na.rm = T),1), ' g/m3')) +
  theme_bw(); g2

g3 <- ggplot(info.df) +
  geom_point(aes(fit_tall, log10(area), col = depth, size = n_obs)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylab('log10 Surface Area') + xlab("Total Time Period RMSE (100%)") +
  ggtitle(paste0('Mean RMSE: ', round(mean(info.df$fit_tall),1), ' g/m3')) +
  theme_bw(); g3

g4 <- ggplot(info.df) +
  geom_point(aes(fit_train, log10(area), col = depth, size = n_obs)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylab('log10 Surface Area') + xlab("Training RMSE (70%)") +
  xlim(0,10)+
  ggtitle(paste0('Median RMSE: ', round(median(info.df$fit_train),1), ' g/m3')) +
  theme_bw(); g4

g5 <- ggplot(info.df) +
  geom_point(aes(fit_test, log10(area), col = depth, size = n_obs)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylab('log10 Surface Area') + xlab("Testing RMSE (30%)") +
  xlim(0,10)+
  ggtitle(paste0('Median RMSE: ', round(median(info.df$fit_test, na.rm = T),1), ' g/m3')) +
  theme_bw(); g5

g6 <- ggplot(info.df) +
  geom_point(aes(fit_tall, log10(area), col = depth, size = n_obs)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylab('log10 Surface Area') + xlab("Total Time Period RMSE (100%)") +
  xlim(0,10)+
  ggtitle(paste0('Median RMSE: ', round(median(info.df$fit_tall),1), ' g/m3')) +
  theme_bw(); g6

g7 <- ggplot(info.df) +
  geom_point(aes(fit_train, fit_test, col = depth, size = n_obs)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylab("Testing RMSE (30%)") + xlab("Training RMSE (70%)") +
  ggtitle(paste0('Training compared to Testing')) +
  theme_bw(); g7


g <-  (g1 + g2) / (  g3 + g7) +plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A', title = paste0('Simulated lakes: ',nrow(info.df))); g
ggsave('analysis/figures/Figure1.png', g, dpi = 600, width =7, height = 5)

g <- ((g1 / g2 / g3)  | (g4 / g5 / g6) )+  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A', title = paste0('Simulated lakes: ',nrow(info.df))); g
ggsave('analysis/figures/process_run.png', g)

g7 <- ggplot(info.df) +
  geom_density(aes(log10(area))) +
  # geom_point(aes(fit_train, log10(area), col = depth, size = n_obs)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylab('Density') + xlab('log10 Surface Area') +
  ggtitle(paste0('Surface area')) +
  theme_bw(); g7

g8 <- ggplot(info.df) +
  geom_density(aes((depth))) +
  # geom_point(aes(fit_train, log10(area), col = depth, size = n_obs)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylab('Density') + xlab('Depth') +
  ggtitle(paste0('Depth')) +
  theme_bw(); g8

g9 <- ggplot(info.df) +
  geom_density(aes((n_obs))) +
  # geom_point(aes(fit_train, log10(area), col = depth, size = n_obs)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylab('Density') + xlab('Observation points') +
  ggtitle(paste0('Data')) +
  theme_bw(); g9

g10 <- ggplot(info.df) +
  geom_density(aes((fit_tall))) +
  # geom_point(aes(fit_train, log10(area), col = depth, size = n_obs)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylab('Density') + xlab('RMSE') +
  ggtitle(paste0('Total Time Period RMSE (100%)')) +
  theme_bw(); g10

g <- ((g10 + g7) / (g8 + g9) )+  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A', title = paste0('Simulated lakes: ',nrow(info.df))); g
ggsave('analysis/figures/distributions_run.png', g)

ggplot(info.df) +
  geom_point(aes(fit_train, fit_test, col = depth, size = n_obs)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylab("Testing RMSE (30%)") + xlab("Training RMSE (70%)") +
  xlim(0,10) + ylim(0,10) +
  ggtitle(paste0('Mean RMSE: ', round(mean(info.df$fit_tall),1), ' g/m3')) +
  theme_bw()

lstm <- read.csv('lstm/all_observations_results.csv')
lstm.sites <- read.csv('lstm/per_site_results.csv')

idx <- match(lstm.sites$site_id, info.df$lake_id)
lstm.sites$depth = info.df$depth[idx]
lstm.sites$area = info.df$area[idx]
lstm.sites$n_obs = info.df$n_obs[idx]

g1=ggplot(lstm.sites) +
  geom_line(aes(n_obs, lstm_model_rmse,  col = 'LSTM')) +
  geom_point(aes(n_obs, lstm_model_rmse,  col = 'LSTM')) +
  geom_line(aes(n_obs, o2_model_rmse,  col = 'PB-metabolism')) +
  geom_point(aes(n_obs, o2_model_rmse,  col = 'PB-metabolism')) +
  # scale_color_gradient(low = "blue", high = "red") +
  ylab("RMSE (g/m3)") + xlab("Obs") +
  ylim(0,20)+
  ggtitle(paste0('Median Obs: ', ceiling(mean(lstm.sites$n_obs)))) +
  theme_bw()

g2=ggplot(lstm.sites) +
  geom_line(aes(log10(area/1e6), lstm_model_rmse,  col = 'LSTM')) +
  geom_point(aes(log10(area/1e6), lstm_model_rmse,  col = 'LSTM')) +
  geom_line(aes(log10(area/1e6), o2_model_rmse,  col = 'PB-metabolism')) +
  geom_point(aes(log10(area/1e6), o2_model_rmse,  col = 'PB-metabolism')) +
  # scale_color_gradient(low = "blue", high = "red") +
  ylab("RMSE (g/m3)") + xlab("log10 Area (km2)") +
  ylim(0,20)+
  ggtitle(paste0('Median Area: ', floor(mean(lstm.sites$area)/1e6), ' km2')) +
  theme_bw()

g3=ggplot(lstm.sites) +
  geom_line(aes(depth, lstm_model_rmse,  col = 'LSTM')) +
  geom_point(aes(depth, lstm_model_rmse,  col = 'LSTM')) +
  geom_line(aes(depth, o2_model_rmse,  col = 'PB-metabolism')) +
  geom_point(aes(depth, o2_model_rmse,  col = 'PB-metabolism')) +
  # scale_color_gradient(low = "blue", high = "red") +
  ylab("RMSE (g/m3)") + xlab("Depth (m)") +
  ylim(0,20)+
  ggtitle(paste0('Median Depth: ', round(mean(lstm.sites$depth),1), ' m')) +
  theme_bw()


# ggplot(lstm.sites ) +
#   geom_density(aes(o2_model_rmse, fill = 'PB-metabolism'), alpha = 0.1) +
#   geom_density(aes(lstm_model_rmse,  fill = 'LSTM'), alpha = 0.1) +
#   ggtitle('Total average DO depletion') +
#   xlim(0,10)+
#   xlab('average Jz (Livingstone) against average SED+NEP (g/m3/d)') +
# # geom_point(aes(Jv, NEP, col = as.factor(year))) +
# # ylim(min(c(coeff$Jv, coeff$NEP), na.rm = T),max(c(coeff$Jv, coeff$NEP), na.rm = T)) +
# # xlim(min(c(coeff$Jv, coeff$NEP), na.rm = T),max(c(coeff$Jv, coeff$NEP), na.rm = T)) +
# theme_bw()

df <- data.frame('RMSE' = append(lstm.sites$o2_model_rmse, lstm.sites$lstm_model_rmse),
                 'Model' = append(rep('PB-metabolism', length(lstm.sites$o2_model_rmse)),
                                 rep('LSTM', length(lstm.sites$lstm_model_rmse))))
g4=ggplot(df) +
  geom_boxplot(aes(Model, RMSE,  fill = Model))  +
  ggtitle(paste0('Median RMSE ', round(median(lstm.sites$o2_model_rmse),2),' and ', round(median(lstm.sites$lstm_model_rmse),2), ' g/m3 for PB and LSTM')) +
  ylim(0,7.5)+ ylab('RMSE (g/m3)')+
  theme_bw()


g <- ((g1 + g2) / (g3 + g4) )+  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A', title = paste0('Simulated lakes: ',nrow(lstm.sites))); g
ggsave('analysis/figures/comparison_LSTM.png', g)



print(mean(lstm.sites$o2_model_rmse))
print(mean(lstm.sites$lstm_model_rmse))

library(tidyverse)
library(mapview)
library(lattice)
library(sf)
library(lubridate)
library("rnaturalearth")
library("rnaturalearthdata")
library(here)
library(urbnmapr)
library(tmap)


# install.packages('RPostgreSQL')
library(RPostgreSQL)
require('RPostgreSQL')
dbDisconnect(con)

################################################################
# Open the database connection
drv <- dbDriver("PostgreSQL") ##this is where you call the name of your ODBC connection

##create connection (I use the server id here but the url name works too)
con <- dbConnect(drv,
                 dbname = "UW_data_science",
                 host = '144.92.62.199',
                 port = 5432,
                 user = "postgres",
                 password = password)
#################################################################


###try and make map of all lakes from Jordan's data
lake_metrics <- dbGetQuery(con,'select * from data.lake_metrics', stringsAsFactors = FALSE)
lake_metrics_reduced <- lake_metrics[!is.na(match(lake_metrics$nhd_lake_id, info.df$lake_id)),]
lake_metrics_reduced$RMSE = info.df$fit_tall


world <- ne_countries(scale = "medium", returnclass = "sf")
us <- map_data("state")
us.sf <- st_as_sf(us)

gmap <- ggplot() +
  # geom_sf(color = "black") +
  geom_polygon(data = us, aes(x=long, y=lat,
                              group = group), color = "black", fill = 'white',
               size =.5, alpha = 0) +
  geom_point(data = lake_metrics_reduced, aes(longitude, latitude, col = RMSE)) +
  scale_color_gradient(low="darkgreen", high="red") +
  coord_sf(xlim = c(-98, -84), ylim = c(42, 50), expand = FALSE) +
  xlab('Longitude') + ylab('Latitude') +
  theme_minimal(); gmap

g <-( ((g1 + g2 + g3) / (g7 + g8 + g9) )  +  plot_layout(guides = 'collect') ) / (gmap) +
  plot_annotation(tag_levels = 'A', title = paste0('Simulated lakes: ',nrow(info.df))); g
ggsave('analysis/figures/map+plots.png', g, dpi = 300, units = 'cm', width = 80, height = 20)


m = mapview(lake_metrics_reduced, xcol = "longitude", ycol = "latitude", zcol = "RMSE",cex = 2, crs = 4269, at = seq(0,10,1),legend = TRUE,grid = FALSE); m
mapshot(m, file = 'analysis/figures/map.png', remove_controls = c("homeButton", "layersControl"), selfcontained = FALSE)
