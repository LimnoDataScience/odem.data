        library(odem.data)

        setwd('~/Documents/DSI/odem.data/')


        if (!exists("password")){
          password <- as.character(read.delim('analysis/scripts/password.txt', header = FALSE, stringsAsFactor = FALSE))
        }

        all.nml <- list.files('inst/extdata/pball_nml/pball_nml_files/')
        all.nml_short <- sub(".*pball_", "", all.nml)
        all.nml_short2 <- sub(".nml.*", "", all.nml_short)

        library(parallel)
        library(MASS)
        library(factoextra)
        library(cluster)


        all.dne <- list.files('analysis/')
        all.dne_all <- all.dne[grepl('nhdhr', all.dne)]


        library(tidyverse)
        library(LakeMetabolizer)
        library(dtw)
        library(zoo)
        library(patchwork)


        anoxDym <- list()
        anoxLab <- c()
        a=1
        lake.list <- all.dne_all
        morph <- data.frame('lake' = NULL,
                            'depth' = NULL,
                            'area' = NULL,
                            'evel' = NULL)
        for (i in lake.list){
          if (file.exists(paste0('analysis/',i,'/modeled_o2.RData')) == FALSE){
            next
          }
          info <- read.csv(paste0('analysis/',i,'/lakeinfo.txt'))
          if (info$fit_tall > 3 | info$n_obs <= 10){ #fit_train, 5
            next
          }
          load(paste0('analysis/',i,'/modeled_o2.RData'))# load('Allequash/Allequash.Rda')
          data <- o2$df_kgml

          if (max(data$o2_hyp) > 20000){
                  next
          }

          nml = glmtools::read_nml(paste0('inst/extdata/pball_nml/pball_nml_files/pball_',i,'.nml'))

          morph <- rbind(morph, data.frame('lake'=i,
                                           'depth'=mean(data$max.d),
                                           'area' = mean(data$area_surface),
                                           'elev' = max(nml$morphometry$H)))

          for (an in unique(data$year)){
            dataAnn = data[ which(data$year %in% an),]
            dataStrat = dataAnn[which(dataAnn$strat == 1),]
            if ((max(dataStrat$doy) - min(dataStrat$doy)) <2){
              next
            }
            dataStrat$Sat_hypo <- o2.at.sat.base(temp = dataStrat$temperature_hypo , altitude = max(nml$morphometry$H)) * 1000
            dataStrat$dist <- (dataStrat$doy - min(dataStrat$doy)) / (max(dataStrat$doy) - min(dataStrat$doy))
            dataStrat$normalDO <- dataStrat$o2_hyp / dataStrat$Sat_hypo
            scaledDOdiff <- approx(dataStrat$dist, dataStrat$normalDO , seq(0,1,0.01))

            if (max(scaledDOdiff$y) > 2){
                    print(i)
                    print(dataStrat$normalDO )
                    print(dataStrat$o2_hyp)
                    print(dataStrat$Sat_hypo)
                    print(scaledDOdiff$y)
                    readline(prompt="Press [enter] to continue")
            }

            anoxDym[[a]] = scaledDOdiff$y
            a=a+1
            anoxLab <- append(anoxLab, paste0(i,'_',an))
          }
        }

        anoxDym2 = lapply(anoxDym,rollmean,k = 10)

        mydata = (as.matrix(as.data.frame(anoxDym2)))
        wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))

        for (i in 2:10) wss[i] <- sum(kmeans(mydata,
                                             centers=i)$withinss)
        plot(1:10, wss, type='b', xlab='Number of Clusters',
             ylab='Within groups sum of squares')

        df = mydata
        g.elb <- fviz_nbclust(df, kmeans, method = 'wss')
        print(g.elb)
        avg_sil <- function(k) {
          km.res <- kmeans(df, centers = k, nstart = 25)
          ss <- silhouette(km.res$cluster, dist(df))
          mean(ss[, 3])
        }
        # Compute and plot wss for k = 2 to k = 15
        k.values <- 2:15
        # extract avg silhouette for 2-15 clusters
        avg_sil_values <- map_dbl(k.values, avg_sil)
        plot(k.values, avg_sil_values,
             type = "b", pch = 19, frame = FALSE,
             xlab = "Number of clusters K",
             ylab = "Average Silhouettes")
        g.sil <- fviz_nbclust(df, kmeans, method = "silhouette")
        print(g.sil)

        distMatrix <- dist(anoxDym2, method= 'euclidean')
        hc <- hclust(distMatrix, method='ward.D')
        plot(hc, main='')
        groups <- cutree(hc, k=2) #k=5) # cut tree into 7 clusters

        rect.hclust(hc, k=2,border='red')#k=7
        indMycl = unique(groups)
        dataGroups = list()
        idz = as.numeric(table(groups))
        for (i in 1:length(indMycl)){
          idy = which(groups == i)
          data = anoxDym2[idy]

          data.df = as.data.frame(do.call(cbind, data))
          dataGroups[[i]] = apply(data.df,1, mean)

          data.long = data.df %>%
            mutate(x = 1:92) %>%
            pivot_longer(-x)

          # Individual Cluster Plots
          p1 = ggplot(data.long, aes(x, value, colour=name)) +
            geom_line() +
            theme(legend.position = "none") +
            labs(title = paste0('Cluster = ',i,' n = ',idz[i]))
          print(p1)
        }

        df = as.data.frame(dataGroups)
        names(df) = c("Heavy consumption", "Low consumption")#c('Semi-bad','Good','Bad')#c('Semi-bad','Good','Bad')#,'Convex')

        nameVec = names(df)
        df$depth = seq(1,nrow(df))
        table(groups)
        lakeinv <- nameVec[unique(which(table(groups) > 0))]# >5

        table(groups)[1]

        # Pivot wide for plotting
        df.long = df %>%
          dplyr::select(lakeinv, depth) %>%
          pivot_longer(lakeinv) %>%
          mutate(name = fct_relevel(name,  "Heavy consumption", "Low consumption"))

        # Cluster lables
        cluster.labels = NA

        order = match(lakeinv, c("Heavy consumption", "Low consumption"))
        for (i in 1:3){
          j = order[i]
          cluster.labels[j] = paste0(lakeinv[i],' (n = ',table(groups)[i],')')
        }

        # Cluster plotting
        g.cluster = ggplot(df.long) +
          geom_line(aes(depth, value, color = name)) +
          scale_color_manual(values = c('red4','lightblue1','gold','red1','red4'), name = 'Cluster',
                             labels = cluster.labels) +
          xlab('Stratification duration [%]') +
          ylab('Ratio of Hypolimnion to \nSaturation DO [-]') +
          theme_minimal(base_size = 8) +
          theme(axis.title.y = element_text(vjust = -10)); g.cluster

        # Grid Plot
        df.grd <-  setNames(data.frame(matrix(ncol = 1+length(seq(1979, 2018,1)), nrow = 207)), c('lake',
                                                                                                  as.character(seq(1979,2018,1))))
        df.grd$lake <- lake.list

        for (i in 1:3) {
          dff = anoxLab[which(groups == i)]
          name = lakeinv[i]
          dmlst <- lake.list[!is.na(match(lake.list, unlist(substr(dff,1,nchar(dff)-5))))]
          for (j in dmlst){
            xrow = which(df.grd$lake == j)
            whpatt = grep(j, dff)
            whyrs = gsub(".*_","",dff[grep(j, dff)])
            df.grd[xrow, !is.na(match(colnames(df.grd),whyrs))] <- name
          }
        }

        m.df.grd <- reshape2::melt(df.grd, id.vars = 'lake')

        g1 <- ggplot(m.df.grd, aes(x = variable, y = lake, fill = as.factor(value))) +
          scale_fill_manual(values = c('red4','lightblue1','gold','red1','red4'), name = 'Cluster',
                            breaks = c("Heavy consumption", "Low consumption")) +
          geom_tile(color = 'black', width = 0.8, height = 0.8, size = 0.5) +
          labs(x = 'Time', y = '') +
          theme_minimal(base_size = 8) +
          theme(axis.text.x = element_text(angle = 45, size = 5),
                axis.title.x = element_blank()); g1

        ggplot(m.df.grd, aes(x = as.numeric(variable), y = as.numeric(as.factor(value)), col = as.factor(lake))) +
          # scale_fill_manual(values = c('red4','gold','lightblue1','red1','red4'), name = 'Cluster',
                            # breaks = c('SemiBad','Good','Bad')) +
          geom_line()+
          labs(x = 'Time', y = '') +
          theme_minimal(base_size = 8) +
          theme(axis.text.x = element_text(angle = 45, size = 5),
                axis.title.x = element_blank(),
                legend.position = 'none')

        g <- g.cluster / g1 + plot_layout(heights = c(1.5,2))  + plot_annotation(tag_levels = 'A'); g
        ggsave(file = 'analysis/figures/cluster.png', g.cluster, dpi = 500, width =6, height = 4)
        ggsave(file = 'analysis/figures/cannual.png', g1, dpi = 600, width =30, height = 20)


        types = m.df.grd %>%
          group_by(lake, variable) %>%
          mutate(type = (factor(value))) %>%
          rename(year = variable) %>%
          summarize(ct = names(which.max(table(type)))) # constant convex linear

        landuse = read_csv('lstm/landuse.csv')
        types$lndu <- factor(landuse$LANDUSE[match(types$lake,landuse$nhdr_id)])
        types$ct = factor(types$ct)

        types$developed <- landuse$developed[match(types$lake,landuse$nhdr_id)]/100
        types$forest <- landuse$forest[match(types$lake,landuse$nhdr_id)]/100
        types$cultivated <- landuse$cultivated[match(types$lake,landuse$nhdr_id)]/100
        types$wetlands <- landuse$wetlands[match(types$lake,landuse$nhdr_id)]/100
        types$water <- landuse$water[match(types$lake,landuse$nhdr_id)]/100
        types$barren <- landuse$barren[match(types$lake,landuse$nhdr_id)]/100
        types$shrubland <- landuse$shrubland[match(types$lake,landuse$nhdr_id)]/100
        types$herbaceous <- landuse$herbaceous[match(types$lake,landuse$nhdr_id)]/100

        morph$depthf <- cut(morph$depth, c(0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 80, 250, Inf))
        morph$areaf <- cut(morph$area, c(0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e9))

        types$depthf = morph$depthf[match(types$lake,morph$lake)]
        types$areaf = morph$areaf[match(types$lake,morph$lake)]

        types$depth = (morph$depth[match(types$lake,morph$lake)])
        types$area = (morph$area[match(types$lake,morph$lake)])

        types$elev = (morph$elev[match(types$lake,morph$lake)])

        # 187 lakes
        ## 75% of the sample size
        # troph_old <- read_csv('analysis/figures/ts_oxygen_join.csv')
        # troph <- read_csv('analysis/figures/complete_lake_predictions.csv')
        troph <- read_csv("lstm/ensemble_preds.csv")
        link <- read_csv('analysis/figures/nhd_hydrolakes_nla.csv')

        types$Hylak_id <- link$Hylak_id[match(types$lake, link$site_id)]

        troph_df <- troph %>%
                rename(eutro = `mean_prob_prob_eu/mixo`,
                       oligo = mean_prob_prob_oligo,
                       dys = mean_prob_prob_dys) %>%
                select(Hylak_id, year, eutro, oligo, dys)

        types <- merge(types, troph_df, by = c("Hylak_id", "year"))

        # sum(unique(match(troph$Hylak_id, link$Hylak_id)), na.rm = T)
        # # 975470004
        # length(unique(link$site_id[match(troph$Hylak_id, link$Hylak_id)]))
        # # 39777
        # sum(!is.na(match(types$lake,link$site_id[match(troph$Hylak_id, link$Hylak_id)])))
        # # 170
        # sum(!is.na( link$Hylak_id[match(types$lake,link$site_id[match(troph$Hylak_id, link$Hylak_id)])]))
        # # 18
        # types$site_id <- link$Hylak_id[match(types$lake,link$site_id[match(troph$Hylak_id, link$Hylak_id)])]
        # types$site_id <- link$Hylak_id[match(link$site_id[match(troph$Hylak_id, link$Hylak_id)],types$lake)]

        # types$eutro <- troph$`mean_prob_prob_eu/mixo`[match(types$Hylak_id,troph$Hylak_id)]
        # types$dys <- troph$mean_prob_prob_dys[match(types$Hylak_id,troph$Hylak_id)]
        # types$oligo <- troph$mean_prob_prob_oligo[match(types$Hylak_id,troph$Hylak_id)]

        # residence times
        library(sf)
        hydLakes <- read_sf(dsn = "inst/extdata/HydroLAKES/HydroLAKES_points_v10_shp/HydroLAKES_points_v10.shp")

        types$RT <- (hydLakes$Res_time[match(troph$Hylak_id[match(types$Hylak_id,troph$Hylak_id)], hydLakes$Hylak_id)])
        types$WshA <- (hydLakes$Wshd_area[match(troph$Hylak_id[match(types$Hylak_id,troph$Hylak_id)], hydLakes$Hylak_id)])
        types$dep_avg <- (hydLakes$Depth_avg[match(troph$Hylak_id[match(types$Hylak_id,troph$Hylak_id)], hydLakes$Hylak_id)])

        col <- read_csv('analysis/figures/limnosat_redux_raw_rel_reflectance_ptl_color.csv')
        # types$NDVI <- (col$Nir_raw[match(troph$Hylak_id[match(types$lake,troph$site_id)], col$Hylak_id)] - col$Red_raw[match(troph$Hylak_id[match(types$lake,troph$site_id)], col$Hylak_id)]) /
        #   (col$Nir_raw[match(troph$Hylak_id[match(types$lake,troph$site_id)], col$Hylak_id)] + col$Red_raw[match(troph$Hylak_id[match(types$lake,troph$site_id)], col$Hylak_id)])



        data = as.data.frame(na.omit(types))

        ################
        library(RPostgreSQL)
        library(rnaturalearth)
        # library(mapview)
        require('RPostgreSQL')
        # dbDisconnect(con)
        drv <- dbDriver("PostgreSQL") ##this is where you call the name of your ODBC connection

        ##create connection (I use the server id here but the url name works too)
        con <- dbConnect(drv,
                         dbname = "UW_data_science",
                         host = '144.92.62.199',
                         port = 5432,
                         user = "postgres",
                         password = password)
        ###try and make map of all lakes from Jordan's data
        lake_metrics <- dbGetQuery(con,'select * from data.lake_metrics', stringsAsFactors = FALSE)
        lake_metrics_reduced <- lake_metrics[!is.na(match(lake_metrics$nhd_lake_id, data$lake)),]


        # data$longitude <- lake_metrics$longitude[!is.na(match(lake_metrics$nhd_lake_id, data$lake))]
        # data$latitude <- lake_metrics$latitude[!is.na(match(lake_metrics$nhd_lake_id, data$lake))]

        data <- data %>%
          mutate(trophic = case_when(eutro > oligo & eutro > dys  ~ 'eutro', # & eutro >= 0.75
                                     oligo > eutro & oligo > dys ~ 'oligo', #  & oligo  >= 0.75
                                     dys > eutro & dys > oligo   ~ 'dys')) %>% #& dys  >= 0.75
          mutate(trophic = ifelse(is.na(trophic), 'gray', trophic))

        write_csv(x = types, file = 'analysis/figures/rawdata_jan19.csv', col_names = T)
        write_csv(x = data.frame('nhdhr' = types$lake), file = 'analysis/figures/my_nhdhr.csv', col_names = T)
        write_csv(x = data, file = 'analysis/figures/data_jan19.csv', col_names = T)

        str(data)
        head(data)
        dim(data)
        summary(data)


        # NEW LANDUSE

        # landuse_nlcd <- read.csv('../landuse/lts_nadp_nlcd_glcp_waste_bathy_glev.csv')
        # landuse <- data.frame('year' = landuse_nlcd$year,
        #                       'Hylak_id' = landuse_nlcd$Hylak_id,
        #                       'developed' =  landuse_nlcd$developed_high_intensity + landuse_nlcd$developed_med_intensity +
        #                               landuse_nlcd$developed_low_intensity + landuse_nlcd$developed_open_space,
        #
        #                       'cultivated' = landuse_nlcd$cultivated_crops)
        #



        # # make plot Figure 2
        # world <- ne_countries(scale = "medium", returnclass = "sf")
        # us <- map_data("state")
        #
        # hydLakes <- read_sf(dsn = "inst/extdata/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp")
        # # data$RT <- hydLakes$Res_time[match(data$Hylak_id, hydLakes$Hylak_id)]
        #
        # lake_shapes <- st_read("inst/extdata/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp")
        #
        # idy <- (match(data$Hylak_id,lake_shapes$Hylak_id))
        # lakes_df <- lake_shapes[idy,]
        # lakes_df$ct <- data$ct
        # lakes_df$lndu <- data$lndu
        #
        # gmap <- ggplot(lakes_df, aes(fill = ct)) +
        #         geom_sf() +
        #         geom_polygon(data = us, aes(x=long, y=lat,
        #                                     group = group), color = "black", fill = 'white',
        #                      size =.5, alpha = 0) +
        #         scale_fill_manual(values= c('red4', 'lightblue1')) +
        #         # geom_point(data = data, aes(longitude, latitude, col = lndu, shape = ct, size = depth)) +
        #         coord_sf(xlim = c(-97.3, -86), ylim = c(42.6, 48.7), expand = FALSE) +
        #         xlab('Longitude') + ylab('Latitude') +
        #         theme_minimal(); gmap
        #
        # g.1 <- g.cluster / gmap + plot_layout(guides = 'collect') +
        #         plot_annotation(tag_levels = 'A') + plot_layout(heights = c(1,3))
        # # ggsave(file = 'analysis/figures/Fig1.png', g.1,  dpi = 600, width =13, height = 15)



        ### FIGURE 1

        data <- read_csv('analysis/figures/data_jan19.csv', col_names = T)

        ## get model performance
        all.dne <- list.files('analysis/')
        all.dne_all <- all.dne[grepl('nhdhr', all.dne)]

        info.df <- c()
        for (idx in all.dne_all){
                if (file.exists(paste0('analysis/',idx,'/lakeinfo.txt'))){
                        info.df <- rbind(info.df,read.csv(paste0('analysis/',idx,'/lakeinfo.txt')))
                }
        }

        data$fit_train = info.df$fit_train[na.omit(match(data$lake,info.df$lake_id))]
        data$fit_test = info.df$fit_test[na.omit(match(data$lake,info.df$lake_id))]
        data$fit_all = info.df$fit_tall[na.omit(match(data$lake,info.df$lake_id))]

        mean_train = data %>%
                pull(fit_train) %>%
                mean() %>%
                signif(3)

        mean_test = data %>%
                pull(fit_test) %>%
                mean(na.rm = T) %>%
                signif(3)

        mean_all = data %>%
                pull(fit_all) %>%
                mean(na.rm = T) %>%
                signif(3)

        plot1 <- ggplot(data) +
                geom_density(aes(x = fit_train, fill = '1_Calibration'), alpha = 0.3) +
                geom_density(aes(x = fit_test, fill = '2_Validation'), alpha = 0.3) +
                geom_density(aes(x = fit_all, fill = '3_Total'), alpha = 0.3) +
                geom_vline(xintercept= mean_train, size=1., col = c("#00AFBB"), linetype = 'dashed') +
                annotate(geom = 'text', x = 3, y = 0.8, label = paste0("Calibration: ",mean_train,
                                                                       '\nValidation: ',mean_test,
                                                                       '\nTotal: ',mean_all),
                         hjust = 0) +
                geom_vline(xintercept= mean_test, size=1.,  col = c("#E7B800"), linetype = 'dashed') +
                geom_vline(xintercept= mean_all, size=1.,  col = c("#FC4E07"), linetype = 'dashed') +
                scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"), name = '') +
                theme_minimal(base_size = 15) +
                xlab('RMSE (g m-3)') + ylab('Density')

        plot4 <- ggplot(df.long) +
                geom_line(aes(depth, value, color = name), size = 1.5) +
                scale_color_manual(values = c('red4','lightblue1','gold','red1','red4'), name = 'Cluster',
                                   labels = cluster.labels) +
                xlab('Stratification duration [%]') +
                ylab('Ratio of Hypolimnion to \nSaturation DO [-]') +
                theme_minimal(base_size = 15)

        # unique(data$lake[which(data$lake[which(data$fit_train < 2)] %in% data$lake[which(data$ct == 'Heavy consumption')])])
        # unique(data$lake[which(data$lake[which(data$fit_train < 2)] %in% data$lake[which(data$ct == 'Low consumption')])])

        id_high = 'nhdhr_120018107' #'nhdhr_120018097'
        id_low = 'nhdhr_120020350'#"nhdhr_120018089"

        # name_high = data %>% filter(lake == id_high) %>% select(Hylak_id)
        # match(unique(name_high),lake_shapes$Hylak_id)
        # lake_shapes %>% filter(Hylak_id > 1041660 & Hylak_id < 1041665)

        plot_time <- function(id, time1, time2, main){
                load(paste0('analysis/',id,'/modeled_o2.RData'))

                o2_data = o2$df_kgml
                idx = which(duplicated(colnames(o2_data)) == T)
                o2_data <- o2_data[, - idx]
                ggplot(o2_data %>% filter(year >= time1 & year <= time2)) +
                        geom_line(aes(datetime, o2_epi/1000, col = 'Epilimnion sim.')) +
                        geom_point(aes(datetime, obs_epi, col = 'Epilimnion obs.'), size = 2) +
                        geom_line(aes(datetime, o2_hyp/1000, col = 'Hypolimnion sim.')) +
                        geom_point(aes(datetime, obs_hyp, col = 'Hypolimnion obs.'), size = 2) +
                        geom_point(aes(datetime, obs_tot, col = 'Total obs.'), size = 2) +
                        # facet_wrap(~year) +
                        ylab(expression("Conc. [g DO"*~m^{-3}*"]")) +
                        scale_color_manual(values = c('red1','red4','lightblue3','lightblue1','gold')) +
                        xlab('') + ggtitle(main) +
                        ylim(c(-2,25)) +  theme_minimal()+
                        theme(legend.text = element_text(size = 11), axis.text.x= element_text(size = 20), plot.title = element_text(size = 20),
                              axis.text.y= element_text(size = 20), text = element_text(size = 20), legend.title = element_blank(), strip.text =element_text(size = 20),
                              legend.position = 'bottom') +
                        guides(col=guide_legend(nrow=3,byrow=TRUE))
        }
        # 2 3 4
        for (id_test in unique(data$lake[which(data$lake[which(data$fit_train < 2)] %in% data$lake[which(data$ct == 'Low consumption')])])){
                p <- plot_time(id_test, 1979, 2018, paste0('Low consumption: ', id_test))
                print(p)
                readline(prompt="Press [enter] to continue")
        }

        plot2 <- plot_time(id_high, 2010, 2013, paste0('High consumption: ', id_high))
        plot3 <- plot_time(id_low, 2005, 2008,  paste0('Low consumption: ', id_low))



        to_plot <- plot1 / (plot2 / plot3 + plot_layout(guides = 'collect') )/ plot4 +
                plot_annotation(tag_levels = 'A')
        ggsave(file = 'analysis/figures/Figure1.png', to_plot,  dpi = 600, width =13, height = 15)



