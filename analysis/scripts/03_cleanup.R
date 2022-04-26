library(feather)
library(odem.data)

setwd('~/Documents/DSI/odem.data/')



all.dne <- list.files('analysis/')
all.dne_all <- all.dne[grepl('nhdhr', all.dne)]

info.df <- c()
for (idx in all.dne_all){
  if (file.exists(paste0('analysis/',idx,'/lakeinfo.txt'))){

    load(paste0('analysis/',idx,'/modeled_o2.RData'))
    feath <- try(df <- read_feather(paste0('analysis/',idx,'/',idx,'.feather')))

    if (grepl("Error", feath[1], fixed = TRUE)){
      df <- o2$df_kgml
      df <- df[!duplicated(colnames(df))]
    }

    # obs.df <- data.frame(df$obs_epi, df$obs_hyp, df$obs_tot)
    # good.rows <- c()
    # for (ix in 1:nrow(obs.df)){
    #   if (all(is.na(obs.df[ix,]))){
    #     good.rows <- append(good.rows, ix)
    #   }
    # }
    # good.obs.df <- obs.df[-good.rows,]
    #
    # smp_size = floor(0.7 * nrow(good.obs.df))
    #
    # df$splitsample <- NA
    # df$splitsample[c(as.numeric(rownames(good.obs.df)[1:smp_size]))] <- 0
    # df$splitsample[c(as.numeric(rownames(good.obs.df)[(smp_size+1):(nrow(good.obs.df))]))] <- 1

    df$o2_epi <- df$o2_epi / 1000
    df$o2_hyp <- df$o2_hyp / 1000
    df$o2_tot <- df$o2_tot / 1000


    write_feather(df, paste0('analysis/',idx,'/',idx,'.feather'))


  }
}

df <- read_feather(paste0('analysis/',idx,'/',idx,'.feather'))
colnames(df)
