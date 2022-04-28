library(odem.data)

# setwd(dirname(rstudioapi::getSourceEditorContext()$path))



all.nml <- list.files('inst/extdata/pball_nml/pball_nml_files/')
all.nml_short <- sub(".*pball_", "", all.nml)
all.nml_short2 <- sub(".nml.*", "", all.nml_short)

library(parallel)
library(MASS)


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
for (i in lake.list){
  load(paste0('analysis/',i,'/modeled_o2.RData'))# load('Allequash/Allequash.Rda')
  data <- o2$df_kgml

  elev = 400

  for (an in unique(data$year)){
    dataAnn = data[ which(data$year %in% an),]
    dataStrat = dataAnn[which(dataAnn$strat == 1),]
    if ((max(dataStrat$doy) - min(dataStrat$doy)) <2){
      next
    }
    dataStrat$Sat_hypo <- o2.at.sat.base(temp = dataStrat$temperature_hypo , altitude = elev) * 1000
    dataStrat$dist <- (dataStrat$doy - min(dataStrat$doy)) / (max(dataStrat$doy) - min(dataStrat$doy))
    dataStrat$normalDO <- dataStrat$o2_hyp / dataStrat$Sat_hypo
    scaledDOdiff <- approx(dataStrat$dist, dataStrat$normalDO , seq(0,1,0.01))

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
distMatrix <- dist(anoxDym2, method= 'euclidean')

hc <- hclust(distMatrix, method='ward.D')
plot(hc, main='')
groups <- cutree(hc, k=3) #k=5) # cut tree into 7 clusters

rect.hclust(hc, k=3,border='red')#k=7
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
names(df) = c('Linear','Constant','Convex')

nameVec = names(df)
df$depth = seq(1,nrow(df))
table(groups)
lakeinv <- nameVec[unique(which(table(groups) > 0))]# >5

table(groups)[1]

# Pivot wide for plotting
df.long = df %>%
  dplyr::select(lakeinv, depth) %>%
  pivot_longer(lakeinv) %>%
  mutate(name = fct_relevel(name,'Linear','Constant','Convex'))

# Cluster lables
cluster.labels = NA
# order = match(lakeinv, c('Hypoxic','Linear','Convex')) # Order in the same way as
order = match(lakeinv, c('Linear','Constant','Convex')) # Order in the same way as
for (i in 1:5){
  j = order[i]
  cluster.labels[j] = paste0(lakeinv[i],' (n = ',table(groups)[i],')')
}

# Cluster plotting
g.cluster = ggplot(df.long) +
  geom_line(aes(depth, value, color = name)) +
  scale_color_manual(values = c('gold','lightblue1','red4','red1','red4'), name = 'Cluster',
                     labels = cluster.labels) +
  xlab('Stratification duration [%]') +
  ylab('Ratio of Hypolimnion to \nSaturation DO [-]') +
  theme_minimal(base_size = 8) +
  theme(axis.title.y = element_text(vjust = -10)); g.cluster

# Grid Plot
df.grd <-  setNames(data.frame(matrix(ncol = 1+length(seq(1979, 2018,1)), nrow = 201)), c('lake',
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
m.df.grd$lake <- factor(m.df.grd$lake, levels= rev(c("Allequash","BigMuskellunge","Crystal","Sparkling", "Trout","Fish","Mendota","Monona")))
g1 <- ggplot(m.df.grd, aes(x = variable, y = lake, fill = as.factor(value))) +
  scale_fill_manual(values = c('gold','lightblue1','red4'), name = 'Cluster',
                    breaks = c('Linear','Constant','Convex')) +
  geom_tile(color = 'black', width = 0.8, height = 0.8, size = 0.5) +
  labs(x = 'Time', y = '') +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, size = 5),
        axis.title.x = element_blank()); g1

g <- g.cluster / g1 + plot_layout(heights = c(1.5,2))  + plot_annotation(tag_levels = 'A'); g
# g <- g.cluster / g1 + plot_annotation(tag_levels = 'A') + plot_layout(heights = c(1.5,2)); g
ggsave(file = 'cluster.png', g.cluster, dpi = 500, width =6, height = 4)
ggsave(file = 'annual.png', g1, dpi = 600, width =30, height = 20)


## Cluster stratification duration boxplots
# clust.strat.dur <- c()
# for (ix in (unique(na.omit(m.df.grd$value)))){
#   da.clust <- m.df.grd[which(ix == m.df.grd$value),]
#   for (i in 1:nrow(da.clust)){
#     year.x <- strat.dur[which(da.clust$variable[i] == strat.dur$year),]
#     lake.x <- year.x[which(da.clust$lake[i] == year.x$lake),]
#     clust.strat.dur <- rbind(clust.strat.dur, data.frame('cluster' = ix, 'stratdur' = lake.x$dur))
#   }
# }
#
# ggplot(clust.strat.dur, aes(x = cluster, y = stratdur)) +
#   geom_boxplot()
