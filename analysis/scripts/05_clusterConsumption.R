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
morph <- data.frame('lake' = NULL,
                    'depth' = NULL,
                    'area' = NULL,
                    'evel' = NULL)
for (i in lake.list){
  if (file.exists(paste0('analysis/',i,'/modeled_o2.RData')) == FALSE){
    next
  }
  load(paste0('analysis/',i,'/modeled_o2.RData'))# load('Allequash/Allequash.Rda')
  data <- o2$df_kgml

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
# kmediods
distMatrix <- dist(anoxDym2, method= 'euclidean')

hc <- hclust(distMatrix, method='ward.D')
plot(hc, main='')
groups <- cutree(hc, k=4) #k=5) # cut tree into 7 clusters

rect.hclust(hc, k=4,border='red')#k=7
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
names(df) = c('Weak','Linear','Constant','Convex')

nameVec = names(df)
df$depth = seq(1,nrow(df))
table(groups)
lakeinv <- nameVec[unique(which(table(groups) > 0))]# >5

table(groups)[1]

# Pivot wide for plotting
df.long = df %>%
  dplyr::select(lakeinv, depth) %>%
  pivot_longer(lakeinv) %>%
  mutate(name = fct_relevel(name,'Weak','Linear','Constant','Convex'))#'Linear','Constant','Convex'))

# Cluster lables
cluster.labels = NA
# order = match(lakeinv, c('Hypoxic','Linear','Convex')) # Order in the same way as
order = match(lakeinv, c('Weak','Linear','Constant','Convex'))#c('Linear','Constant','Convex')) # Order in the same way as
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
df.grd <-  setNames(data.frame(matrix(ncol = 1+length(seq(1979, 2018,1)), nrow = 207)), c('lake',
                                                                                        as.character(seq(1979,2018,1))))
df.grd$lake <- lake.list

for (i in 1:4) {
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
# m.df.grd$lake <- factor(m.df.grd$lake, levels= rev(c("Allequash","BigMuskellunge","Crystal","Sparkling", "Trout","Fish","Mendota","Monona")))
g1 <- ggplot(m.df.grd, aes(x = variable, y = lake, fill = as.factor(value))) +
  scale_fill_manual(values = c('gold','lightblue1','red4','green'), name = 'Cluster',
                    breaks = c('Weak','Linear','Constant','Convex')) + #c('Linear','Constant','Convex')) +
  geom_tile(color = 'black', width = 0.8, height = 0.8, size = 0.5) +
  labs(x = 'Time', y = '') +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, size = 5),
        axis.title.x = element_blank()); g1

g <- g.cluster / g1 + plot_layout(heights = c(1.5,2))  + plot_annotation(tag_levels = 'A'); g
# g <- g.cluster / g1 + plot_annotation(tag_levels = 'A') + plot_layout(heights = c(1.5,2)); g
ggsave(file = 'analysis/figures/cluster.png', g.cluster, dpi = 500, width =6, height = 4)
ggsave(file = 'analysis/figures/cannual.png', g1, dpi = 600, width =30, height = 20)


types = m.df.grd %>%
  group_by(lake) %>%
  mutate(type = (factor(value))) %>%
  summarize(ct = names(which.max(table(type)))) # constant convex linear

landuse = read_csv('lstm/landuse.csv')
types$lndu <- factor(landuse$LANDUSE[match(types$lake,landuse$nhdr_id)])
types$ct = factor(types$ct)

types$developed <- landuse$developed[match(types$lake,landuse$nhdr_id)]
types$forest <- landuse$forest[match(types$lake,landuse$nhdr_id)]
types$cultivated <- landuse$cultivated[match(types$lake,landuse$nhdr_id)]
types$wetlands <- landuse$wetlands[match(types$lake,landuse$nhdr_id)]

morph$depthf <- cut(morph$depth, c(0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 80, 250, Inf))
morph$areaf <- cut(morph$area, c(0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e9))

types$depthf = morph$depthf[match(types$lake,morph$lake)]
types$areaf = morph$areaf[match(types$lake,morph$lake)]

types$depth = log10(morph$depth[match(types$lake,morph$lake)])
types$area = log10(morph$area[match(types$lake,morph$lake)])

write_csv(x = types, file = 'analysis/figures/model.csv', col_names = T)

# 187 lakes
## 75% of the sample size
troph <- read_csv('analysis/figures/ts_oxygen_join.csv')
link <- read_csv('analysis/figures/nhd_hydrolakes_nla.csv')

types$eutro <- troph$eu.mixo[match(types$lake,troph$site_id)]
types$dys <- troph$dys[match(types$lake,troph$site_id)]
types$oligo <- troph$oligo[match(types$lake,troph$site_id)]

data = na.omit(types)
smp_size <- floor(0.75 * nrow(data))

## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind, ]
test <- data[-train_ind, ]


library(nnet)


model <- multinom(ct ~ dys * lndu * depth * area * eutro * oligo , data = train)

summary(model)

pscl::pR2(model)["McFadden"]

library(pROC)
library(multiROC)

# predicted data
prediction <- predict(model, test, type="probs")

# create roc curve
roc_object <- multiclass.roc( test$ct, prediction)



truerel <-data.frame(var = test$ct, value = 1) %>%
  rownames_to_column() %>%
  pivot_wider(names_from = var, values_from = value) %>%
  replace_na(list(Linear = 0,
                  Constant = 0,
                  Convex = 0,
                  Weak = 0)) %>%
  select(-rowname) %>%
  rename(Linear_true = Linear,
         Constant_true = Constant,
         Convex_true = Convex,
         Weak_true = Weak) %>%
  data.frame() %>%
  bind_cols(., prediction %>%
              data.frame() %>%
              select(Linear, Constant, Convex, Weak) %>%
              rename(Linear_pred_1 = Linear,
                     Constant_pred_1 = Constant,
                     Convex_pred_1 = Convex,
                     Weak_pred_1 = Weak))

roc_object <- multi_roc(data.frame(truerel))

plot_roc_df <- plot_roc_data(roc_object)

aucs <- plot_roc_df %>%
  select(AUC, Method, Group) %>%
  filter(!Group %in% c('Micro','Macro'))%>%
  distinct()

print(aucs)

g <- plot_roc_df %>%
  filter(!Group %in% c('Micro','Macro'))%>%
  ggplot(., aes(x=1-Specificity,y=Sensitivity,color = Method)) +
  geom_step() +
  # geom_text(data=aucs[aucs$Method=='logit',], aes(x=.6,y=.3, label=sprintf('AUC = %.3f',AUC))) +
  # geom_text(data=aucs[aucs$Method=='mlp',], aes(x=.6,y=.5, label=sprintf('AUC = %.3f',AUC))) +
  # geom_text(data=aucs[aucs$Method=='xg',], aes(x=.6,y=.7, label=sprintf('AUC = %.3f',AUC))) +
  scale_color_viridis_d(end=.8) +
  facet_wrap(~Group) +
  theme_bw(); g
ggsave(file = 'analysis/figures/auc.png', g, dpi = 600, width =6, height = 3)


# calculate area under curve
auc( roc_object )

plot(roc_object)
