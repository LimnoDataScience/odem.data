library(odem.data)

# setwd('~/Documents/DSI/odem.data/')



library(parallel)
library(MASS)
library(factoextra)
library(cluster)
library(tidyverse)
library(LakeMetabolizer)
library(dtw)
library(zoo)
library(patchwork)
library(rnaturalearth)
# library(mapview)
library(data.table)
library(mlr)
library(nnet)
library(glmulti)
library(performance)# checks and compares quality of models
library(effects)
library(flextable)
library(vegan)
library(pROC)
library(multiROC)


data <- read_csv('analysis/figures/data_nov18.csv', col_names = T)

world <- ne_countries(scale = "medium", returnclass = "sf")
us <- map_data("state")

hydLakes <- read_sf(dsn = "inst/extdata/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp")
data$RT <- hydLakes$Res_time[match(data$Hylak_id, hydLakes$Hylak_id)]

gmap <- ggplot() +
  geom_polygon(data = us, aes(x=long, y=lat,
                              group = group), color = "black", fill = 'white',
               size =.5, alpha = 0) +
  geom_point(data = data, aes(longitude, latitude, col = lndu, shape = ct, size = depth)) +
  coord_sf(xlim = c(-98, -84), ylim = c(42, 50), expand = FALSE) +
  xlab('Longitude') + ylab('Latitude') +
  theme_light(); gmap
ggsave(file = 'analysis/figures/map.png', gmap, dpi = 600, width =15, height = 10)

# check Lake Mendota
data[which(data$lake == 'nhdhr_143249470'),]

df_data = data[,c("ct", "developed", "forest", "cultivated", "wetlands",
                  "depth", "area", "elev", "eutro", "dys", "oligo", "RT", "WshA",
                  "dep_avg", "trophic")]
df_data_num = df_data[, c("developed", "forest", "cultivated", "wetlands",
                  "depth", "area", "elev", "eutro", "dys", "oligo", "RT", "WshA",
                  "dep_avg")]
df_data_num  = apply(df_data_num, 2, function(x)  scale(x))

df_data = cbind(df_data_num, data[,c("ct","trophic")])


normalize <- function(x) {
  num <- x - min(x)
  denom <- max(x) - min(x)
  return (num/denom)
}


smp_size_1 <- floor(0.75 * nrow(df_data[which(df_data$ct == 'Bad'),]))
smp_size_2 <- floor(0.75 * nrow(df_data[which(df_data$ct == 'Good'),]))
# smp_size_3 <- floor(0.75 * nrow(data[which(data$ct == 'SemiBad'),]))

## set the seed to make your partition reproducible
set.seed(123)
train_ind_1 <- sample(seq_len(nrow(df_data[which(df_data$ct == 'Bad'),])), size = smp_size_1)
train_ind_2 <- sample(seq_len(nrow(df_data[which(df_data$ct == 'Good'),])), size = smp_size_2)
# train_ind_3 <- sample(seq_len(nrow(data[which(data$ct == 'SemiBad'),])), size = smp_size_3)

train <- df_data[c(train_ind_1,train_ind_2), ]
test <- df_data[-c(train_ind_1,train_ind_2), ]


# https://jkzorz.github.io/2020/04/04/NMDS-extras.html

ext = df_data[, c("developed", "forest", "cultivated", "wetlands",
                  "depth", "area", "RT", "WshA",
                  "dep_avg", "eutro", "oligo")]
ext = df_data[, c('developed' , 'cultivated' , 'depth' , 'oligo' ,'eutro' , 'RT')]
int = df_data[, c("ct")]
ext$trophic <- as.numeric(as.factor(ext$trophic))

ext_com = as.matrix(ext)

#nmds code
set.seed(123)
nmds = metaMDS(ext, distance = "euclidean", k =2)
nmds

en = envfit(nmds, int, permutations = 999, na.rm = TRUE)
en

plot(nmds)
plot(en)

data.scores = as.data.frame(scores(nmds))
data.scores$trophic <- data$ct

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = data.scores, aes(colour = trophic), size = 3, alpha = 0.5)

adonis2(as.numeric(as.factor(df_data$ct)) ~ cultivated + depth + oligo, data = df_data,
        permutations = 999, method = "euclidean", by = "margin")

g.nmds_1d <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = data.scores, aes(colour = trophic), size = 3, alpha = 0.5)# +
  # # scale_colour_manual(values = c("orange", "steelblue", "red4"))  +
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
  #              data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  # geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2),
  #            shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  # geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04),
  #           label = row.names(en_coord_cat), colour = "navy", fontface = "bold") +
  # geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30",
  #           fontface = "bold", label = row.names(en_coord_cont)) +
  # theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"),
  #       panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"),
  #       axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
  #       legend.title = element_text(size = 10, face = "bold", colour = "grey30"),
  #       legend.text = element_text(size = 9, colour = "grey30")) +
  # labs(colour = "Landuse")

g.nmds_2d <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = data.scores, aes(colour = trophic), size = 3, alpha = 0.5)# +
  # # scale_colour_manual(values = c("orange", "steelblue", "red4"))  +
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
  #              data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  # geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2),
  #            shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  # geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04),
  #           label = row.names(en_coord_cat), colour = "navy", fontface = "bold") +
  # geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30",
  #           fontface = "bold", label = row.names(en_coord_cont)) +
  # theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"),
  #       panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"),
  #       axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
  #       legend.title = element_text(size = 10, face = "bold", colour = "grey30"),
  #       legend.text = element_text(size = 9, colour = "grey30")) +
  # labs(colour = "Landuse")

g.nmds_3d <- ggplot(data = data.scores, aes(x = NMDS2, y = NMDS3)) +
  geom_point(data = data.scores, aes(colour = trophic), size = 3, alpha = 0.5) #+
  # # scale_colour_manual(values = c("orange", "steelblue", "red4"))  +
  # geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
  #              data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  # geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2),
  #            shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  # geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04),
  #           label = row.names(en_coord_cat), colour = "navy", fontface = "bold") +
  # geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30",
  #           fontface = "bold", label = row.names(en_coord_cont)) +
  # theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"),
  #       panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"),
  #       axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
  #       legend.title = element_text(size = 10, face = "bold", colour = "grey30"),
  #       legend.text = element_text(size = 9, colour = "grey30")) +
  # labs(colour = "Landuse")

g.nmds = g.nmds_1d + g.nmds_2d + g.nmds_3d; g.nmds
ggsave(file = 'analysis/figures/nmds.png', g.nmds, dpi = 600, width =14, height = 5)



# library(sjPlot)

glmulti(ct   ~ lndu + developed + forest + cultivated + wetlands + depth + area + elev + eutro + dys + oligo + log10(RT) + log10(WshA) + log10(dep_avg) +
          trophic,
        data   = data,
        # crit   = aicc,       # AICC corrected AIC for small samples
        level  = 1,          # 2 with interactions, 1 without
        method = "d",        # "d", or "h", or "g"
        # family = gaussian,
        fitfunction = multinom,   # Type of model (LM, GLM etc.)
        confsetsize = 100)   # Keep 100 best models

glmulti(ct   ~ developed + forest + cultivated  + depth + area +  eutro + dys + oligo,
        data   = data,
        # crit   = aicc,       # AICC corrected AIC for small samples
        level  = 2,          # 2 with interactions, 1 without
        method = "g",        # "d", or "h", or "g"
        # family = gaussian,
        fitfunction = multinom,   # Type of model (LM, GLM etc.)
        confsetsize = 100)   # Keep 100 best models

test_h <- glmulti(as.factor(ct)   ~ developed + forest + cultivated + wetlands + depth + area + elev +
                    eutro + dys + oligo + (RT) + (WshA) + (dep_avg) + as.factor(trophic),# + trophic,
                  data   = df_data,
                  # crit   = aicc,       # AICC corrected AIC for small samples
                  level  = 1,          # 2 with interactions, 1 without
                  method = "h",        # "d", or "h", or "g"
                  family = binomial,
                  fitfunction = glm,   # Type of model (LM, GLM etc.)
                  confsetsize = 100)   # Keep 100 best models

test_h <- glmulti(as.factor(ct)   ~ developed  + cultivated  + depth  +
                    eutro + dys + oligo + (RT)  ,# + trophic,
                  data   = df_data,
                  # crit   = aicc,       # AICC corrected AIC for small samples
                  level  = 1,          # 2 with interactions, 1 without
                  method = "h",        # "d", or "h", or "g"
                  family = binomial,
                  fitfunction = glm,   # Type of model (LM, GLM etc.)
                  confsetsize = 100)   # Keep 100 best models

optimal_model_glmulti_exhaustive <- test_h@objects[[1]]
print(optimal_model_glmulti_exhaustive)



png(file='analysis/figures/effects.png', width = 25, height = 20, units = "cm",res = 600)
plot(effects::allEffects(test_h@objects[[6]]),
     lines = list(multiline = T),
     confint = list(style = "auto"))
dev.off()

plot(test_h)

weightable(test_h)[1:6,] %>%
  regulartable() %>%       # beautifying tables
  autofit()

plot(test_h, type = "s")

best_model <- test_h@objects[[2]]

car::Anova(best_model)


model <- multinom(ct ~ developed  + cultivated  + depth  +
                    eutro + dys + oligo + (RT) , data = train)

model1 <- multinom(ct ~  developed + depth + eutro, data = df_data)
model2 <- multinom(ct ~  developed + depth + oligo, data = df_data)
model3 <- multinom(ct ~  developed + cultivated + depth + eutro, data = df_data)
model4 <- multinom(ct ~  developed + cultivated  + depth + eutro, data = df_data)
model5 <- multinom(ct ~  developed + cultivated + depth, data = df_data)
model6 <- multinom(ct ~  developed + depth + eutro + RT, data = df_data)

pscl::pR2(model1)["McFadden"]
pscl::pR2(model2)["McFadden"]
pscl::pR2(model3)["McFadden"]
pscl::pR2(model4)["McFadden"]
pscl::pR2(model5)["McFadden"]
pscl::pR2(model6)["McFadden"]

summary(model)
summary(model)$AIC
z <- summary(model)$coefficients/summary(model)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2
print(p)

# https://datasciencebeginners.com/2018/12/20/multinomial-logistic-regression-using-r/
# https://stats.oarc.ucla.edu/r/dae/multinomial-logistic-regression/
## extracting coefficients from the model and exponentiate
exp(coef(model))

head(probability.table <- fitted(model))

# Predicting the values for train dataset
train$precticed <- predict(model, newdata = train, "class")

# Building classification table
ctable <- table(train$ct, train$precticed)

# Calculating accuracy - sum of diagonal elements divided by total obs
round((sum(diag(ctable))/sum(ctable))*100,2)

# Predicting the values for train dataset
test$precticed <- predict(model, newdata = test, "class")

# Building classification table
ctable <- table(test$ct, test$precticed)

# Calculating accuracy - sum of diagonal elements divided by total obs
round((sum(diag(ctable))/sum(ctable))*100,2)

pscl::pR2(model)["McFadden"]



# predicted data
prediction <- predict(model, test, type="probs")


# create roc curve
roc_object <- multiclass.roc( test$ct, prediction)



truerel <-data.frame(var = test$ct, value = 1) %>%
  rownames_to_column() %>%
  pivot_wider(names_from = var, values_from = value) %>%
  replace_na(list(Bad = 0,
                  Good = 0,
                  #  Convex = 0,
                  SemiBad = 0)) %>%
  select(-rowname) %>%
  rename(Bad_true = Bad,
         Good_true = Good,
         # Convex_true = Convex,
         SemiBad_true = SemiBad) %>%
  data.frame() %>%
  bind_cols(., prediction %>%
              data.frame() %>%
              select(Bad, Good, SemiBad) %>%
              rename(Bad_pred_1 = Bad,
                     Good_pred_1 = Good,
                     #   Convex_pred_1 = Convex,
                     SemiBad_pred_1 = SemiBad))

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
