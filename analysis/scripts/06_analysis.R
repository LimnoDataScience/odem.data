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
library(sf)
library(ggExtra)
library(ggmosaic)
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

hydLakes <- read_sf(dsn = "inst/extdata/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp")
# data$RT <- hydLakes$Res_time[match(data$Hylak_id, hydLakes$Hylak_id)]

lake_shapes <- st_read("inst/extdata/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp")

idy <- (match(data$Hylak_id,lake_shapes$Hylak_id))
lakes_df <- lake_shapes[idy,]
lakes_df$ct <- data$ct
lakes_df$lndu <- data$lndu

ggplot(lakes_df, aes(fill = ct)) +
  theme_minimal() +
  geom_sf() +
  scale_fill_brewer(type = "qual")

gmap <- ggplot(lakes_df, aes(fill = ct)) +
  geom_sf() +
  geom_polygon(data = us, aes(x=long, y=lat,
                              group = group), color = "black", fill = 'white',
               size =.5, alpha = 0) +
  scale_fill_manual(values= c('red4', 'lightblue1')) +
  # geom_point(data = data, aes(longitude, latitude, col = lndu, shape = ct, size = depth)) +
  coord_sf(xlim = c(-97.3, -86), ylim = c(42.6, 48.7), expand = FALSE) +
  xlab('Longitude') + ylab('Latitude') +
  theme_minimal(); gmap

## get model performance
all.dne <- list.files('analysis/')
all.dne_all <- all.dne[grepl('nhdhr', all.dne)]

info.df <- c()
for (idx in all.dne_all){
  if (file.exists(paste0('analysis/',idx,'/lakeinfo.txt'))){
    info.df <- rbind(info.df,read.csv(paste0('analysis/',idx,'/lakeinfo.txt')))
  }
}

data$fit = info.df$fit_tall[na.omit(match(data$lake,info.df$lake_id))]

g.rmse.depth <- ggplot(data, aes(fit, depth, col = trophic)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values= c('darkgreen', 'gold')) +
  xlab("RMSE (g/m3)") + ylab('Depth (m)') +
  theme_minimal(); g.rmse.depth
# g.rmse.depth <- ggMarginal(g.rmse.depth)

g.wsh.area <- ggplot(data, aes(log10(WshA), log10(area), col = trophic)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values= c('darkgreen', 'gold')) +
  xlab("Watershed area (log10 m2)") + ylab('Lake area (log10 m2)') +
  theme_minimal(); g.wsh.area
# g.wsh.area <- ggMarginal(g.wsh.area)

g.mos <- ggplot(data) +
  geom_mosaic( aes( x = product(ct, trophic), fill = ct)) +
  xlab("RMSE (g/m3)") + ylab('Depth (m)') +
  scale_fill_manual(values= c('red4', 'lightblue1')) +
  theme_minimal()+
  theme(legend.position = 'none'); g.mos

g.lndu.wsh <- ggplot(data, aes(lndu, log10(WshA))) +
  geom_boxplot() +
  xlab("Land use") + ylab('Watershed area (log10 m2)') +
  geom_jitter(color = 'black', size = 0.4, alpha = 0.9) +
  theme_minimal(); g.lndu.wsh

g.density <- ggplot(data, aes(depth, fill = ct)) +
  geom_density(alpha = 0.5) +
  xlab("Depth (m)") + ylab('Density (-)') +
  scale_fill_manual(values= c('red4', 'lightblue1')) +
  theme_minimal() +
  theme(legend.position = 'none'); g.density

g <- ggplot(data, aes(depth, fill = trophic)) +
  geom_density(alpha = 0.25) +
  xlab("Depth (m)") + ylab('Density (-)') +
  theme_minimal(); g

fig1 <- gmap / (g.rmse.depth + g.density + g.mos + g.wsh.area) + plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A')

ggsave(file = 'analysis/figures/Figure2.png', fig1, dpi = 600, width =13, height = 15)

df_data = data[,c("ct", "developed", "forest", "cultivated", "wetlands", "water", "barren",
                  "shrubland", "herbaceous",
                  "depth", "area", "elev", "eutro", "dys", "oligo", "RT", "WshA",
                  "dep_avg", "trophic")]
df_data_num = df_data[, c("developed", "forest", "cultivated", "wetlands", "water", "barren",
                          "shrubland", "herbaceous",
                  "depth", "area", "elev", "eutro", "dys", "oligo", "RT",
                  "dep_avg")]

data$ct <- as.numeric(as.factor(data$ct)) - 1

data_new <- data %>%
  mutate(human_impact = developed + cultivated,
         vegetated = forest + wetlands + shrubland,
         depth = log10(depth))


models_exhaust <- glmulti(ct ~ human_impact + log10(area) + log10(depth) +
          eutro + log10(RT),
        data   = data_new,
        # crit   = aicc,       # AICC corrected AIC for small samples
        level  = 2,          # 2 with interactions, 1 without
        method = "h",        # "d", or "h", or "g"
        family = "binomial",
        fitfunction = glm,   # Type of model (LM, GLM etc.)
        confsetsize = 100)   # Keep 100 best models

model_averaged <- model.avg(object = models_exhaust@objects[c(1:24)])



png(file='analysis/figures/Figure3.png', width = 25, height = 20, units = "cm",res = 600)
plot(effects::allEffects(test_h@objects[[1:2]]),
     lines = list(multiline = T),
     confint = list(style = "auto"))
dev.off()

plot(test_h)

weightable(test_h)[1:10,] %>%
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
train$precticed <- predict(model_averaged, newdata = train, "class")

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
data_new$prediction <- stats::predict(model_averaged, type = "response")


# create roc curve
roc_object <- pROC::roc(data_new$ct, data_new$prediction)

train_fraction <- 0.7
reduced_model_data_interation <- 100
# Run the specified number of model averaging iterations
for(rdmi in 1:reduced_data_model_iteration){

  # Define data subset for modeling and sample a fraction for model training
  temp_lake_dat <- data_new %>%
    sample_frac(size = train_fraction,
                weight = as.factor(ct))


  AIC <- rep(0, length(model_result@formulas[c(1:24)]))
  MODEL <- rep(NA, length(model_result@formulas[c(1:24)]))
  AUC <- rep(0, length(model_result@formulas[c(1:24)]))
  RSQUARED <- rep(0, length(model_result@formulas[c(1:24)]))
  ACCURACY <- rep(0, length(model_result@formulas[c(1:24)]))
  PVALUE <- rep(0, length(model_result@formulas[c(1:24)]))
  for(i in 1:length(model_result@formulas[c(1:24)])){
    fit <- glm(paste(as.character(model_result@formulas[i])),
               data = temp_lake_dat,
               family = binomial)
    MODEL[i] <- paste(as.character(model_result@formulas[i]))
    AIC[i] <- fit$aic
    predictpr <- predict(fit, type = "response")
    ROC <- pROC::roc(temp_lake_dat$ct ~ predictpr)
    temp_lake_dat$PREDICTION <- predictpr
    AUC[i] <- pROC::auc(ROC)
    RSQUARED[i] <- 1 - (fit$deviance/fit$null.deviance)

    temp_lake_dat_permute <- temp_lake_dat %>%
      dplyr::mutate(PREDICTION = ifelse(as.numeric(PREDICTION) < 0.65, 0, 1))
    table <- table(Reality = temp_lake_dat_permute$ct,
                   Prediction = temp_lake_dat_permute$PREDICTION)
    ACCURACY[i] <- (table[1,1]+table[2+2])/sum(table)

    j <- 1
    nreps <- 1000
    AUC.repo <- rep(0, nreps)

    for(j in 1:nreps) {
      temp_lake_dat_permute$ct <- sample(temp_lake_dat_permute$ct,
                                                   size = length(temp_lake_dat_permute$ct),
                                                   replace = FALSE)
      predictpr <- predict(fit, type = "response")
      ROC <- pROC::roc(temp_lake_dat_permute$ct ~ predictpr, quiet = TRUE,
                       levels = c(0,1), direction = "<")
      AUC.repo[j] <- pROC::auc(ROC)
    }
    PVALUE[i] <- length(AUC.repo[AUC.repo > AUC[i]])/length(AUC.repo)

  }
  INDEX <- seq(1:length(model_result@formulas[c(1:24)]))
  lake_fits <- data.frame(INDEX, MODEL, AIC, RSQUARED, AUC, ACCURACY, PVALUE)
  lake_fits$MODEL <- as.character(lake_fits$MODEL)
  lake_fits$AIC <- as.numeric(lake_fits$AIC)
  lake_fits$RSQUARED <- as.numeric(lake_fits$RSQUARED)
  lake_fits$AUC <- as.numeric(lake_fits$AUC)
  lake_fits$ACCURACY <- as.numeric(lake_fits$ACCURACY)
  lake_fits$PVALUE <- as.numeric(lake_fits$PVALUE)

  lake_top_mod$MODEL <- gsub(pattern = "log00", replacement = "log10",
                             x = lake_top_mod$MODEL)

  lake_top_mod <- lake_fits %>%
    filter(AIC <= (min(AIC)+2)) %>%
    filter(RSQUARED >= median(RSQUARED),
           AUC >= median(AUC))

  lake_mod_fits <- map(.x = lake_top_mod$MODEL,
                      .f = ~ glm(formula = .x,
                                 family = "binomial",
                                 data = temp_lake_dat))

}

