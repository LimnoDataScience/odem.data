library(odem.data)

# setwd('~/Documents/DSI/odem.data/')

library(MuMIn)
library(broom)
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
  scale_color_manual(values= c('brown', 'darkgreen', 'lightblue')) +
  xlab("RMSE (g/m3)") + ylab('Depth (m)') +
  theme_minimal(); g.rmse.depth
# g.rmse.depth <- ggMarginal(g.rmse.depth)

g.wsh.area <- ggplot(data, aes(log10(WshA), log10(area), col = trophic)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values= c('brown', 'darkgreen', 'lightblue')) +
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
  scale_color_manual(values= c('brown', 'darkgreen', 'lightblue')) +
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
        level  = 1,          # 2 with interactions, 1 without
        method = "h",        # "d", or "h", or "g"
        family = "binomial",
        fitfunction = glm,   # Type of model (LM, GLM etc.)
        confsetsize = 100)   # Keep 100 best models

model_averaged <- model.avg(object = models_exhaust@objects[c(1:2)])

# predicted data
data_new$prediction <- stats::predict(model_averaged, type = "response")

# create roc curve
roc_object <- pROC::roc(data_new$ct, data_new$prediction)

train_fraction <- 0.7
reduced_model_data_interation <- 5000

model_result <- models_exhaust
# Run the specified number of model averaging iterations
for(rdmi in 1:reduced_model_data_interation){

  # Define data subset for modeling and sample a fraction for model training
  temp_lake_dat <- data_new %>%
    sample_frac(size = train_fraction,
                weight = as.factor(ct))

  # model_result <- glmulti(ct ~ human_impact + log10(area) + log10(depth) +
  #                           eutro + log10(RT),
  #                         data   = temp_lake_dat,
  #                         # crit   = aicc,       # AICC corrected AIC for small samples
  #                         level  = 1,          # 2 with interactions, 1 without
  #                         method = "h",        # "d", or "h", or "g"
  #                         family = "binomial",
  #                         fitfunction = glm,   # Type of model (LM, GLM etc.)
  #                         confsetsize = 100)   # Keep 100 best models

  #num_models <- length(model_result@crits[model_result@crits <= min(model_result@crits)+2])

  num_models <- 24

  AIC <- rep(0, length(model_result@formulas[c(1:num_models)]))
  MODEL <- rep(NA, length(model_result@formulas[c(1:num_models)]))
  AUC <- rep(0, length(model_result@formulas[c(1:num_models)]))
  RSQUARED <- rep(0, length(model_result@formulas[c(1:num_models)]))
  ACCURACY <- rep(0, length(model_result@formulas[c(1:num_models)]))
  PVALUE <- rep(0, length(model_result@formulas[c(1:num_models)]))
  for(i in 1:length(model_result@formulas[c(1:num_models)])){
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
  INDEX <- seq(1:length(model_result@formulas[c(1:num_models)]))
  lake_fits <- data.frame(INDEX, MODEL, AIC, RSQUARED, AUC, ACCURACY, PVALUE)
  lake_fits$MODEL <- as.character(lake_fits$MODEL)
  lake_fits$AIC <- as.numeric(lake_fits$AIC)
  lake_fits$RSQUARED <- as.numeric(lake_fits$RSQUARED)
  lake_fits$AUC <- as.numeric(lake_fits$AUC)
  lake_fits$ACCURACY <- as.numeric(lake_fits$ACCURACY)
  lake_fits$PVALUE <- as.numeric(lake_fits$PVALUE)

  lake_top_mod <- lake_fits %>%
    filter(AIC <= (min(AIC)+2)) %>%
    filter(RSQUARED >= median(RSQUARED),
           AUC >= median(AUC))

  lake_top_mod$MODEL <- gsub(pattern = "log00", replacement = "log10",
                             x = lake_top_mod$MODEL)


  lake_mod_fits <- map(.x = lake_top_mod$MODEL,
                      .f = ~ glm(formula = .x,
                                 family = "binomial",
                                 data = temp_lake_dat))

  out_path <- "permute_odem_model"
  # Create export folder if it doesn't exist
  ifelse(!dir.exists(file.path(out_path)),
         dir.create(file.path(out_path), recursive = TRUE), FALSE)

  # Export results
  if(length(lake_mod_fits) == 1) {

    results <- tidy(lake_mod_fits[[1]])

    write.csv(x = results,
              file = paste0(out_path, "/run_",
                            rdmi, ".csv"),
              row.names = FALSE)

  } else if (length(lake_mod_fits) == 0){

    tryCatch(write.csv(x = results,
                       file = paste0(out_path, "/run_",
                                     rdmi, ".csv"),
                       row.names = FALSE), error = function(e) NULL)

  } else {

    all_average <- model.avg(lake_mod_fits)

    results <- data.frame(summary(all_average)$coefmat.subset) %>%
      rename(estimate = Estimate,
             std.error = Std..Error,
             statistic = z.value,
             p.value = Pr...z..)
    results$term <- rownames(results)

    write.csv(x = results,
              file = paste0(out_path, "/run_",
                            rdmi, ".csv"),
              row.names = FALSE)
  }

}

run_files <- list.files(path = "permute_odem_model",
                        pattern = "run_", full.names = TRUE)

run_results <- map_df(.x = run_files,
                      .f = ~ read_csv(.x) %>%
                        select(estimate, term) %>%
                        mutate(run_number = .x,
                               run_number = gsub(pattern = "permute_odem_model/",
                                                 replacement = "", x = run_number),
                               run_number = gsub(pattern = ".csv",
                                                 replacement = "", x = run_number),
                               term = gsub(pattern = ".1",
                                           replacement = "", x = term)))


paramter_counts <- run_results %>%
  # filter(estimate >= -20,
  #        estimate <= 20) %>%
  unique() %>%
  filter(term != "(Intercept)") %>%
  mutate(term = str_replace(pattern = ":", replacement = ".", string = term),
         term = ifelse(term == "lo0(area)", "area", term),
         term = ifelse(term == "lo0(RT)", "RT", term),
         term = ifelse(term == "lo0(depth)", "depth", term)) %>%
  group_by(term) %>%
  count() #%>%
  #mutate(label = paste("Number of models:", n))

run_results_filtered <- run_results %>%
  ungroup() %>%
  # filter(estimate >= -10,
  #        estimate <= 10) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = str_replace(pattern = ":", replacement = ".", string = term),
         term = ifelse(term == "lo0(area)", "area", term),
         term = ifelse(term == "lo0(RT)", "RT", term),
         term = ifelse(term == "lo0(depth)", "depth", term),
         #term = str_replace(pattern = "", replacement = "", string = term),
         est_prob = exp(estimate)/(1+exp(estimate))) %>%
  unique() %>%
  inner_join(x = .,
             y = paramter_counts) %>%
  mutate(facet_label = paste(term, "(n =", n, ")"))

complete_results <- (model_averaged$coefficients)%>%
  data.frame() %>%
  rownames_to_column() %>%
  filter(rowname == "full") %>%
  rename("depth" = log10.depth.,
         "RT" = log10.RT.,
         "area" = log10.area.) %>%
  pivot_longer(cols = c(human_impact:RT), names_to = "term", values_to = "estimate")

all_plot <- ggplot() +
  geom_histogram(data = run_results_filtered, aes(x = est_prob)) +
  #geom_label(data = paramter_counts, aes(label = label, x = 0, y = 1000)) +
  #geom_vline(data = complete_results, aes(xintercept = estimate)) +
  ggtitle("Distribution of subsampled estimate values") +
  #xlim(-25, 25) +
  facet_wrap(vars(term), scales = "free_x")

ggsave(file = 'analysis/figures/histogram_probs_params.png',
        dpi = 600, height = 6, width = 8)

