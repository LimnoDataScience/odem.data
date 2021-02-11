library(odem.data)

setwd('~/Documents/DSI/odem.data/')

all.nml <- list.files('inst/extdata/pball_nml/pball_nml_files/')
all.nml_short <- sub(".*pball_", "", all.nml)
all.nml_short2 <- sub(".nml.*", "", all.nml_short)

library(parallel)
library(MASS)

# for (lake_id in all.nml_short2){
#   source("analysis/scripts/01_data_merge_sql.R")
# }

run_all <- function(x){
  lake_id <<- x
  source("analysis/scripts/01_data_merge_sql.R")
}

numCores <- detectCores()

# system.time(
#   results <- mclapply(all.nml_short2, run_all, mc.cores = numCores)
# )

all.dne <- list.files('analysis/')
all.dne_all <- all.dne[grepl('nhdhr', all.dne)]

idx <- match(all.dne_all, all.nml_short2)
idy <- seq(1:length(all.nml_short2))

if (length(idx) == 0){
  all.nml_missing <- all.nml_short2[idy]
} else {
  all.nml_missing <- all.nml_short2[idy[-c(idx)]]
}


system.time(
  results <- mclapply(all.nml_missing, run_all, mc.cores = numCores)
)
