## Helper functions

#' Calculate water density from temperature
#'
#' Calculate water density from water temperature using the formula from (Millero & Poisson, 1981).
#'
#' @param wtemp vector or matrix; Water temperatures
#' @return vector or matrix; Water densities in kg/m3
#' @export
calc_dens <- function(wtemp){
  dens = 999.842594 + (6.793952 * 10^-2 * wtemp) - (9.095290 * 10^-3 *wtemp^2) + (1.001685 * 10^-4 * wtemp^3) - (1.120083 * 10^-6* wtemp^4) + (6.536336 * 10^-9 * wtemp^5)
  return(dens)
}


#' Extract time and space information
#'
#' Extracts time (from date column) and space (aka depth) information
#'
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return list of datetimes and depths
#' @export
#' @examples \dontrun{
#' wtemp <- data.frame(
#'   date = c("1979-04-02", "1979-04-02"),
#'   temp_0 = c(0, 0),
#'   temp_1 = c(1, 1),
#'   stringsAsFactors = FALSE)
#' extract_time_space(wtemp)
#' }
extract_time_space <- function(wtemp){
  time <- as.Date(as.character(wtemp$date))
  depth <- sub("^[^_]*_", "", colnames(wtemp[2:ncol(wtemp)]))
  return(list('datetime' = time, 'depth' = depth))
}

#' Calculate thermocline depth
#'
#' Calculate planar thermocline depth by checking the highest density gradient over depth.
#'
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return vector of thermocline depths in m
#' @importFrom stats na.omit
#' @importFrom rLakeAnalyzer thermo.depth center.buoyancy
#' @importFrom zoo na.approx
#' @importFrom lubridate year yday
#' @export
#'
calc_td_depth <- function(wtemp){

  grd.info <- extract_time_space(wtemp)
  temp <- as.matrix(wtemp[,-c(1)])
  dens <- calc_dens(temp)

  cbuoy.depth <- rep(NA, length(grd.info$datetime))
  thermocline_depth <- rep(NA, length(grd.info$datetime))

  condition<- apply(temp, 1, FUN=min,na.rm=TRUE) > 4

  for (ii in 1:length(cbuoy.depth)){
    idx = !is.na(temp[ii,])
    dens_data = dens[ii,idx]
    dens.diff = rev(dens_data)[1] - dens_data[1]

    if (condition[ii] && abs(dens.diff) > 0.05){
    cbuoy.depth[ii] <- center.buoyancy(temp[ii, idx], as.numeric(grd.info$depth[idx]))
    thermocline_depth[ii] <- thermo.depth(temp[ii, idx], as.numeric(grd.info$depth[idx]))
    }
  }

  zdeps <- as.numeric(grd.info$depth)
  wlm.depth <- rep(NA, length(grd.info$datetime))

  for (ii in 1:length(wlm.depth)){
    idx = !is.na(temp[ii,])
    dens_data = dens[ii,idx]
    dens.diff = rev(dens_data)[1] - dens_data[1]

    if (condition[ii] && abs(dens.diff) > 0.05){

      Ch <- rep(NA, length = length(dens_data))
      for (jj in 1:(length(dens_data)-1)){
        Ah = 1/(zdeps[jj+1]) * sum(dens_data[1:jj])
        Bh = 1/(zdeps[length(zdeps)] -zdeps[jj]) * sum(dens_data[(jj + 1): length(dens_data)])

        diffAh = sum( (dens_data[jj:(jj+1)] - Ah)^2 )
        diffBh = sum( (dens_data[jj:(jj+1)] - Bh)^2 )

        Ch[jj] = diffAh + diffBh
      }
      clineDep = zdeps[which.min(na.omit(Ch))]
      wlm.depth[ii] <- clineDep
    }
  }

  # plot(cbuoy.depth)
  # points(thermocline_depth, col='red')
  #
  #
  test <- data.frame('year' = year(grd.info$datetime), 'doy' = yday(grd.info$datetime),
                     'depth' = cbuoy.depth) #wlm.depth

  for (kk in unique(test$year)){
    idx <- which(kk == test$year)

    dx <- test[idx,3]
    dx[which(dx == (max(zdeps-1)))] = NA
    dx[which(dx == 0)] = NA

    NonNAindex <- which(!is.na(dx))
    if (length(na.omit(dx)) != 0){
      firstNonNA <- min(NonNAindex)
      lastNonNA <- max(NonNAindex)
      dx[firstNonNA:lastNonNA] =na.approx(dx)
    }
    test[idx,3] <-  dx
  }
  #

  return(test$depth)#return(cbuoy.depth)#return(test$depth)
}

#'
#' Calculate planar thermocline depth by checking the highest density gradient over depth.
#'
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return vector of thermocline depths in m
#' @importFrom rLakeAnalyzer meta.depths
#' @export
#'
calc_metalimdepth <- function(wtemp){

  grd.info <- extract_time_space(wtemp)
  temp <- as.matrix(wtemp[,-c(1)])
  dens <- calc_dens(temp)

  cbuoy.depth <- rep(NA, length(grd.info$datetime))
  metalimn.depth <- matrix(c(NA,NA), ncol=length(grd.info$datetime), nrow=2)#rep(NA, length(grd.info$datetime))

  condition<- apply(temp, 1, FUN=min,na.rm=TRUE) > 4

  for (ii in 1:length(cbuoy.depth)){
    idx = !is.na(temp[ii,])
    dens_data = dens[ii,idx]
    dens.diff = rev(dens_data)[1] - dens_data[1]

    if (condition[ii] && abs(dens.diff) > 0.05){
      cbuoy.depth[ii] <- center.buoyancy(temp[ii, idx], as.numeric(grd.info$depth[idx]))
      metalimn.depth[,ii] <- meta.depths(temp[ii, idx], as.numeric(grd.info$depth[idx]),
                                        slope = 0.1, seasonal = TRUE, mixed.cutoff = 1)
    }
  }


  #

  return(metalimn.depth)#return(cbuoy.depth)#return(test$depth)
}

#'
#' Calculate mean epilimnion and hypolimnion surface temperature,using the thermocline depth
#'
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @param thermocline_depth matrix; thermocline depth
#' @param H the depth info of the lake
#' @return list of temperatures vector
calc_epil_hypo_temp<-function(wtemp,thermocline_depth,H){
  grd.info <- extract_time_space(wtemp)
  temp <- as.matrix(wtemp[,-c(1)])
  depth_data = as.double(grd.info$depth)

  epil_temp <-  rep(NA, length(thermocline_depth))
  hypo_temp <- rep(NA, length(thermocline_depth))
  total_temp<- rep(NA, length(thermocline_depth))

  total<- rep(NA, length(thermocline_depth))
  hypo<- rep(NA, length(thermocline_depth))
  total<- rep(NA, length(thermocline_depth))

  td_not_exist <- is.na(thermocline_depth)


  for (ii in (1:length(thermocline_depth))){
    idx = !is.na(temp[ii,])
    temp_data = as.numeric(temp[ii,idx])
    #total_temp[ii]<-max(sum(temp_data)/length(temp_data),4)
    total_temp[ii]<-sum(temp_data)/length(temp_data)
    if(!td_not_exist[ii]){
      td_idx <- max(which(thermocline_depth[ii]>=depth_data))
      epil_temp[ii] <- mean(temp_data[1:td_idx])
      if (td_idx >= length(temp_data)) {td_idx = length(temp_data)-1}
      hypo_temp[ii] <- mean(temp_data[(td_idx+1):length(temp_data)])
      # if (is.na(hypo_temp[ii]) && !is.na(epil_temp[ii])){
      #   break
      # }
    }
  }


  return(list('t_epil' = epil_temp,'t_hypo' = hypo_temp,'t_total' = total_temp))

}

#' Calculate water total volume,using the thermocline depth
#'
#' @param H depths
#' @param A areas
#'
#' @return list of temperatures
#' @importFrom pracma trapz
calc_volume_total<-function(H, A){
  if (length(H)==1){
    volume_total <- 1/3.0 * A * H
  }else{
    volume_total <- trapz(rev(H),rev(A))
  }
  return (volume_total)
}

#'
#' Calculate water epilimnion and hypolimnine volume,using the thermocline depth
#'
#' @param H vector; the depth info of the lake
#' @param A vector; the area info of the lake for each depth
#' @param thermocline_depth matrix; thermocline depth
#' @param volume_total number;the total volume of the lake
#' @importFrom stats approx
#' @return matrix of epil and hypo volume
calc_epil_hypo_vol <- function(H,A,thermocline_depth,volume_total){
  vol_data <- matrix(NA, nrow = length(thermocline_depth), ncol = 2)
  colnames(vol_data) <- c("volume_epi","volume_hypo")

  td_not_exist<-is.na(thermocline_depth)
  for (ii in 1:length(thermocline_depth)){
    if(!td_not_exist[ii]){
      h_idx <- min(which(thermocline_depth[ii]>=H))
      approx_td.area<-approx(H, A, c(0, thermocline_depth[ii]))$y[2]
      if (is.na(approx_td.area)){approx_td.area = min(A)}
      H_with_td<-c(thermocline_depth[ii],H[h_idx:length(H)])
      A_with_td<-c(approx_td.area,A[h_idx:length(A)])
      vol_data[ii,1] <- trapz(rev(H_with_td),rev(A_with_td))##epil
      vol_data[ii,2] <- volume_total - vol_data[ii,1]##hypo
      if(vol_data[ii,2] <= 0) {
        vol_data[ii,2]= min(A) * 0.5
        vol_data[ii,1] = volume_total - vol_data[ii,2]
      }
      # if (is.na( vol_data[ii,1]) && !is.na( vol_data[ii,2])){
      #   break
      # }
    }
  }
  return (vol_data)
}


#' Create temperature and volume input values for oxygen model
#'
#' Calculate mean epilimnion and hypolimnion surface temperature, as well as volumes.
#'
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @param H depths
#' @param A areas
#' @return vector of thermocline depths in m
#' @export
#'
input <- function(wtemp, H, A){
  grd.info       <- extract_time_space(wtemp)
  thermocline_depth       <- calc_td_depth(wtemp)
  metalimn.depth <- calc_metalimdepth(wtemp)
  td_area        <- approx(H, A, thermocline_depth)$y
  area_surface      <- rep(max(A), length(grd.info$datetime))
  temp_out       <- calc_epil_hypo_temp(wtemp, thermocline_depth, H) # epi T hypo T
  volume_total      <- calc_volume_total(H, A) # epi V hypo V
  vol            <- calc_epil_hypo_vol(H, A, thermocline_depth, volume_total)
  max.d          <- max(H)

  return(data.frame(
    datetime = as.POSIXct(grd.info$datetime),
    thermocline_depth,
    temperature_epi   = temp_out$t_epil,
    temperature_hypo   = temp_out$t_hypo,
    temperature_total  = temp_out$t_total,
    volume_total, vol,
    td_area, area_surface,
    upper.metalim = metalimn.depth[1,], lower.metalim = metalimn.depth[2,],
    max.d
    ))
}



###

#' Calculate the mass of dissolved oxgen in epil and hypolimnion(during stratified period)
#' and the total dissovled oxygen in the lake
#' @param wtemp matrix; Water temperatures (rows correspond to time, cols to depth)
#' @return list of datetimes and depths
#' @export
#'
calc_do<-function(input.values,fsed_stratified_epi,fsed_stratified_hypo,fsed_not_stratified,nep_stratified,nep_not_stratified,
                  min_stratified, min_not_stratified, wind = NULL, khalf = NULL, startdate = NULL, enddate = NULL){

  ##initialize matrix
  o2_data <- matrix(NA, nrow = length(input.values$thermocline_depth), ncol = 18)
  colnames(o2_data) <- c("o2_epil","o2_hypo","o2_total",
                         'Fsed_total', "NEP_total", "Fatm_total", "Mineral_total",
                         'Fsed_epi', "NEP_epi", 'Fatm_epi', 'Entrain_epi',
                         'Fsed_hypo', 'Mineral_hypo','Entrain_hypo',
                         'sat_o2_epil', 'sat_o2_hypo', 'sat_o2_total', 'massbal')
  o2_data <- as.data.frame(o2_data)

  init_o2sat <- o2.at.sat.base(temp=input.values$temperature_total[1],altitude = 300)*1000 # returns mg O2/L

  theta<-1.08

  if (is.null(khalf)){
    khalf <- 4800.
  }

  td_not_exist <- is.na(input.values$thermocline_depth)

  if (is.null(startdate)){
    startdate = 1
  }
  if (is.null(enddate)){
    enddate =  length(input.values$thermocline_depth)
  }

  o2_data$o2_total[startdate] <- init_o2sat # returns mg O2 (m3 = 1000 L)

  for(day in (startdate + 1):enddate){

    K600<-ifelse(is.null(wind), k.cole.base(2),k.cole.base(wind[day]))

    ## not stratified period, only consider o2_total dynamics
    if(td_not_exist[day]){
      theta_total <- theta^(input.values$temperature_total[day-1]-20)

      NEP <- (nep_not_stratified * theta_total *(o2_data[day-1,"o2_total"])/(khalf + o2_data[day-1,"o2_total"]))

      MINER <- (min_not_stratified * theta_total*(o2_data[day-1,"o2_total"])/(khalf + o2_data[day-1,"o2_total"]))

      kO2 <- k600.2.kGAS.base(k600=K600,temperature=input.values$temperature_total[day-1],gas='O2') # velocity value m/d?
      o2sat<-o2.at.sat.base(temp=input.values$temperature_total[day-1],altitude = 300)*1000 # mg O2/L -> mg/m3

      Fatm <- (kO2*(o2sat - o2_data[day-1,"o2_total"])/max(H))

      Fsed <- (fsed_not_stratified * (o2_data[day-1,"o2_total"])/(khalf + o2_data[day-1,"o2_total"]) * theta_total * input.values$area_surface[day-1] / input.values$volume_total[day-1])

      o2_data[day,"o2_total"] <-( o2_data[day-1,"o2_total"]+ valid((Fsed + NEP + Fatm + MINER),o2_data[day-1,"o2_total"]) ) * ( input.values$volume_total[day-1]/ input.values$volume_total[day])

      if (valid((Fsed + NEP + Fatm + MINER), o2_data[day-1,"o2_total"]) != (Fsed + NEP + Fatm + MINER)){
        Fsed <- Fsed/(Fsed + NEP + Fatm + MINER) * o2_data[day-1,"o2_total"]
        NEP <- NEP/(Fsed + NEP + Fatm + MINER) * o2_data[day-1,"o2_total"]
        Fatm <- Fatm/(Fsed + NEP + Fatm + MINER)* o2_data[day-1,"o2_total"]
        MINER <- MINER/(Fsed + NEP + Fatm + MINER) * o2_data[day-1,"o2_total"]
      }

      o2_data[day,"Fsed_total"] <- Fsed
      o2_data[day,"NEP_total"] <- NEP
      o2_data[day,"Fatm_total"] <- Fatm
      o2_data[day,"Mineral_total"] <- MINER
      o2_data[day,"sat_o2_total"] <- (100. * o2_data[day,"o2_total"] )/ o2sat

    }
    # the day it turns to stratified, need to reassign the o2 to hypo and epil
    else if(is.na(input.values$thermocline_depth[day-1])){
      # o2_data[day,"o2_epil"] <- (o2_data[day-1,"o2_total"]*input.values$volume_epi[day])/input.values$volume_total[day]
      # o2_data[day,"o2_hypo"] <- (o2_data[day-1,"o2_total"]*input.values$volume_hypo[day])/input.values$volume_total[day]
      # o2_data[day,"o2_total"] <- (o2_data[day,"o2_hypo"]* input.values$volume_hypo[day] + o2_data[day,"o2_epil"] * input.values$volume_epi[day])/input.values$volume_total[day]
      o2_data[day,"o2_epil"] <- o2_data[day-1,"o2_total"]#((o2_data[day-1,"o2_total"]*input.values$volume_epi[day])/input.values$volume_total[day]) #/ input.values$volume_epi[day] #/input.values$temperature_total[day]*input.values$volume_epi[day]
      o2_data[day,"o2_hypo"] <- o2_data[day-1,"o2_total"]#((o2_data[day-1,"o2_total"]*input.values$volume_hypo[day])/input.values$volume_total[day])  #/ input.values$volume_hypo[day] #/input.values$temperature_total[day]*input.values$volume_hypo[day]
      o2_data[day,"o2_total"]<- (o2_data[day,"o2_epil"]+o2_data[day,"o2_hypo"]) /2 #/ input.values$volume_total[day]

      o2sat<-o2.at.sat.base(temp=input.values$temperature_epi[day],altitude = 300)*1000
      o2_data[day,"sat_o2_epil"] <- (100. * o2_data[day,"o2_epil"])/ o2sat
      o2sat<-o2.at.sat.base(temp=input.values$temperature_hypo[day],altitude = 300)*1000
      o2_data[day,"sat_o2_hypo"] <- (100. * o2_data[day,"o2_hypo"]  )/ o2sat
      o2sat<-o2.at.sat.base(temp=input.values$temperature_total[day],altitude = 300)*1000
      o2_data[day,"sat_o2_total"] <- (100. * o2_data[day,"o2_total"] )/ o2sat
    }else{
      theta_epil <- theta^(input.values$temperature_epi[day-1]-20)
      theta_hypo <-  theta^(input.values$temperature_hypo[day-1]-20)

      NEP_epil <- (nep_stratified * theta_epil *(o2_data[day-1,"o2_epil"])/(khalf + o2_data[day-1,"o2_epil"]))

      kO2_epil <- k600.2.kGAS.base(k600=K600,temperature=input.values$temperature_epi[day-1],gas='O2')
      o2sat_epil<-o2.at.sat.base(temp=input.values$temperature_epi[day-1],altitude = 300)*1000

      Fatm_epil <- (kO2_epil*(o2sat_epil-o2_data[day-1,"o2_epil"] )/input.values$thermocline_depth[day-1])


      volumechange_epi = input.values$volume_epi[day]-input.values$volume_epi[day-1]
      volumechange_epi_proportion =  volumechange_epi/input.values$volume_epi[day-1]
      if (volumechange_epi_proportion >= 0){
        x_do <- o2_data[day - 1,"o2_hypo"] #( o2_data[day - 1,"o2_hypo"] * abs(volumechange_epi)) / input.values$volume_hypo[day-1]
      } else {
        x_do <- o2_data[day - 1,"o2_epil"] #( o2_data[day - 1,"o2_epil"] * abs(volumechange_epi)) / input.values$volume_epi[day-1]
      }

      Fepi <-  (volumechange_epi_proportion* x_do)

      Fsed_epi <- (fsed_stratified_epi * (o2_data[day-1,"o2_epil"])/(khalf + o2_data[day-1,"o2_epil"]) * theta_epil * input.values$area_surface[day-1] /input.values$volume_epi[day-1])

      o2_data[day,"o2_epil"] <- (o2_data[day-1,"o2_epil"] + valid((NEP_epil+Fatm_epil + Fepi + Fsed_epi) , o2_data[day-1,"o2_epil"]) ) * (input.values$volume_epi[day-1]/input.values$volume_epi[day])

      if (valid((NEP_epil+Fatm_epil + Fepi + Fsed_epi), o2_data[day-1,"o2_epil"]) != (NEP_epil+Fatm_epil + Fepi + Fsed_epi)){
        Fsed_epi <- Fsed_epi/(NEP_epil+Fatm_epil + Fepi + Fsed_epi) * o2_data[day-1,"o2_epil"]
        NEP_epil <- NEP_epil/(NEP_epil+Fatm_epil + Fepi + Fsed_epi) * o2_data[day-1,"o2_epil"]
        Fatm_epil <- Fatm_epil/(NEP_epil+Fatm_epil + Fepi + Fsed_epi) * o2_data[day-1,"o2_epil"]
        Fepi <- Fepi/(NEP_epil+Fatm_epil + Fepi + Fsed_epi) * o2_data[day-1,"o2_epil"]
      }

      o2_data[day,"Fsed_epi"] <- Fsed_epi
      o2_data[day,"NEP_epi"] <- NEP_epil
      o2_data[day,"Fatm_epi"] <- Fatm_epil
      o2_data[day,"Entrain_epi"] <- Fepi

      volumechange_hypo = input.values$volume_hypo[day]-input.values$volume_hypo[day-1]  #in m^3
      volumechange_hypo_proportion =  volumechange_hypo/input.values$volume_hypo[day-1]

      if (volumechange_hypo_proportion >= 0){
        x_do <-  o2_data[day - 1,"o2_epil"] #( o2_data[day - 1,"o2_epil"] * abs(volumechange_hypo)) / input.values$volume_epi[day-1]
      } else {
        x_do <- o2_data[day - 1,"o2_hypo"] #( o2_data[day - 1,"o2_hypo"] * abs(volumechange_hypo)) / input.values$volume_hypo[day-1]
      }

      Fhypo <- (volumechange_hypo_proportion* x_do )

      MINER_hypo <- (min_stratified * theta_hypo  * (o2_data[day-1,"o2_hypo"])/(khalf + o2_data[day-1,"o2_hypo"]))


      Fsed_hypo <- (fsed_stratified_hypo * (o2_data[day-1,"o2_hypo"])/(khalf + o2_data[day-1,"o2_hypo"]) * theta_hypo * input.values$ area_thermocline[day-1] / input.values$volume_hypo[day-1])


      o2_data[day,"o2_hypo"] <- (o2_data[day-1,"o2_hypo"] + valid((Fhypo + Fsed_hypo + MINER_hypo), o2_data[day-1,"o2_hypo"]) ) * (input.values$volume_hypo[day-1]/input.values$volume_hypo[day])

      if (valid((Fhypo + Fsed_hypo + MINER_hypo), o2_data[day-1,"o2_hypo"]) != (Fhypo + Fsed_hypo + MINER_hypo)){
        Fsed_hypo <- Fsed_hypo/(Fhypo + Fsed_hypo + MINER_hypo) * o2_data[day-1,"o2_hypo"]
        Fhypo <- Fhypo/(Fhypo + Fsed_hypo + MINER_hypo) * o2_data[day-1,"o2_hypo"]
        MINER_hypo <- MINER_hypo/(Fhypo + Fsed_hypo + MINER_hypo) * o2_data[day-1,"o2_hypo"]
      }


      o2_data[day,"o2_total"] <- (o2_data[day,"o2_hypo"]* input.values$volume_hypo[day] + o2_data[day,"o2_epil"] * input.values$volume_epi[day])/input.values$volume_total[day]
      o2_data[day,"Fsed_hypo"] <- Fsed_hypo
      o2_data[day,"Mineral_hypo"] <- MINER_hypo
      o2_data[day,"Entrain_hypo"] <- Fhypo


      o2sat<-o2.at.sat.base(temp=input.values$temperature_epi[day],altitude = 300)*1000
      o2_data[day,"sat_o2_epil"] <- (100. * o2_data[day,"o2_epil"] )/ o2sat
      o2sat<-o2.at.sat.base(temp=input.values$temperature_hypo[day],altitude = 300)*1000
      o2_data[day,"sat_o2_hypo"] <- (100. * o2_data[day,"o2_hypo"]  )/ o2sat
      o2sat<-o2.at.sat.base(temp=input.values$temperature_total[day],altitude = 300)*1000
      o2_data[day,"sat_o2_total"] <- (100. * o2_data[day,"o2_total"] )/ o2sat

      # mass balance
      mass_thr <- sum(c(o2_data$o2_epil[day-1] * input.values$volume_epi[day-1], o2_data$o2_hypo[day-1] * input.values$volume_hypo[day-1])) +
        ( (NEP_epil + Fatm_epil + Fepi + Fsed_epi) * (input.values$volume_epi[day-1]) + (MINER_hypo + Fhypo + Fsed_hypo ) * (input.values$volume_hypo[day-1]))
      mass_balance <-  sum(c(o2_data$o2_epil[day]*input.values$volume_epi[day], o2_data$o2_hypo[day]* input.values$volume_hypo[day])) - mass_thr
      o2_data[day,"massbal"] <- mass_balance
      # print(mass_balance)
    }
    if (is.na(o2_data[day, "o2_total"])) {
      break
      print('RED ALERT!')
    }
  }

  return (o2_data)
}


odem_static<-function(input.values,
                      sed,
                      nep,
                      min,
                      wind = NULL,
                      khalf = 500,
                      elev = 450,
                      startdate = NULL, enddate = NULL,
                      field.values){

  ##initialize matrix
  o2_data <- matrix(NA, nrow = length(input.values$thermocline_depth), ncol = 16)
  colnames(o2_data) <- c("o2_epi","o2_hyp","o2_tot",
                         "NEP_mgm3d",
                         "SED_mgm2d",
                         "MIN_mgm3d",
                         'fnep',
                         'fmineral',
                         'fsed',
                         'fatm',
                         'fentr_epi',
                         'fentr_hyp',
                         'sat_o2_epi', 'sat_o2_hyp', 'sat_o2_tot', 'massbal')
  o2_data <- as.data.frame(o2_data)
  ##when is the lake mixed/stratified?
  input.values$strat <- ifelse(is.na(input.values$thermocline_depth),0,1)
  strat.pos <- c()
  for (ii in 1:length(input.values$strat)){
    if (input.values$strat[ii] == 1 && input.values$strat[ii-1] == 0){
      strat.pos <- append(strat.pos, ii)
    }
  }

  theta0 = 1.08^(input.values$temperature_total - 20)
  theta1 = 1.08^(input.values$temperature_epi - 20)
  theta2 = 1.08^(input.values$temperature_hypo - 20)
  k600t = k600.2.kGAS.base(k.cole.base(input.values$wind),temperature = input.values$temperature_total, gas = "O2")
  o2satt = o2.at.sat.base(temp = input.values$temperature_total, altitude = elev) * 1000
  k600 = k600.2.kGAS.base(k.cole.base(input.values$wind),temperature = input.values$temperature_epi, gas = "O2")
  o2sat = o2.at.sat.base(temp = input.values$temperature_epi, altitude = elev) * 1000
  o2sattt = o2.at.sat.base(temp = input.values$temperature_hypo, altitude = elev) * 1000
  volume_epi = input.values$volume_epi
  volume_tot = input.values$volume_total
  area_epi = input.values$area_surface
  volume_hyp = input.values$volume_hypo
  area_hyp = input.values$area_thermocline
  tddepth = input.values$thermocline_depth
  wtr_epi = input.values$temperature_epi
  wtr_hyp = input.values$temperature_hypo
  wtr_tot = input.values$temperature_total
  khalf = khalf # New, was 3000
  DO_epi_init = 15 * 1000 #simdata$DO_obs[1],
  DO_hyp_init = 15 * 1000
  DO_tot_init = 15 * 1000
  stratified = input.values$strat
  strat_pos = strat.pos
  len_strat_pos = length(strat.pos)
  d_strat_pos = length(strat.pos)

  airtemp = input.values$airtemp
  delvol_epi = c(diff(input.values$volume_epi),0)/c(input.values$volume_epi)
  delvol_hyp =  c(diff(input.values$volume_hypo),0)/c(input.values$volume_hypo)
  delvol_epi[strat.pos] = 0
  delvol_hyp[strat.pos] = 0
  delvol_epi[is.na(delvol_epi)] = 0
  delvol_hyp[is.na(delvol_hyp)] = 0
  # simdata$DO_obs_epi = simdata$DO_obs_epi * 1.5
  DO_obs_epi = field.values[3,]
  DO_obs_hyp = field.values[4,]
  DO_obs_tot = field.values[2,]
  k600t[which(airtemp <= 0 & input.values$temperature_total <= 4)] = 1e-5
  mean_depth = volume_tot[1]/area_epi[1];

  delvol = rep(NA, nrow(o2_data))
  delvol[1] = 0;
  for(i in 2:nrow(o2_data)) {
    delvol[i] = volume_epi[i]-volume_epi[i-1];
  }


  o2_data$o2_epi[1] = DO_epi_init;
  o2_data$o2_hyp[1] = DO_hyp_init;
  o2_data$o2_tot[1] = DO_tot_init;
  o2_data$NEP_mgm3d = nep
  o2_data$SED_mgm2d = sed
  o2_data$MIN_mgm3d = min

  first_day = 0;

  for(i in 2:nrow(o2_data)) {
    first_day = 0; #determine if it is first of stratification change
    for (k in 1:len_strat_pos){
      if (strat_pos[k] == i){
        first_day = 1;
      }
    }
    if (stratified[i] == 0) {
      o2_data$fnep[i] =   o2_data$NEP_mgm3d[i] * theta0[i-1];
      o2_data$fmineral[i] = o2_data$MIN_mgm3d[i] * theta0[i-1];
      o2_data$fsed[i] = o2_data$SED_mgm2d[i] *  (max((o2_data$o2_tot[i-1]),1e-06)/(khalf + max((o2_data$o2_tot[i-1]),1e-06))) * theta0[i-1] / mean_depth; # THETA BUG
      o2_data$fatm[i] =  k600t[i-1]  *  (o2satt[i-1] - o2_data$o2_tot[i-1]) / mean_depth; # EPI O2 BUG
      o2_data$o2_tot[i] =  o2_data$o2_tot[i-1] + o2_data$fnep[i] -
              o2_data$fsed[i] + o2_data$fatm[i] + o2_data$fmineral[i];
      o2_data$o2_hyp[i] = o2_data$o2_tot[i] ;
      o2_data$o2_epi[i] = o2_data$o2_tot[i] ;
    } else if (stratified[i] == 1){

      if(first_day != 1){
        if(delvol[i]>0) {
          x_do1 = o2_data$o2_hyp[i-1];
        } else {
          x_do1 = o2_data$o2_epi[i-1];
        }
      }

      if(first_day == 1){
        o2_data$fnep[i] =  o2_data$fnep[i-1];
        o2_data$fatm[i]  = o2_data$fatm[i-1];
        o2_data$o2_epi[i] =  (o2_data$o2_epi[i-1]);
      } else {
        o2_data$fnep[i] =  o2_data$NEP_mgm3d[i] *theta1[i-1];
        o2_data$fatm[i] = k600[i-1] *  (o2sat[i-1] - o2_data$o2_epi[i-1]) / tddepth[i-1]; #if DO is negative...use 0.1 instead of inferred DO
        o2_data$fentr_epi[i] =
        o2_data$o2_epi[i] =  ((o2_data$o2_epi[i-1] + o2_data$fnep[i] + o2_data$fatm[i])*volume_epi[i-1] + (delvol[i]*x_do1))/volume_epi[i];
      }
      if(first_day == 1) {
        o2_data$fsed[i] = o2_data$fsed[i-1];
        o2_data$fmineral[i] = o2_data$fmineral[i-1];
        o2_data$o2_hyp[i] = o2_data$o2_hyp[i-1];
      } else {
        # mineral[i] = MIN[i_Param[i]] * (fmax((DO_hyp[i-1]*tau+mu),1e-06)/(khalf + fmax((DO_hyp[i-1]*tau+mu),1e-06))) * theta2[i-1];
        o2_data$fmineral[i] = o2_data$MIN_mgm3d[i] * theta2[i-1];
        o2_data$fsed[i] = (o2_data$SED_mgm2d[i] *  (max((o2_data$o2_hyp[i-1]),1e-06)/(khalf + max((o2_data$o2_hyp[i-1]),1e-06))) * theta2[i-1]) / max((volume_hyp[i-1]/area_hyp[i-1]),1);
        o2_data$o2_hyp[i] =  ((o2_data$o2_hyp[i-1] - o2_data$fsed[i] + o2_data$fmineral[i])*volume_hyp[i-1] - (delvol[i]*x_do1))/volume_hyp[i];
      }
      o2_data$o2_tot[i] = (o2_data$o2_epi[i] *volume_epi[i] + o2_data$o2_hyp[i] *volume_hyp[i])/volume_tot[i];
    }
  }

  o2_data$sat_o2_epi <- (100. * o2_data$o2_epi )/ o2sat
  o2_data$sat_o2_hyp <- (100. * o2_data$o2_hyp )/ o2sattt
  o2_data$sat_o2_tot <- (100. * o2_data$o2_tot )/ o2satt

  # mass balance
  o2_data$massbal <-( o2_data$o2_epi * volume_epi + o2_data$o2_hyp * volume_hyp) - (
    ( (o2_data$fnep + o2_data$fatm ) * (volume_epi) +
        (o2_data$fsed * (-1) + o2_data$fmineral  ) * (volume_hyp)))


  idx.obs <- match(as.Date(obs_weigh_df$Date),as.Date(input.values$datetime))

  obs <- cbind(as.numeric(obs_weigh_df$Epi[!is.na(idx.obs)]),
               as.numeric(obs_weigh_df$Hyp[!is.na(idx.obs)]))
  mod <- cbind(o2_data$o2_epi[idx.obs[!is.na(idx.obs)]]/1000,
               o2_data$o2_hyp[idx.obs[!is.na(idx.obs)]]/1000)

  fit <- sqrt(mean((obs-mod)**2,na.rm = TRUE))

  o2_data$obs_tot <- NaN
  o2_data$obs_epi <- NaN
  o2_data$obs_hyp <- NaN
  o2_data$obs_tot[idx.obs[!is.na(idx.obs)]] <- as.numeric(obs_weigh_df$Tot[!is.na(idx.obs)])
  o2_data$obs_epi[idx.obs[!is.na(idx.obs)]] <- as.numeric(obs_weigh_df$Epi[!is.na(idx.obs)])
  o2_data$obs_hyp[idx.obs[!is.na(idx.obs)]] <- as.numeric(obs_weigh_df$Hyp[!is.na(idx.obs)])
  o2_data$date <- input.values$datetime
  o2_data$doy <- yday(o2_data$date)
  o2_data$year <- year(o2_data$date)

  plot <- ggplot(o2_data) +
    geom_line(aes(doy, o2_epi/1000, col = 'Epilimnion sim.')) +
    geom_point(aes(doy, obs_epi, col = 'Epilimnion obs.'), size = 2) +
    geom_line(aes(doy, o2_hyp/1000, col = 'Hypolimnion sim.')) +
    geom_point(aes(doy, obs_hyp, col = 'Hypolimnion obs.'), size = 2) +
    geom_point(aes(doy, obs_tot, col = 'Total obs.'), size = 2) +
    facet_wrap(~year) +
    ylab(expression("Conc. [g DO"*~m^{-3}*"]")) +
    scale_color_manual(values = c('red1','red4','lightblue3','lightblue1','gold')) +
    xlab('') +
    ylim(c(-2,25)) +  theme_minimal()+
    theme(legend.text = element_text(size = 11), axis.text.x= element_text(size = 20), plot.title = element_text(size = 20),
          axis.text.y= element_text(size = 20), text = element_text(size = 20), legend.title = element_blank(), strip.text =element_text(size = 20),
          legend.position = 'bottom'); plot

  return(list('df' = o2_data,
              'fit' = fit,
              'plot' = plot))
}

optim_odem_static <- function(p, input.values, nep = 1000, min = 100, sed = 3000,
                     wind, khalf = 500, elev = NULL, verbose,  startdate = NULL, enddate = NULL, field.values){

  p <- lb+(ub - lb)/(10)*(p)

  o2 <- odem_static(input.values = input.values,
                   nep = p[1],
                   min = p[2],
                   sed = p[3],
                   wind = wind,
                   khalf = p[4],
               startdate = startdate, enddate = enddate,
               field.values = obs_weigh_df, elev = elev)
  print(o2$fit)
  return(o2$fit)
}


#' check whether the flux is valid and return a valid flux
#' @param flux to be checked
#' @param pool cuurent existing o2 in the lake
#' @return a valid flux
#' @export
#'
valid<-function(flux,pool){
  if(abs(flux)>abs(pool)){
    if(flux<0){
      return (-pool)
    }else{
      return (flux)
    }

  }else{
    return (flux)
  }
}

#' preprocesses observed data and area-weighs them
#' @param obs observed data
#' @param pool input matrix of for instance thermocline depth
#' @param H depths
#' @param A areas
#' @return matched and weighted-averaged data
#' @export
#'
preprocess_obs <- function(obs, input.values, H, A){
  deps <- seq(round(max(H),4), round(min(H),4), by = -0.5)

  if (max(H) > max(deps)){
    deps <- c(max(H), deps)
  }
  if (min(H) < min(deps)){
    deps <- c(deps, min(H))
  }

  areas <- approx(round(H,4), round(A,4), round(deps,4))$y

  apprObs <- matrix(NA, nrow= length(deps), ncol =length(unique(zoo::as.Date(obs$ActivityStartDate))))
  ts.apprObs <- matrix(NA, nrow= 3, ncol =length(unique(zoo::as.Date(obs$ActivityStartDate))))
  idx <- c()
  for (jj in unique(zoo::as.Date(obs$ActivityStartDate))){

    idy =  (match(zoo::as.Date(obs$ActivityStartDate),zoo::as.Date(jj)))
    idy <- which(!is.na(idy))
    dat <- obs[idy,]


    if (sd(dat$ActivityDepthHeightMeasure.MeasureValue) == 0 | length(dat$ActivityDepthHeightMeasure.MeasureValue) <= 1 |
        length(na.omit( round(dat$ResultMeasureValue,2))) <= 1){
      next} else {

        if (max(deps) > max(round(dat$ActivityDepthHeightMeasure.MeasureValue,2))) {
          dat <- rbind(dat, data.frame('ActivityStartDate' = zoo::as.Date(jj),
                                       'ActivityDepthHeightMeasure.MeasureValue' = max(deps),
                                       'ResultMeasureValue' = dat$ResultMeasureValue[nrow(dat)]))
        }

        if (any(is.na(approx(round(dat$ActivityDepthHeightMeasure.MeasureValue,2), round(dat$ResultMeasureValue,2),
                             deps)$y))){
          intvec <- (approx(round(dat$ActivityDepthHeightMeasure.MeasureValue,2), round(dat$ResultMeasureValue,2),
                            deps)$y)
          intvec[ which(is.na(intvec))] <- intvec[( which(is.na(intvec)))+1]
        } else {
          intvec <- approx(round(dat$ActivityDepthHeightMeasure.MeasureValue,2), round(dat$ResultMeasureValue,2),
                           deps)$y
        }
        apprObs[,match(jj, unique(zoo::as.Date(obs$ActivityStartDate)))] <-intvec

        if (zoo::as.Date(jj) < min(zoo::as.Date(input.values$datetime)) | zoo::as.Date(jj) > max(zoo::as.Date(input.values$datetime))){
          next
        } else {
          idx <- append(idx,  match(zoo::as.Date(jj), zoo::as.Date(input.values$datetime)))
          idz <-  which(zoo::as.Date(jj) == zoo::as.Date(input.values$datetime))
          if (is.na(input.values$thermocline_depth[abs(idz)])){
            dz.areas <- (1*areas)/sum(areas, na.rm= TRUE)#(areas - min(areas)) / (max(areas) - min(areas))
            ts.apprObs[1, match(zoo::as.Date(jj), unique(zoo::as.Date(obs$ActivityStartDate)))] <- weighted.mean(apprObs[,match(zoo::as.Date(jj), unique(as.Date(obs$ActivityStartDate)))],
                                                                                                                 dz.areas, na.rm = TRUE)
          } else{
            z.td <- which(abs(input.values$thermocline_depth[abs(idz)] - deps) == (min(abs(input.values$thermocline_depth[abs(idz)] - deps))[1]))
            dz.hypo <- (1*areas[1:z.td])/sum(areas[1:z.td])
            dz.epi <- (1*areas[(z.td+1):length(areas)])/sum(areas[(z.td+1):length(areas)])

            ts.apprObs[2, match(jj, unique(zoo::as.Date(obs$ActivityStartDate)))] <- weighted.mean(apprObs[(z.td+1):length(areas),match(zoo::as.Date(jj), unique(as.Date(obs$ActivityStartDate)))],
                                                                                                   dz.epi, na.rm = TRUE)

            ts.apprObs[3, match(jj, unique(zoo::as.Date(obs$ActivityStartDate)))] <- weighted.mean(apprObs[1:z.td,match(jj, unique(zoo::as.Date(obs$ActivityStartDate)))],
                                                                                                   dz.hypo, na.rm = TRUE)
          }
        }
      }
  }
  check.na <- c()
  for (p in 1:ncol(ts.apprObs)){
    if (all(is.na(ts.apprObs[,p]))){
      check.na <- append(check.na, p)
    }
  }
  # idx <- match(unique(as.Date(obs$ActivityStartDate)), as.Date(input.values$datetime))
  if (is.null(check.na)){
    return(rbind(idx, ts.apprObs))
  } else {
    return(rbind(idx, ts.apprObs[,-(check.na)]))
  }

}

#' preprocesses observed data and area-weighs them w/o interpolation
#' @param obs observed data
#' @param input.values input matrix of for instance thermocline depth
#' @param H depths
#' @param A areas
#' @return matched and weighted-averaged data
#' @export
#'
weigh_obs <- function(obs, input.values, H, A){

  data_long <- obs %>% arrange(ActivityStartDate)
  data_long$Area <- approx(H, A, data_long$ActivityDepthHeightMeasure.MeasureValue)$y

  idx <- match(zoo::as.Date(data_long$ActivityStartDate), zoo::as.Date(input.values$datetime))
  data_long$Layer <- data_long$ActivityDepthHeightMeasure.MeasureValue <= input.values$thermocline_depth[idx]

  data_long$Layer[which(data_long$Layer == TRUE)] = 'EPILIMNION'
  data_long$Layer[which(data_long$Layer == FALSE)] = 'HYPOLIMNION'
  data_long$Layer[which(is.na(data_long$Layer))] = 'TOTAL'

  data_long$WeightValue <- rep(NA, nrow(data_long))
  weight_obs <- matrix(NA, nrow = 4, ncol = length(unique(zoo::as.Date(data_long$ActivityStartDate))))

  for (ii in unique(zoo::as.Date(data_long$ActivityStartDate))){
    # print(zoo::as.Date(ii))
    idx <- which(zoo::as.Date(ii) == zoo::as.Date(data_long$ActivityStartDate))
    idz <- match(zoo::as.Date(ii), zoo::as.Date(input.values$datetime))
    thdepth <- input.values$thermocline_depth[idz]
    data <- data_long[idx, ]

    weight_obs[1, match(ii, unique(zoo::as.Date(data_long$ActivityStartDate)))] <- idz

    if (all(data$Layer == 'TOTAL')){
      total_areas <- approx(H, A, seq(from = max(H), to = 0, by = -0.5))$y
      perc <- (1 * data$Area) / sum(data$Area, na.rm = TRUE)#max(total_areas)
      data_long$WeightValue[idx] <- data$ResultMeasureValue * perc
      data$WeightValue<- data$ResultMeasureValue * perc

      # print(paste(
      #   mean(data$WeightValue[which(data$Layer == 'TOTAL')])))

      weight_obs[2, match(ii, unique(zoo::as.Date(data_long$ActivityStartDate)))] <- sum(data$WeightValue[which(data$Layer == 'TOTAL')], na.rm = TRUE)
    } else {
      idy = which(data$Layer == 'EPILIMNION')
      idy <- idy[1: floor(length(idy)*0.75)]
      epi_areas <- approx(H, A, seq(from = round(thdepth,1), to = 0, by = -0.1))$y
      epi_perc <- (1 * data$Area[idy]) / sum(data$Area[idy], na.rm = TRUE)#max(epi_areas)

      idt = which(data$Layer == 'HYPOLIMNION')
      idt <- idt[length(idt): (length(idt)-floor(length(idt)*0.75))]
      idt <- rev(idt)
      hypo_areas <- approx(H, A, seq(from = max(H), to = round(thdepth,1), by = -0.1))$y
      hypo_perc <- (1 * data$Area[idt]) / sum(data$Area[idt], na.rm = TRUE)#max(hypo_areas)

      data_long$WeightValue[idx] <- rep(NA, length(idx))
      data_long$WeightValue[idx[idy]] <- data$ResultMeasureValue[idy] * epi_perc
      data_long$WeightValue[idx[idt]] <-   data$ResultMeasureValue[idt] * hypo_perc
      # data$WeightValue <- c(data$ResultMeasureValue[idy] * epi_perc,
      #                       data$ResultMeasureValue[idt] * hypo_perc)
      data$WeightValue<- rep(NA, length(idx))
      data$WeightValue[idy] <- data$ResultMeasureValue[idy] * epi_perc
      data$WeightValue[idt] <-   data$ResultMeasureValue[idt] * hypo_perc
      #
      # print(paste(
      #   mean(data$WeightValue[which(data$Layer == 'EPILIMNION')]),
      #   mean(data$WeightValue[which(data$Layer == 'HYPOLIMNION')])
      # ))
      # trapz(x = data$ActivityDepthHeightMeasure.MeasureValue[which(data$Layer == 'EPILIMNION')], y = data$WeightValue[which(data$Layer == 'EPILIMNION')])
      # trapz(x = data$ActivityDepthHeightMeasure.MeasureValue[which(data$Layer == 'HYPOLIMNION')], y = data$WeightValue[which(data$Layer == 'HYPOLIMNION')])
      weight_obs[3, match(ii, unique(zoo::as.Date(data_long$ActivityStartDate)))] <- sum(data$WeightValue[which(data$Layer == 'EPILIMNION')], na.rm = TRUE)
      weight_obs[4, match(ii, unique(zoo::as.Date(data_long$ActivityStartDate)))] <- sum(data$WeightValue[which(data$Layer == 'HYPOLIMNION')], na.rm = TRUE)

    }
  }
  return(list(data_long, weight_obs))
}

