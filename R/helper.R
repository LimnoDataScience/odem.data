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
  td.depth <- rep(NA, length(grd.info$datetime))

  condition<- apply(temp, 1, FUN=min,na.rm=TRUE) > 4

  for (ii in 1:length(cbuoy.depth)){
    idx = !is.na(temp[ii,])
    dens_data = dens[ii,idx]
    dens.diff = rev(dens_data)[1] - dens_data[1]

    if (condition[ii] && abs(dens.diff) > 0.05){
    cbuoy.depth[ii] <- center.buoyancy(temp[ii, idx], as.numeric(grd.info$depth[idx]))
    td.depth[ii] <- thermo.depth(temp[ii, idx], as.numeric(grd.info$depth[idx]))
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
  # points(td.depth, col='red')
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
#' @param td.depth matrix; thermocline depth
#' @param H the depth info of the lake
#' @return list of temperatures vector
calc_epil_hypo_temp<-function(wtemp,td.depth,H){
  grd.info <- extract_time_space(wtemp)
  temp <- as.matrix(wtemp[,-c(1)])
  depth_data = as.double(grd.info$depth)

  epil_temp <-  rep(NA, length(td.depth))
  hypo_temp <- rep(NA, length(td.depth))
  total_temp<- rep(NA, length(td.depth))

  total<- rep(NA, length(td.depth))
  hypo<- rep(NA, length(td.depth))
  total<- rep(NA, length(td.depth))

  td_not_exist <- is.na(td.depth)


  for (ii in (1:length(td.depth))){
    idx = !is.na(temp[ii,])
    temp_data = as.numeric(temp[ii,idx])
    #total_temp[ii]<-max(sum(temp_data)/length(temp_data),4)
    total_temp[ii]<-sum(temp_data)/length(temp_data)
    if(!td_not_exist[ii]){
      td_idx <- max(which(td.depth[ii]>=depth_data))
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
calc_vol_total<-function(H, A){
  if (length(H)==1){
    vol_total <- 1/3.0 * A * H
  }else{
    vol_total <- trapz(rev(H),rev(A))
  }
  return (vol_total)
}

#'
#' Calculate water epilimnion and hypolimnine volume,using the thermocline depth
#'
#' @param H vector; the depth info of the lake
#' @param A vector; the area info of the lake for each depth
#' @param td.depth matrix; thermocline depth
#' @param vol_total number;the total volume of the lake
#' @importFrom stats approx
#' @return matrix of epil and hypo volume
calc_epil_hypo_vol <- function(H,A,td.depth,vol_total){
  vol_data <- matrix(NA, nrow = length(td.depth), ncol = 2)
  colnames(vol_data) <- c("vol_epil","vol_hypo")

  td_not_exist<-is.na(td.depth)
  for (ii in 1:length(td.depth)){
    if(!td_not_exist[ii]){
      h_idx <- min(which(td.depth[ii]>=H))
      approx_td.area<-approx(H, A, c(0, td.depth[ii]))$y[2]
      if (is.na(approx_td.area)){approx_td.area = min(A)}
      H_with_td<-c(td.depth[ii],H[h_idx:length(H)])
      A_with_td<-c(approx_td.area,A[h_idx:length(A)])
      vol_data[ii,1] <- trapz(rev(H_with_td),rev(A_with_td))##epil
      vol_data[ii,2] <- vol_total - vol_data[ii,1]##hypo
      if(vol_data[ii,2] <= 0) {
        vol_data[ii,2]= min(A) * 0.5
        vol_data[ii,1] = vol_total - vol_data[ii,2]
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
  td.depth       <- calc_td_depth(wtemp)
  metalimn.depth <- calc_metalimdepth(wtemp)
  td_area        <- approx(H, A, td.depth)$y
  surf_area      <- rep(max(A), length(grd.info$datetime))
  temp_out       <- calc_epil_hypo_temp(wtemp, td.depth, H) # epi T hypo T
  vol_total      <- calc_vol_total(H, A) # epi V hypo V
  vol            <- calc_epil_hypo_vol(H, A, td.depth, vol_total)
  max.d          <- max(H)

  return(data.frame(
    datetime = as.POSIXct(grd.info$datetime),
    td.depth,
    t.epil   = temp_out$t_epil,
    t.hypo   = temp_out$t_hypo,
    t.total  = temp_out$t_total,
    vol_total, vol,
    td_area, surf_area,
    upper.metalim = metalimn.depth[1,], lower.metalim = metalimn.depth[2,],
    max.d
    ))
}
