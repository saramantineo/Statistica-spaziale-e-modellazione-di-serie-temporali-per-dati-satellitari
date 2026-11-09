# ==============================================================================
# 03_SERIE_TEMPORALE.R
# ==============================================================================

library(terra)
library(lubridate)
library(dplyr)
library(stringr)
source("00_config.R")

cat("\n--- COSTRUZIONE SERIE TEMPORALE ---\n")
cat("N_LAG:",N_LAG,"pixel | distanza nominale:",DISTANZA_CAMPIONAMENTO_M,"m\n")
cat("Esclusione bordi:",ESCLUDI_BORDI,"\n")

lista_files <- list.files(DIR_DATI,pattern="^S2_.*\\.tif$",full.names=TRUE)
if(length(lista_files)==0) stop("Nessun raster S2 disponibile. Eseguire prima 01_download.R.")

date_files <- estrai_data_da_nome(lista_files)
keep <- !is.na(date_files)&date_files>=DATA_INIZIO&date_files<=DATA_FINE
lista_files <- lista_files[keep]
date_files <- date_files[keep]
lista_files <- lista_files[order(date_files)]

risultati <- data.frame()
riferimento_geom <- NULL

for(file_path in lista_files){
  nome_file <- basename(file_path)
  data_acquisizione <- estrai_data_da_nome(nome_file)
  if(is.na(data_acquisizione)) next

  r <- terra::rast(file_path)
  if("NDVI_Masked" %in% names(r)) ndvi_layer <- r[["NDVI_Masked"]]
  else if(terra::nlyr(r)>=6) ndvi_layer <- r[[6]]
  else{
    warning("NDVI_Masked assente in ",nome_file,": file ignorato.")
    next
  }

  # Controllo tecnico: tutti i raster devono riferirsi alla stessa griglia.
  if(is.null(riferimento_geom)) riferimento_geom <- ndvi_layer
  else if(!isTRUE(terra::compareGeom(riferimento_geom,ndvi_layer,stopOnError=FALSE))){
    warning("Geometria diversa in ",nome_file,": file ignorato.")
    next
  }

  mat_ndvi <- as.matrix(ndvi_layer,wide=TRUE)
  n_pixel_totali <- sum(is.finite(mat_ndvi))
  nrow_tot <- nrow(mat_ndvi)
  ncol_tot <- ncol(mat_ndvi)

  if(ESCLUDI_BORDI){
    if(nrow_tot<3||ncol_tot<3) next
    start_row <- 2; end_row <- nrow_tot-1
    start_col <- 2; end_col <- ncol_tot-1
  }else{
    start_row <- 1; end_row <- nrow_tot
    start_col <- 1; end_col <- ncol_tot
  }

  idx_rows <- seq(start_row,end_row,by=N_LAG)
  idx_cols <- seq(start_col,end_col,by=N_LAG)
  pixel_campionati <- as.vector(mat_ndvi[idx_rows,idx_cols,drop=FALSE])
  pixel_campionati <- pixel_campionati[is.finite(pixel_campionati)]
  n_effettivo <- length(pixel_campionati)

  media_camp <- if(n_effettivo>0) mean(pixel_campionati) else NA_real_
  var_camp <- if(n_effettivo>1) var(pixel_campionati) else NA_real_
  sd_camp <- if(n_effettivo>1) sd(pixel_campionati) else NA_real_
  se_camp <- if(n_effettivo>1) sd_camp/sqrt(n_effettivo) else NA_real_

  riga <- data.frame(
    file_name=nome_file,
    data=data_acquisizione,
    anno=lubridate::year(data_acquisizione),
    mese=lubridate::month(data_acquisizione),
    doy=lubridate::yday(data_acquisizione),
    media_ndvi=media_camp,
    X_t=media_camp,
    varianza_ndvi=var_camp,
    sd_spaziale=sd_camp,
    se_spaziale=se_camp,
    n_pixel_totali=n_pixel_totali,
    n_pixel_effettivi=n_effettivo,
    t_num=as.numeric(data_acquisizione-ORIGINE_TEMPO)
  )
  risultati <- bind_rows(risultati,riga)
  cat(sprintf("Processato: %s | Data: %s | n_eff: %d | Media: %.3f\n",
              nome_file,data_acquisizione,n_effettivo,media_camp))
}

risultati <- risultati %>% arrange(data) %>% distinct(data,.keep_all=TRUE)
if(nrow(risultati)==0) stop("La serie temporale è vuota.")
write.csv(risultati,FILE_SERIE,row.names=FALSE)

cat("\nEstrazione completata. File salvato in:",FILE_SERIE,"\n")
cat("Date totali:",nrow(risultati),"\n")
cat("Date con media NDVI disponibile:",sum(is.finite(risultati$media_ndvi)),"\n")
cat("Date con n_effettivo <",MIN_PIXEL_WARNING,":",
    sum(risultati$n_pixel_effettivi<MIN_PIXEL_WARNING,na.rm=TRUE),"\n")
