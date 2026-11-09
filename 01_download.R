# ==============================================================================
# 01_DOWNLOAD.R
# ==============================================================================

library(CDSE)
library(sf)
library(terra)
library(dplyr)
library(stringr)
source("00_config.R")

OAuthClient <- get_oauth_client()
aoi <- crea_aoi()

# ==============================================================================
# 1. CONTROLLO CACHE
# ==============================================================================
cat("\n--- CONTROLLO CONFIGURAZIONE DOWNLOAD ---\n")
raster_esistenti <- list.files(DIR_DATI,pattern="^S2_.*\\.tif$",full.names=TRUE)
firma_corrente <- firma_download()
rigenera_tutto <- FALSE

if(file.exists(FILE_DOWNLOAD_SIGNATURE)){
  firma_precedente <- paste(readLines(FILE_DOWNLOAD_SIGNATURE,warn=FALSE),collapse="")
  if(!identical(trimws(firma_precedente),trimws(firma_corrente))){
    cat(">>> Configurazione di preprocessing modificata: rigenerazione completa.\n")
    rigenera_tutto <- TRUE
  }else cat("Configurazione invariata: riuso dei raster già presenti.\n")
}else if(length(raster_esistenti)>0){
  cat(">>> Raster presenti senza firma: rigenerazione completa.\n")
  rigenera_tutto <- TRUE
}else cat("Prima esecuzione: nessun raster precedente trovato.\n")

if(rigenera_tutto){
  cat("\n--- PULIZIA OUTPUT PRECEDENTI ---\n")
  file_raster <- list.files(DIR_DATI,pattern="^(S2_|temp_).*\\.tif$",full.names=TRUE)
  if(length(file_raster)>0) file.remove(file_raster)
  for(f in c(FILE_SERIE,FILE_ANOMALIE,FILE_PARAMETRI_SPAZIALI,FILE_MORAN)){
    if(file.exists(f)) file.remove(f)
  }
}

writeLines(AOI_WKT,FILE_AOI_USATA)
writeLines(firma_corrente,FILE_DOWNLOAD_SIGNATURE)

# ==============================================================================
# 2. RICERCA ACQUISIZIONI
# ==============================================================================
cat("\n--- RICERCA ACQUISIZIONI SENTINEL-2 ---\n")
immagini_disp <- SearchCatalog(aoi=aoi,from=DATA_INIZIO,to=DATA_FINE,
                               collection="sentinel-2-l2a",client=OAuthClient)
if(is.null(immagini_disp)||nrow(immagini_disp)==0) stop("Nessuna acquisizione Sentinel-2 trovata.")

immagini_ordinate <- immagini_disp %>%
  mutate(data_acquisizione=as.Date(acquisitionDate)) %>%
  filter(!is.na(data_acquisizione),!is.na(tileCloudCover),tileCloudCover<=MAX_TILE_CLOUD) %>%
  arrange(data_acquisizione,tileCloudCover) %>%
  distinct(data_acquisizione,.keep_all=TRUE) %>%
  arrange(data_acquisizione)

if(nrow(immagini_ordinate)==0) stop("Nessuna acquisizione supera il filtro di copertura nuvolosa.")
cat("Acquisizioni candidate:",nrow(immagini_ordinate),"\n")

# ==============================================================================
# 3. EVALSCRIPT
# ==============================================================================
mio_evalscript <- "
//VERSION=3
function setup(){
  return{input:['B04','B03','B02','B08','SCL'],output:{id:'default',bands:5,sampleType:'UINT16'}};
}
function evaluatePixel(sample){
  return[sample.B04*10000,sample.B03*10000,sample.B02*10000,sample.B08*10000,sample.SCL];
}
"

# ==============================================================================
# 4. DOWNLOAD E PREPROCESSING
# ==============================================================================
n_salvate <- 0L
n_scartate <- 0L
n_errori <- 0L

for(i in seq_len(nrow(immagini_ordinate))){
  data_acquisizione <- immagini_ordinate$data_acquisizione[i]
  file_finale <- file.path(DIR_DATI,paste0("S2_",data_acquisizione,".tif"))
  if(file.exists(file_finale)){
    cat("Già presente:",basename(file_finale),"\n")
    next
  }

  file_temp <- file.path(DIR_DATI,paste0("temp_",data_acquisizione,".tif"))
  if(file.exists(file_temp)) file.remove(file_temp)
  cat("\nDownload:",as.character(data_acquisizione),"\n")

  ok <- tryCatch({
    GetImage(aoi=aoi,time_range=data_acquisizione,collection="sentinel-2-l2a",
             script=mio_evalscript,resolution=RISOLUZIONE_M,format="image/tiff",
             client=OAuthClient,file=file_temp)
    TRUE
  },error=function(e){
    message("Errore nel download del ",data_acquisizione,": ",e$message)
    FALSE
  })

  if(!ok||!file.exists(file_temp)){
    if(file.exists(file_temp)) file.remove(file_temp)
    n_errori <- n_errori+1L
    next
  }

  raster_dati <- tryCatch(terra::rast(file_temp),error=function(e) NULL)
  if(is.null(raster_dati)||terra::nlyr(raster_dati)!=5){
    message("Raster non valido per la data ",data_acquisizione)
    if(file.exists(file_temp)) file.remove(file_temp)
    n_errori <- n_errori+1L
    next
  }

  names(raster_dati) <- c("Rosso","Verde","Blu","NIR","SCL")
  denominatore <- raster_dati[["NIR"]]+raster_dati[["Rosso"]]
  ndvi_grezzo <- terra::ifel(denominatore!=0,
                             (raster_dati[["NIR"]]-raster_dati[["Rosso"]])/denominatore,NA)
  maschera_valida <- crea_maschera_scl(raster_dati[["SCL"]])
  ndvi_pulito <- terra::mask(ndvi_grezzo,maschera_valida,maskvalues=FALSE)
  names(ndvi_pulito) <- "NDVI_Masked"

  n_pixel_totali <- terra::ncell(ndvi_pulito)
  n_validi <- terra::global(!is.na(ndvi_pulito),"sum",na.rm=TRUE)[1,1]

  if(is.na(n_pixel_totali)||!is.finite(n_pixel_totali)||n_pixel_totali<=0){
    message("Numero totale di pixel non valido per la data ",data_acquisizione)
    if(file.exists(file_temp)) file.remove(file_temp)
    n_errori <- n_errori+1L
    next
  }

  prop_validi <- n_validi/n_pixel_totali
  soglia_pixel <- ceiling(PROP_MIN_RASTER*n_pixel_totali)
  cat("Pixel validi:",n_validi,"su",n_pixel_totali,
      sprintf("(%.1f%%)",prop_validi*100),"| soglia:",soglia_pixel,
      sprintf("(%.1f%%)",PROP_MIN_RASTER*100),"\n")

  if(is.na(n_validi)||!is.finite(n_validi)||!is.finite(prop_validi)||prop_validi<PROP_MIN_RASTER){
    cat("Scartata",as.character(data_acquisizione),"- copertura valida insufficiente:",
        sprintf("%.1f%%",prop_validi*100),"\n")
    n_scartate <- n_scartate+1L
    if(file.exists(file_temp)) file.remove(file_temp)
    next
  }

  raster_completo <- c(raster_dati,ndvi_pulito)
  terra::writeRaster(raster_completo,file_finale,overwrite=TRUE)
  file.remove(file_temp)
  n_salvate <- n_salvate+1L
  cat("Salvata:",basename(file_finale),"\n")
}

cat("\n============================================================\n")
cat("RIEPILOGO DOWNLOAD\n")
cat("============================================================\n")
cat("Acquisizioni candidate:",nrow(immagini_ordinate),"\n")
cat("Nuove acquisizioni salvate:",n_salvate,"\n")
cat("Acquisizioni scartate per qualità:",n_scartate,"\n")
cat("Errori di download/lettura:",n_errori,"\n")
cat(sprintf("Percentuale minima pixel validi: %.1f%%\n",PROP_MIN_RASTER*100))
cat("Maschera SCL utilizzata:",paste(SCL_VALIDE,collapse=","),"\n")
cat("============================================================\n")
