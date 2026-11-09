# ==============================================================================
# 00_CONFIG.R
# ==============================================================================

# ---- Percorsi ----------------------------------------------------------------
DIR_DATI <- "/run/media/ale/DATI2/Sara/frana"
DIR_SPAZIALE <- "/run/media/ale/DATI2/Sara/analisi_spaziale"
dir.create(DIR_DATI,recursive=TRUE,showWarnings=FALSE)
dir.create(DIR_SPAZIALE,recursive=TRUE,showWarnings=FALSE)

FILE_SERIE <- file.path(DIR_DATI,"NDVI_Time_Series_GAMLSS.csv")
FILE_PARAMETRI_SPAZIALI <- file.path(DIR_SPAZIALE,"parametri_spaziali.csv")
FILE_MORAN <- file.path(DIR_SPAZIALE,"serie_storica_moran.csv")
FILE_ANOMALIE <- file.path(DIR_DATI,"anomalie_GAMLSS.csv")
FILE_AOI_USATA <- file.path(DIR_DATI,"AOI_utilizzata.txt")
FILE_DOWNLOAD_SIGNATURE <- file.path(DIR_DATI,"download_signature.txt")

# ---- AOI ---------------------------------------------------------------------
AOI_WKT <- "POLYGON((137.090240 37.434309,137.094760 37.434309,137.094760 37.437913,137.090240 37.437913,137.090240 37.434309))"
CRS_AOI <- 4326
CRS_METRICO <- 32653

# AOI alternative
# AOI_WKT <- "POLYGON((136.906712 37.350276,136.917997 37.350093,136.917767 37.341084,136.906484 37.341267,136.906712 37.350276))"
# AOI_WKT <- "POLYGON((136.933474 37.339597,136.944757 37.339412,136.944525 37.330403,136.933243 37.330588,136.933474 37.339597))"
# AOI_WKT <- "POLYGON((137.041133 37.408221,137.049038 37.408084,137.048866 37.401779,137.040962 37.401915,137.041133 37.408221))"
# AOI_WKT <- "POLYGON((137.090240 37.434309,137.094760 37.434309,137.094760 37.437913,137.090240 37.437913,137.090240 37.434309))"
# AOI_WKT <- "POLYGON((137.122740 37.442198,137.127260 37.442198,137.127260 37.445802,137.122740 37.445802,137.122740 37.442198))"
# AOI_WKT <- "POLYGON((136.911118 37.341960,136.914504 37.341960,136.914504 37.344663,136.911118 37.344663,136.911118 37.341960))"
# AOI_WKT <- "POLYGON((136.937307 37.333648,136.940693 37.333648,136.940693 37.336352,136.937307 37.336352,136.937307 37.333648))"

# ---- Periodo -----------------------------------------------------------------
# Periodo usato anche dal GAMLSS
DATA_INIZIO <- as.Date("2017-01-01")
DATA_FINE <- as.Date("2025-12-31")
DATA_EVENTO <- as.Date("2024-01-01")
ORIGINE_TEMPO <- DATA_INIZIO

# Finestra usata solo dal modulo spaziale supplementare
DATA_VARIOGRAMMA_DA <- as.Date("2018-06-01")
DATA_VARIOGRAMMA_A <- as.Date("2018-08-31")

# ---- Sentinel-2 / qualità ----------------------------------------------------
MAX_TILE_CLOUD <- 30
RISOLUZIONE_M <- 10L
SCL_VALIDE <- c(4,5) # 4 = vegetazione; 5 = non vegetato/suolo nudo

# Percentuale minima di pixel NDVI validi richiesta per mantenere una scena
PROP_MIN_RASTER <- 0.70
if(!is.finite(PROP_MIN_RASTER)||PROP_MIN_RASTER<=0||PROP_MIN_RASTER>1){
  stop("PROP_MIN_RASTER deve essere compreso tra 0 e 1.")
}

# ---- Spatial thinning --------------------------------------------------------
# 1 pixel ogni 5, con esclusione dei bordi
N_LAG <- 5L
ESCLUDI_BORDI <- TRUE
MIN_PIXEL_WARNING <- 30L
DISTANZA_CAMPIONAMENTO_M <- N_LAG*RISOLUZIONE_M
if(!is.finite(N_LAG)||N_LAG<1L) stop("N_LAG deve essere almeno 1 pixel.")

# ---- Analisi spaziale supplementare ------------------------------------------
MAX_PUNTI_VARIOGRAMMA <- 1500L
SEED <- 123L
N_PERM_KS <- 4999L
ESEGUI_ANALISI_SPAZIALE <- TRUE

# ---- GAMLSS / anomalie -------------------------------------------------------
# Soglia: coda inferiore all'1%
ALPHA_CODA <- 0.01
# Primi 3 casi post-evento
N_CASI_POST <- 3L

# ---- Funzioni condivise ------------------------------------------------------
crea_aoi <- function(){
  sf::st_as_sfc(AOI_WKT,crs=CRS_AOI)
}

crea_maschera_scl <- function(scl_raster){
  m <- scl_raster==SCL_VALIDE[1]
  if(length(SCL_VALIDE)>1){
    for(v in SCL_VALIDE[-1]) m <- m|(scl_raster==v)
  }
  m
}

estrai_data_da_nome <- function(filename){
  date_str <- stringr::str_extract(basename(filename),"[0-9]{4}-[0-9]{2}-[0-9]{2}")
  as.Date(date_str,format="%Y-%m-%d")
}

get_oauth_client <- function(){
  mio_id <- "sh-65f52c32-5103-44be-b39c-290b993f671c"
  mio_secret <- "EzI11ViLLmcg7N5a2L3dllr0boLzWCOy"
  OAuthClient <- CDSE::GetOAuthClient(id=mio_id,secret=mio_secret)
  return(OAuthClient)
}

firma_download <- function(){
  paste(trimws(AOI_WKT),as.character(DATA_INIZIO),as.character(DATA_FINE),MAX_TILE_CLOUD,
        paste(SCL_VALIDE,collapse=","),PROP_MIN_RASTER,RISOLUZIONE_M,sep="|")
}
