# ==============================================================================
# 99_RUN_PIPELINE.R
# ==============================================================================

trova_cartella_script <- function(){
  args <- commandArgs(trailingOnly=FALSE)
  file_arg <- grep("^--file=",args,value=TRUE)
  if(length(file_arg)>0){
    percorso <- sub("^--file=","",file_arg[1])
    return(dirname(normalizePath(percorso,mustWork=TRUE)))
  }
  frames <- sys.frames()
  for(i in rev(seq_along(frames))){
    if(!is.null(frames[[i]]$ofile)){
      percorso <- frames[[i]]$ofile
      if(file.exists(percorso)) return(dirname(normalizePath(percorso,mustWork=TRUE)))
    }
  }
  getwd()
}

CARTELLA_PIPELINE <- trova_cartella_script()
setwd(CARTELLA_PIPELINE)
cat("\n============================================================\n")
cat("PIPELINE NDVI - GAMLSS\n")
cat("============================================================\n")
cat("Cartella pipeline:",getwd(),"\n")

script_core <- c("00_config.R","01_download.R","03_serie_temporale.R","04_gamlss.R")
mancanti <- script_core[!file.exists(script_core)]
if(length(mancanti)>0) stop("Script core mancanti: ",paste(mancanti,collapse=", "))

esegui_script <- function(nome,obbligatorio=TRUE){
  cat("\n============================================================\n")
  cat("ESECUZIONE:",nome,"\n")
  cat("============================================================\n")
  t0 <- Sys.time()
  ok <- tryCatch({source(nome,local=FALSE); TRUE},error=function(e){
    cat("ERRORE:",conditionMessage(e),"\n")
    if(obbligatorio) stop("Pipeline interrotta in ",nome,".",call.=FALSE)
    FALSE
  })
  if(ok) cat(sprintf("Completato: %s | %.1f secondi\n",nome,as.numeric(difftime(Sys.time(),t0,units="secs"))))
  invisible(ok)
}

tempo_pipeline <- Sys.time()
esegui_script("00_config.R")
esegui_script("01_download.R")

# 02 non determina N_LAG e non deve bloccare il GAMLSS.
if(isTRUE(ESEGUI_ANALISI_SPAZIALE)&&file.exists("02_analisi_spaziale.R")){
  esegui_script("02_analisi_spaziale.R",obbligatorio=FALSE)
}else cat("\n02_analisi_spaziale.R non eseguito.\n")

esegui_script("03_serie_temporale.R")
esegui_script("04_gamlss.R")