# ==============================================================================
# 04_GAMLSS.R
# ==============================================================================

library(gamlss)
library(tseries)
library(forecast)
library(dplyr)
library(ggplot2)
library(lubridate)
source("00_config.R")

# ==============================================================================
# 1. CARICAMENTO E FILTRO TEMPORALE
# ==============================================================================
if(!file.exists(FILE_SERIE)) stop("Serie temporale non trovata. Eseguire prima 03_serie_temporale.R.")
dati <- read.csv(FILE_SERIE)
dati$data <- as.Date(dati$data)

# Compatibilità con eventuali CSV precedenti
if(!"media_ndvi" %in% names(dati)){
  if("X_t" %in% names(dati)) dati$media_ndvi <- dati$X_t
  else stop("Nel CSV non esiste né media_ndvi né X_t.")
}

dati <- dati %>%
  filter(!is.na(data),data>=DATA_INIZIO,data<=DATA_FINE,is.finite(media_ndvi)) %>%
  arrange(data) %>%
  distinct(data,.keep_all=TRUE) %>%
  mutate(doy=lubridate::yday(data))

if(nrow(dati)<10) stop("Pochi dati validi nel periodo di analisi.")

# ==============================================================================
# 2. TRAIN PRE-EVENTO E PRIMI 3 CASI POST-EVENTO
# ==============================================================================
dati_train <- dati %>% filter(data<DATA_EVENTO)
dati_post <- dati %>% filter(data>=DATA_EVENTO)
if(nrow(dati_train)<30) stop("Meno di 30 osservazioni pre-evento: training insufficiente.")
if(nrow(dati_post)<N_CASI_POST) stop("Sono disponibili meno di ",N_CASI_POST," osservazioni post-evento.")
dati_test <- dati_post %>% slice_head(n=N_CASI_POST)

cat("\n--- SPLIT TEMPORALE ---\n")
cat("Periodo analisi:",as.character(DATA_INIZIO),"-",as.character(DATA_FINE),"\n")
cat("Training pre-evento:",nrow(dati_train),"osservazioni\n")
cat("Casi post-evento:",nrow(dati_test),"\n")
cat("Data evento:",as.character(DATA_EVENTO),"\n")
print(dati_test$data)

# ==============================================================================
# 3. ESPLORAZIONE TEMPORALE DEL TRAINING
# ==============================================================================
lag_grezzo <- min(100,nrow(dati_train)-1)
par(mfrow=c(1,2))
forecast::Acf(dati_train$media_ndvi,main="ACF - NDVI Grezzo (Train)",lag.max=lag_grezzo)
forecast::Pacf(dati_train$media_ndvi,main="PACF - NDVI Grezzo (Train)",lag.max=lag_grezzo)
par(mfrow=c(1,1))

intervalli_temporali <- as.numeric(diff(dati_train$data))
serie_regolare <- length(intervalli_temporali)>0&&length(unique(intervalli_temporali))==1
if(serie_regolare){
  cat("\n--- TEST ADF ---\n")
  print(tseries::adf.test(dati_train$media_ndvi,alternative="stationary"))
}else cat("\nADF non eseguito: le acquisizioni non sono temporalmente equidistanti.\n")

# ==============================================================================
# 4. MODELLAZIONE GAMLSS
# ==============================================================================
mod_gamlss <- gamlss(
  formula=media_ndvi~pbc(doy),
  sigma.formula=~pbc(doy),
  family=NO,
  data=dati_train,
  control=gamlss.control(n.cyc=100,trace=FALSE)
)
if(!isTRUE(mod_gamlss$converged)) warning("Il GAMLSS non risulta convergente.")
cat("\n--- SUMMARY GAMLSS ---\n")
summary(mod_gamlss)

# ==============================================================================
# 5. DIAGNOSTICA RESIDUI
# ==============================================================================
residui_modello <- resid(mod_gamlss)
lag_residui <- min(50,length(residui_modello)-1)
par(mfrow=c(1,2))
forecast::Acf(residui_modello,main="ACF - Residui GAMLSS",lag.max=lag_residui)
forecast::Pacf(residui_modello,main="PACF - Residui GAMLSS",lag.max=lag_residui)
par(mfrow=c(1,1))
lag_ljung <- min(10,floor(length(residui_modello)/5))
if(lag_ljung>=1){
  cat("\n--- LJUNG-BOX RESIDUI ---\n")
  print(Box.test(residui_modello,lag=lag_ljung,type="Ljung-Box",fitdf=0))
}
wp(mod_gamlss)

# ==============================================================================
# 6. FIT IN-SAMPLE E IC95% DELLA MEDIA
# ==============================================================================
dati_train$mu_fit <- as.numeric(fitted(mod_gamlss,"mu"))
dati_train$sigma_fit <- as.numeric(fitted(mod_gamlss,"sigma"))
dati_train$tipo_dato <- "Storico (Fit)"

pred_mu_train <- predict(mod_gamlss,what="mu",type="response",se.fit=TRUE)
if(is.list(pred_mu_train)&&!is.null(pred_mu_train$se.fit)){
  dati_train$se_mu <- as.numeric(pred_mu_train$se.fit)
  z95 <- qnorm(0.975)
  dati_train$IC95_lower <- dati_train$mu_fit-z95*dati_train$se_mu
  dati_train$IC95_upper <- dati_train$mu_fit+z95*dati_train$se_mu
}else{
  warning("IC95% della media non disponibile.")
  dati_train$IC95_lower <- NA_real_
  dati_train$IC95_upper <- NA_real_
}

# ==============================================================================
# 7. PREVISIONE OUT-OF-SAMPLE
# ==============================================================================
previsioni_test <- predictAll(mod_gamlss,newdata=dati_test,type="response",data=dati_train)
dati_test$mu_fit <- as.numeric(previsioni_test$mu)
dati_test$sigma_fit <- as.numeric(previsioni_test$sigma)
dati_test$tipo_dato <- "Nuova Previsione (Test)"

if(any(!is.finite(dati_test$mu_fit))||any(!is.finite(dati_test$sigma_fit))||any(dati_test$sigma_fit<=0)){
  stop("Previsioni post-evento non valide.")
}

# ==============================================================================
# 8. SOGLIA DI ALLARME: CODA INFERIORE 1%
# ==============================================================================
if(!is.finite(ALPHA_CODA)||ALPHA_CODA<=0||ALPHA_CODA>=0.5) stop("ALPHA_CODA non valido.")
z_lower <- qnorm(ALPHA_CODA)
z_upper <- qnorm(1-ALPHA_CODA)

for(nome in c("dati_train","dati_test")){
  x <- get(nome)
  x$lower_bound <- x$mu_fit+z_lower*x$sigma_fit
  x$upper_bound <- x$mu_fit+z_upper*x$sigma_fit
  x$p_lower <- pnorm(x$media_ndvi,mean=x$mu_fit,sd=x$sigma_fit)
  x$frana_flag <- x$media_ndvi<x$lower_bound
  assign(nome,x)
}

cat("\n--- RISULTATI OUT-OF-SAMPLE ---\n")
for(i in seq_len(nrow(dati_test))){
  cat("\nCASO POST-EVENTO",i,"\n")
  cat("Data:",as.character(dati_test$data[i]),"\n")
  cat(sprintf("NDVI osservato: %.4f\n",dati_test$media_ndvi[i]))
  cat(sprintf("Valore atteso (mu): %.4f\n",dati_test$mu_fit[i]))
  cat(sprintf("Sigma prevista: %.4f\n",dati_test$sigma_fit[i]))
  cat(sprintf("Soglia inferiore 1%%: %.4f\n",dati_test$lower_bound[i]))
  cat(sprintf("p-value coda inferiore: %.5f\n",dati_test$p_lower[i]))
  cat("FRANA RILEVATA:",ifelse(dati_test$frana_flag[i],"SI","NO"),"\n")
}

# ==============================================================================
# 9. PREPARAZIONE DATI PER IL GRAFICO
# ==============================================================================
dati_completi <- bind_rows(dati_train,dati_test) %>% arrange(data)

ultimo_train <- tail(dati_train,1)
primo_test <- head(dati_test,1)

if(DATA_EVENTO>ultimo_train$data && DATA_EVENTO<primo_test$data){
  intervallo_totale <- as.numeric(primo_test$data-ultimo_train$data)
  intervallo_evento <- as.numeric(DATA_EVENTO-ultimo_train$data)
  peso_evento <- intervallo_evento/intervallo_totale
  interpola <- function(pre,post) pre+peso_evento*(post-pre)
  punto_evento <- data.frame(
    data=DATA_EVENTO,
    mu_fit=interpola(ultimo_train$mu_fit,primo_test$mu_fit),
    lower_bound=interpola(ultimo_train$lower_bound,primo_test$lower_bound),
    upper_bound=interpola(ultimo_train$upper_bound,primo_test$upper_bound)
  )
}else if(DATA_EVENTO==primo_test$data){
  punto_evento <- data.frame(
    data=DATA_EVENTO,
    mu_fit=primo_test$mu_fit,
    lower_bound=primo_test$lower_bound,
    upper_bound=primo_test$upper_bound
  )
}else{
  stop("DATA_EVENTO non compresa tra ultima osservazione pre-evento e prima post-evento.")
}

# Intervallo di previsione PRE-EVENTO: banda blu
dati_ribbon_train <- bind_rows(
  dati_train %>% select(data,mu_fit,lower_bound,upper_bound),
  punto_evento
) %>% arrange(data)

# Intervallo di previsione POST-EVENTO: banda arancione
dati_ribbon_test <- bind_rows(
  punto_evento,
  dati_test %>% select(data,mu_fit,lower_bound,upper_bound)
) %>%
  arrange(data) %>%
  distinct(data,.keep_all=TRUE)

# IC 95% della media nel training
dati_ribbon_ic <- dati_train

# Linea della media stimata
dati_linea <- bind_rows(
  dati_train %>% select(data,mu_fit),
  punto_evento %>% select(data,mu_fit),
  dati_test %>% select(data,mu_fit)
) %>%
  arrange(data) %>%
  distinct(data,.keep_all=TRUE)

# ==============================================================================
# 10. SALVATAGGIO RISULTATI
# ==============================================================================
tabella_anomalie <- dati_test %>%
  select(data,media_ndvi,mu_fit,sigma_fit,lower_bound,upper_bound,p_lower,frana_flag,
         n_pixel_effettivi,everything())
write.csv(tabella_anomalie,FILE_ANOMALIE,row.names=FALSE)
cat("\nRisultati salvati in:",FILE_ANOMALIE,"\n")

# ==============================================================================
# 11. GRAFICO FINALE
# ==============================================================================
grafico_gamlss <- ggplot() +
  geom_ribbon(
    data=dati_ribbon_train,
    aes(x=data,ymin=lower_bound,ymax=upper_bound),
    fill="blue",alpha=0.18
  ) +
  geom_ribbon(
    data=dati_ribbon_ic,
    aes(x=data,ymin=IC95_lower,ymax=IC95_upper),
    fill="blue",alpha=0.15
  ) +
  geom_ribbon(
    data=dati_ribbon_test,
    aes(x=data,ymin=lower_bound,ymax=upper_bound),
    fill="darkorange",alpha=0.4
  ) +
  geom_vline(
    xintercept=as.numeric(DATA_EVENTO),
    linetype="dashed",
    color="gray60",
    linewidth=0.8
  ) +
  geom_line(
    data=dati_linea,
    aes(x=data,y=mu_fit),
    color="blue",
    linewidth=0.8
  ) +
  geom_point(
    data=dati_completi,
    aes(x=data,y=media_ndvi,color=frana_flag,shape=tipo_dato),
    size=2.5
  ) +
  scale_color_manual(values=c("FALSE"="black","TRUE"="red")) +
  scale_shape_manual(values=c("Storico (Fit)"=16,"Nuova Previsione (Test)"=17)) +
  labs(
    title="Rilevamento Anomalie NDVI tramite GAMLSS",
    subtitle=paste0("La linea tratteggiata indica la data dell'evento: ",format(DATA_EVENTO,"%d/%m/%Y")),
    x="Data",
    y="NDVI Medio",
    color="Anomalia Rilevata",
    shape="Fase dell'Analisi"
  ) +
  theme_minimal()

print(grafico_gamlss)

cat("\n============================================================\n")
cat("RIEPILOGO GAMLSS\n")
cat("============================================================\n")
cat("Training pre-evento:",nrow(dati_train),"\n")
cat("Casi post-evento:",nrow(dati_test),"\n")
cat("Anomalie pre-evento:",sum(dati_train$frana_flag,na.rm=TRUE),"su",nrow(dati_train),"\n")
cat("Anomalie post-evento:",sum(dati_test$frana_flag,na.rm=TRUE),"su",nrow(dati_test),"\n")
cat(sprintf("Soglia frana: coda inferiore %.1f%%\n",ALPHA_CODA*100))
cat("Modello media: media_ndvi ~ pbc(doy)\n")
cat("Modello sigma: ~ pbc(doy)\n")
cat("Famiglia: Normale (NO)\n")
cat("============================================================\n")