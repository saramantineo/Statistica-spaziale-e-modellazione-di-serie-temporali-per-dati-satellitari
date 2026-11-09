# ==============================================================================
# 02_ANALISI_SPAZIALE.R
# ==============================================================================

library(sf)
library(terra)
library(gstat)
library(ggplot2)
library(dplyr)
library(trend)
library(stringr)

source("00_config.R")

# ==============================================================================
# 1. FILE DISPONIBILI
# ==============================================================================

lista_files <- list.files(DIR_DATI, pattern = "^S2_.*\\.tif$", full.names = TRUE)

if (length(lista_files) == 0) {
  stop("Nessun raster S2 trovato. Eseguire prima 01_download.R.")
}

tab_files <- data.frame(
  file = lista_files,
  data = estrai_data_da_nome(lista_files),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(data), data >= DATA_INIZIO, data <= DATA_FINE) %>%
  arrange(data)

if (nrow(tab_files) == 0) {
  stop("Nessun raster valido nell'intervallo temporale specificato.")
}

# ==============================================================================
# 2. SCENA PRE-EVENTO PER SEMIVARIOGRAMMA
# ==============================================================================

candidati <- tab_files %>%
  filter(data >= DATA_VARIOGRAMMA_DA, data <= DATA_VARIOGRAMMA_A)

if (nrow(candidati) == 0) {
  stop("Nessun raster disponibile nella finestra scelta per il semivariogramma.")
}

conta_validi <- function(file_path) {
  r <- terra::rast(file_path)
  
  if (!("NDVI_Masked" %in% names(r)) || !("SCL" %in% names(r))) {
    return(NA_real_)
  }
  
  ndvi <- r[["NDVI_Masked"]]
  veg_valida <- !is.na(ndvi) & r[["SCL"]] == 4
  
  terra::global(veg_valida, "sum", na.rm = TRUE)[1, 1]
}

candidati$n_validi <- sapply(candidati$file, conta_validi)

candidati <- candidati %>%
  filter(!is.na(n_validi)) %>%
  arrange(desc(n_validi), data)

if (nrow(candidati) == 0) {
  stop("Nessun raster valido per il semivariogramma.")
}

file_scelto <- candidati$file[1]
data_comune <- candidati$data[1]

cat("\n--- SCENA USATA PER IL SEMIVARIOGRAMMA ---\n")
cat("Data:", as.character(data_comune), "\n")
cat("Pixel SCL 4 validi:", candidati$n_validi[1], "\n")

# ==============================================================================
# 3. RASTER PER ANALISI SPAZIALE: SOLO VEGETAZIONE SCL 4
# ==============================================================================

raster_dati <- terra::rast(file_scelto)

if (!("NDVI_Masked" %in% names(raster_dati)) || !("SCL" %in% names(raster_dati))) {
  stop("Il raster selezionato non contiene NDVI_Masked oppure SCL.")
}

ndvi_base <- raster_dati[["NDVI_Masked"]]
maschera_vegetazione <- raster_dati[["SCL"]] == 4

ndvi_pulito <- terra::mask(
  ndvi_base,
  maschera_vegetazione,
  maskvalues = FALSE
)

names(ndvi_pulito) <- "NDVI"

n_pixel_vegetazione <- terra::global(
  !is.na(ndvi_pulito),
  fun = "sum",
  na.rm = TRUE
)[1, 1]

cat("Pixel SCL 4 disponibili per semivariogramma e K-S:",
    n_pixel_vegetazione, "\n")

if (is.na(n_pixel_vegetazione) || n_pixel_vegetazione < 30) {
  stop("Troppi pochi pixel SCL 4 per eseguire le analisi spaziali.")
}

if (is.na(terra::crs(ndvi_pulito)) || terra::crs(ndvi_pulito) == "") {
  stop("Il raster non contiene un CRS valido.")
}

# ==============================================================================
# 4. CAMPIONAMENTO SPAZIALE
# ==============================================================================

cat("\n--- CAMPIONAMENTO SPAZIALE ---\n")

ndvi_utm <- terra::project(
  ndvi_pulito,
  paste0("EPSG:", CRS_METRICO),
  method = "bilinear"
)

aoi_utm <- sf::st_transform(
  crea_aoi(),
  CRS_METRICO
)

set.seed(SEED)

punti_geo <- sf::st_sample(
  aoi_utm,
  size = MAX_PUNTI_VARIOGRAMMA,
  type = "random"
)

punti_patch <- sf::st_sf(
  ID = seq_along(punti_geo),
  geometry = punti_geo
)

punti_patch$NDVI <- terra::extract(
  ndvi_utm,
  terra::vect(punti_patch),
  ID = FALSE
)[[1]]

punti_patch <- punti_patch[
  is.finite(punti_patch$NDVI),
]

if (nrow(punti_patch) < 30) {
  stop("Troppi pochi punti SCL 4 validi per stimare il semivariogramma.")
}

cat("Punti utilizzati per il semivariogramma:", nrow(punti_patch), "\n")

# ==============================================================================
# 5. CONTROLLO DEL TREND SPAZIALE
# ==============================================================================

coordinate <- sf::st_coordinates(punti_patch)

punti_patch$X_coord <- coordinate[, 1]
punti_patch$Y_coord <- coordinate[, 2]

modello_trend_spaziale <- lm(
  NDVI ~ X_coord + Y_coord,
  data = sf::st_drop_geometry(punti_patch)
)

cat("\n--- CONTROLLO TREND SPAZIALE ---\n")
print(summary(modello_trend_spaziale))

punti_patch$residuo_spaziale <- residuals(modello_trend_spaziale)
R2_trend <- summary(modello_trend_spaziale)$r.squared

cat("R^2 trend spaziale:", round(R2_trend, 4), "\n")

# ==============================================================================
# 6. SEMIVARIOGRAMMA NDVI
# ==============================================================================

cat("\n--- CALCOLO SEMIVARIOGRAMMA E RANGE ---\n")

cutoff_vario <- 150
width_vario <- RISOLUZIONE_M

v_empirico <- gstat::variogram(
  NDVI ~ 1,
  data = punti_patch,
  cutoff = cutoff_vario,
  width = width_vario
)

if (nrow(v_empirico) < 5) {
  stop("Troppi pochi bin nel semivariogramma empirico.")
}

max_g <- max(v_empirico$gamma, na.rm = TRUE)
max_d <- max(v_empirico$dist, na.rm = TRUE)

if (!is.finite(max_g) || max_g <= 0) {
  stop("Semivarianza empirica non valida.")
}

if (!is.finite(max_d) || max_d <= 0) {
  stop("Distanza massima del semivariogramma non valida.")
}

modello_iniziale <- gstat::vgm(
  psill = max_g * 0.7,
  model = "Sph",
  range = max_d * 0.5,
  nugget = max_g * 0.3
)

modello_fittato <- gstat::fit.variogram(
  v_empirico,
  modello_iniziale
)

cat("\n--- MODELLO DI SEMIVARIOGRAMMA FITTATO ---\n")
print(modello_fittato)

righe_struttura <- which(modello_fittato$model != "Nug")

if (length(righe_struttura) == 0) {
  stop("Il modello del semivariogramma non contiene una componente spaziale.")
}

range_metri <- modello_fittato$range[righe_struttura[1]]

if (!is.finite(range_metri) || range_metri <= 0) {
  stop("Range stimato non valido.")
}

if (isTRUE(attr(modello_fittato, "singular"))) {
  warning("Fit del semivariogramma singolare: interpretare il range con cautela.")
}

lag_pixel <- N_LAG

cat(sprintf("\n>>> Portata spaziale stimata dal semivariogramma: %.1f m <<<\n", range_metri))
cat(sprintf(">>> Distanza di campionamento fissata: %.1f m <<<\n", DISTANZA_CAMPIONAMENTO_M))
cat(sprintf(">>> N_LAG operativo fisso: %d pixel <<<\n\n", lag_pixel))

parametri_spaziali <- data.frame(
  data_riferimento = data_comune,
  range_m = range_metri,
  n_lag = N_LAG,
  distanza_campionamento_m = DISTANZA_CAMPIONAMENTO_M,
  cutoff_m = cutoff_vario,
  width_m = width_vario,
  n_punti = nrow(punti_patch),
  r2_trend_spaziale = R2_trend
)

write.csv(
  parametri_spaziali,
  FILE_PARAMETRI_SPAZIALI,
  row.names = FALSE
)

cat("Parametri spaziali salvati in:", FILE_PARAMETRI_SPAZIALI, "\n")

titolo_variogramma <- sprintf(
  "Semivariogramma NDVI (Patch Frana)\nRange: %.1f m",
  range_metri
)

print(
  plot(
    v_empirico,
    model = modello_fittato,
    main = titolo_variogramma
  )
)

# ==============================================================================
# 7. K-S INTERNO: EST vs OVEST
# ==============================================================================

cat("\n--- K-S EST vs OVEST ---\n")

bbox_aoi <- sf::st_bbox(aoi_utm)
x_centro <- mean(c(bbox_aoi[["xmin"]], bbox_aoi[["xmax"]]))

coordinate_punti <- sf::st_coordinates(punti_patch)
punti_patch$X <- coordinate_punti[, 1]

punti_patch$Settore <- ifelse(
  punti_patch$X >= x_centro,
  "Settore EST",
  "Settore OVEST"
)

seleziona_punti_indipendenti <- function(sf_punti, min_dist, seed = 123L) {
  coordinate <- sf::st_coordinates(sf_punti)
  
  set.seed(seed)
  ordine <- sample(seq_len(nrow(sf_punti)))
  selezionati <- integer(0)
  
  for (idx in ordine) {
    if (length(selezionati) == 0) {
      selezionati <- idx
      next
    }
    
    dx <- coordinate[selezionati, 1] - coordinate[idx, 1]
    dy <- coordinate[selezionati, 2] - coordinate[idx, 2]
    distanze <- sqrt(dx^2 + dy^2)
    
    if (all(distanze >= min_dist)) {
      selezionati <- c(selezionati, idx)
    }
  }
  
  sf_punti[selezionati, ]
}

punti_indipendenti <- seleziona_punti_indipendenti(
  punti_patch,
  min_dist = DISTANZA_CAMPIONAMENTO_M,
  seed = SEED
)

cat("Punti prima del thinning:", nrow(punti_patch), "\n")
cat("Punti dopo il thinning:", nrow(punti_indipendenti), "\n")
cat(sprintf("Distanza minima utilizzata: %.1f m (%d pixel)\n",
            DISTANZA_CAMPIONAMENTO_M, N_LAG))

ndvi_est <- punti_indipendenti$NDVI[
  punti_indipendenti$Settore == "Settore EST"
]

ndvi_ovest <- punti_indipendenti$NDVI[
  punti_indipendenti$Settore == "Settore OVEST"
]

cat("N EST:", length(ndvi_est), "\n")
cat("N OVEST:", length(ndvi_ovest), "\n")

if (length(ndvi_est) < 5 || length(ndvi_ovest) < 5) {
  
  warning("Campioni troppo piccoli dopo il thinning: K-S non eseguito.")
  ks_result <- NULL
  
} else {
  
  ks_result <- suppressWarnings(
    ks.test(
      ndvi_est,
      ndvi_ovest,
      exact = FALSE
    )
  )
  
  cat("\n--- K-S CLASSICO ---\n")
  print(ks_result)
  
  calcola_D_ks <- function(x, y) {
    valori <- sort(unique(c(x, y)))
    F_x <- stats::ecdf(x)
    F_y <- stats::ecdf(y)
    
    max(abs(F_x(valori) - F_y(valori)))
  }
  
  D_osservato <- calcola_D_ks(ndvi_est, ndvi_ovest)
  
  pool_ndvi <- c(ndvi_est, ndvi_ovest)
  n_est <- length(ndvi_est)
  
  set.seed(SEED)
  
  D_perm <- replicate(
    N_PERM_KS,
    {
      indici_est <- sample(
        seq_along(pool_ndvi),
        size = n_est,
        replace = FALSE
      )
      
      calcola_D_ks(
        pool_ndvi[indici_est],
        pool_ndvi[-indici_est]
      )
    }
  )
  
  p_perm <- (
    1 + sum(D_perm >= D_osservato)
  ) / (
    N_PERM_KS + 1
  )
  
  cat("\n--- K-S PERMUTAZIONALE ---\n")
  cat(sprintf("D osservato = %.4f\n", D_osservato))
  cat(sprintf("Numero permutazioni = %d\n", N_PERM_KS))
  cat(sprintf("p-value permutazionale = %.5f\n", p_perm))
  
  if (p_perm < 0.05) {
    cat(
      "Interpretazione: esiste evidenza statistica di una differenza ",
      "tra le distribuzioni NDVI dei settori EST e OVEST.\n"
    )
  } else {
    cat(
      "Interpretazione: non emerge evidenza statistica sufficiente ",
      "per distinguere le distribuzioni NDVI dei due settori.\n"
    )
  }
  
  print(
    ggplot(
      sf::st_drop_geometry(punti_indipendenti),
      aes(x = NDVI, fill = Settore, color = Settore)
    ) +
      geom_density(alpha = 0.35, linewidth = 1) +
      labs(
        title = "Confronto distribuzione NDVI - EST vs OVEST",
        subtitle = paste0(
          "K-S permutazionale: D = ",
          round(D_osservato, 3),
          " | p-value = ",
          format.pval(p_perm, digits = 3)
        ),
        x = "NDVI",
        y = "Densità"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  )
}

# ==============================================================================
# 8. MORAN
# ==============================================================================

cat("\n--- CALCOLO MORAN ---\n")

matrice_pesi_queen <- matrix(
  c(
    1, 1, 1,
    1, 0, 1,
    1, 1, 1
  ),
  nrow = 3,
  byrow = TRUE
)

matrice_pesi_rook <- matrix(
  c(
    0, 1, 0,
    1, 0, 1,
    0, 1, 0
  ),
  nrow = 3,
  byrow = TRUE
)

serie_storica_moran <- data.frame(
  Data = as.Date(character()),
  Moran_Queen = numeric(),
  Moran_Rook = numeric()
)

for (file_path in tab_files$file) {
  
  data_acq <- estrai_data_da_nome(file_path)
  r <- terra::rast(file_path)
  
  if (!("NDVI_Masked" %in% names(r))) next
  
  ndvi <- r[["NDVI_Masked"]]
  
  valore_queen <- tryCatch(
    terra::autocor(
      ndvi,
      w = matrice_pesi_queen,
      method = "moran",
      global = TRUE,
      standardize = TRUE
    ),
    error = function(e) NA_real_
  )
  
  valore_rook <- tryCatch(
    terra::autocor(
      ndvi,
      w = matrice_pesi_rook,
      method = "moran",
      global = TRUE,
      standardize = TRUE
    ),
    error = function(e) NA_real_
  )
  
  serie_storica_moran <- bind_rows(
    serie_storica_moran,
    data.frame(
      Data = data_acq,
      Moran_Queen = as.numeric(valore_queen),
      Moran_Rook = as.numeric(valore_rook)
    )
  )
}

serie_storica_moran <- serie_storica_moran %>%
  filter(
    !is.na(Data),
    is.finite(Moran_Queen),
    is.finite(Moran_Rook)
  ) %>%
  arrange(Data)

write.csv(
  serie_storica_moran,
  FILE_MORAN,
  row.names = FALSE
)

cat("Serie Moran salvata in:", FILE_MORAN, "\n")
cat("Numero acquisizioni valide per Moran:",
    nrow(serie_storica_moran), "\n")

# ==============================================================================
# 9. PETTITT SULLA SERIE DI MORAN
# ==============================================================================

if (nrow(serie_storica_moran) >= 10) {
  
  test_pettitt_queen <- trend::pettitt.test(
    serie_storica_moran$Moran_Queen
  )
  
  cat("\n--- TEST DI PETTITT SU MORAN (QUEEN) ---\n")
  print(test_pettitt_queen)
  
  idx_cambio_q <- as.integer(test_pettitt_queen$estimate)
  
  if (
    is.finite(idx_cambio_q) &&
    idx_cambio_q >= 1 &&
    idx_cambio_q < nrow(serie_storica_moran)
  ) {
    
    data_cambio_q <- serie_storica_moran$Data[idx_cambio_q]
    
    media_pre_q <- mean(
      serie_storica_moran$Moran_Queen[1:idx_cambio_q],
      na.rm = TRUE
    )
    
    media_post_q <- mean(
      serie_storica_moran$Moran_Queen[
        (idx_cambio_q + 1):nrow(serie_storica_moran)
      ],
      na.rm = TRUE
    )
    
    p_value_pettitt <- test_pettitt_queen$p.value
    
    par(mar = c(6, 5, 4, 2) + 0.1)
    
    plot(
      serie_storica_moran$Data,
      serie_storica_moran$Moran_Queen,
      type = "n",
      xlab = "Data di acquisizione",
      ylab = "Indice di Moran",
      main = "Test di Pettitt",
      ylim = range(
        serie_storica_moran$Moran_Queen,
        na.rm = TRUE
      )
    )
    
    grid(col = "grey85", lty = "dotted")
    
    lines(
      serie_storica_moran$Data,
      serie_storica_moran$Moran_Queen,
      col = "#2C7FB8",
      lwd = 1.8
    )
    
    points(
      serie_storica_moran$Data,
      serie_storica_moran$Moran_Queen,
      col = "#2C7FB8",
      pch = 19,
      cex = 1
    )
    
    abline(
      v = data_cambio_q,
      col = "#D73027",
      lwd = 2,
      lty = 2
    )
    
    lines(
      serie_storica_moran$Data[1:idx_cambio_q],
      rep(media_pre_q, idx_cambio_q),
      col = "#1A9850",
      lwd = 3
    )
    
    lines(
      serie_storica_moran$Data[
        (idx_cambio_q + 1):nrow(serie_storica_moran)
      ],
      rep(
        media_post_q,
        nrow(serie_storica_moran) - idx_cambio_q
      ),
      col = "#FF7F00",
      lwd = 3
    )
    
    legend(
      "bottomleft",
      legend = c(
        "Indice di Moran",
        paste0("Punto di rottura (", data_cambio_q, ")"),
        paste0("Media PRE (I = ", round(media_pre_q, 3), ")"),
        paste0("Media POST (I = ", round(media_post_q, 3), ")")
      ),
      col = c(
        "#2C7FB8",
        "#D73027",
        "#1A9850",
        "#FF7F00"
      ),
      lty = c(1, 2, 1, 1),
      lwd = c(2, 2, 3, 3),
      pch = c(19, NA, NA, NA),
      bty = "n",
      cex = 0.9
    )
    
    mtext(
      paste0(
        "Rottura stimata il: ",
        data_cambio_q,
        " | p-value (approssimazione) = ",
        format(
          p_value_pettitt,
          scientific = TRUE,
          digits = 2
        )
      ),
      side = 1,
      line = 4.5,
      cex = 0.9
    )
    
    cat(sprintf(
      "\nPettitt: change-point %s | p-value = %.6g\n",
      data_cambio_q,
      p_value_pettitt
    ))
    
    cat(sprintf("Media Moran PRE: %.3f\n", media_pre_q))
    cat(sprintf("Media Moran POST: %.3f\n", media_post_q))
    cat(sprintf(
      "Differenza POST - PRE: %.3f\n",
      media_post_q - media_pre_q
    ))
    
  } else {
    warning("Change-point di Pettitt non valido.")
  }
  
} else {
  warning("Meno di 10 osservazioni valide: Pettitt non eseguito.")
}

# ==============================================================================
# 10. RIEPILOGO FINALE
# ==============================================================================

cat("\n============================================================\n")
cat("RIEPILOGO ANALISI SPAZIALE\n")
cat("============================================================\n")
cat("Scena semivariogramma:", as.character(data_comune), "\n")
cat(sprintf("Range NDVI diagnostico: %.1f m\n", range_metri))
cat(sprintf(
  "N_LAG operativo fisso: %d pixel (%.1f m)\n",
  N_LAG,
  DISTANZA_CAMPIONAMENTO_M
))
cat(sprintf("R^2 trend spaziale: %.4f\n", R2_trend))
cat("Punti semivariogramma:", nrow(punti_patch), "\n")
cat("Punti dopo thinning:", nrow(punti_indipendenti), "\n")

if (exists("D_osservato") && exists("p_perm")) {
  cat(sprintf(
    "K-S EST vs OVEST: D = %.3f | p-permutazionale = %.5f\n",
    D_osservato,
    p_perm
  ))
}

cat("Numero acquisizioni Moran:", nrow(serie_storica_moran), "\n")

if (exists("data_cambio_q")) {
  cat("Change-point Pettitt:", as.character(data_cambio_q), "\n")
  cat(sprintf(
    "Moran PRE = %.3f | POST = %.3f\n",
    media_pre_q,
    media_post_q
  ))
}

cat("============================================================\n")