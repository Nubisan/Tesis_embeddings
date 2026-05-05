# ==============================================================================
# HCAKC.R — Clustering con Restricciones de Cardinalidad (HCAKC)
# ==============================================================================

library(cluster)
library(proxy)
library(dplyr)

if (!requireNamespace("aricode", quietly = TRUE)) install.packages("aricode")
library(aricode)

# -----------------------------
# Registro interno de hiperparámetros (solo en memoria)
# Consultar: .hcakc_hyperparams$global / .hcakc_hyperparams$runs[["nombre"]]
# -----------------------------
.hcakc_hyperparams <- list(
  global = list(
    algorithm          = "HCAKC",
    distance_method    = "cosine",
    silhouette_metric  = "euclidean",
    constraint_matrices = "M_zero_C_zero",
    centroid_update    = "colMeans",
    init_method        = "sequential_block",
    objective          = "IS_cross_product_with_CUCMC"
  ),
  runs = list()
)

compute_metrics <- function(y_true, y_predict, X) {
  y_true_int    <- as.integer(as.factor(y_true))
  y_predict_int <- as.integer(y_predict)
  ARI_val <- tryCatch(aricode::ARI(y_true_int, y_predict_int),  error = function(e) NA_real_)
  AMI_val <- tryCatch(aricode::AMI(y_true_int, y_predict_int),  error = function(e) NA_real_)
  NMI_val <- tryCatch(aricode::NMI(y_true_int, y_predict_int),  error = function(e) NA_real_)
  Sil_val <- tryCatch({
    if (length(unique(y_predict_int)) < 2) return(NA_real_)
    d   <- proxy::dist(as.matrix(X), method = "cosine")
    sil <- cluster::silhouette(y_predict_int, d)
    mean(sil[, "sil_width"], na.rm = TRUE)
  }, error = function(e) NA_real_)
  list(ARI = ARI_val, AMI = AMI_val, NMI = NMI_val, Silhouette_mean = Sil_val)
}

# Salida
dir.create("predictions", showWarnings = FALSE, recursive = FALSE)
.pred_hcakc_csv <- "predictions/pred_HCAKC.csv"

# -----------------------------
# Helpers robustos
# -----------------------------
make_numeric_X <- function(X) {
  X <- as.data.frame(X)
  
  X_num <- lapply(X, function(col) {
    if (is.numeric(col) || is.integer(col)) return(as.numeric(col))
    if (is.logical(col)) return(as.numeric(col))
    if (is.factor(col))  return(as.numeric(col))
    if (is.character(col)) return(as.numeric(as.factor(col)))
    return(NULL)
  })
  
  keep <- !vapply(X_num, is.null, logical(1))
  X_num <- X_num[keep]
  if (length(X_num) == 0) return(NULL)
  
  X_num <- as.data.frame(X_num)
  
  X_num <- X_num[, vapply(X_num, function(z) length(unique(z)) > 1, logical(1)), drop = FALSE]
  if (ncol(X_num) == 0) return(NULL)
  
  X_num
}

prepare_data <- function(dataset) {
  dataset <- as.data.frame(dataset)
  if (ncol(dataset) < 2) return(NULL)
  
  if ("class" %in% colnames(dataset)) {
    y <- dataset$class
    X <- dataset[, setdiff(colnames(dataset), "class"), drop = FALSE]
  } else {
    y <- dataset[[ncol(dataset)]]
    X <- dataset[, -ncol(dataset), drop = FALSE]
  }
  
  X_num <- make_numeric_X(X)
  if (is.null(X_num)) return(NULL)
  
  list(X = X_num, y = as.factor(y))
}

# -----------------------------
# HCAKC
# OPT: IS vectorizada — calcula distancias de todos los puntos a todos los
#      centroides en una sola operación matricial, en vez de sapply punto por punto
# OPT: CUCMC vectorizada con indexación lógica directa en vez de doble for
# -----------------------------
run_HCAKC_algo <- function(X, y, target_cardinality) {
  if (is.null(X) || is.null(y) || is.null(target_cardinality)) stop("Datos nulos detectados.")
  X <- as.matrix(X)
  
  if (nrow(X) != length(y)) stop("nrow(X) != length(y)")
  if (any(is.na(target_cardinality))) stop("target_cardinality contiene NA")
  if (sum(target_cardinality) != nrow(X)) stop("Suma de cardinalidades != nrow(X)")
  
  K <- length(target_cardinality)
  if (K < 2) stop("K < 2")
  
  M <- matrix(0, K, K)
  C <- matrix(0, K, K)
  
  # OPT: IS vectorizada para todos los puntos de una vez
  # Calcula (b - a) / max(a, b) para cada punto
  IS_vectorized <- function(data, clusters, centroids) {
    n <- nrow(data)
    # Matriz de distancias euclidianas de cada punto a cada centroide (n x K)
    # Cada fila i: sqrt(sum((data[i,] - centroids[j,])^2)) para j=1..K
    diff_sq <- matrix(0, n, nrow(centroids))
    for (j in seq_len(nrow(centroids))) {
      diff_sq[, j] <- sqrt(rowSums((data - matrix(centroids[j, ], n, ncol(data), byrow = TRUE))^2))
    }
    # a = distancia al propio centroide
    a <- diff_sq[cbind(seq_len(n), clusters)]
    # b = distancia mínima a cualquier otro centroide
    # Poner Inf en la columna del propio cluster
    diff_sq_other <- diff_sq
    diff_sq_other[cbind(seq_len(n), clusters)] <- Inf
    b <- apply(diff_sq_other, 1, min)
    
    (b - a) / pmax(a, b)
  }
  
  # OPT: CUCMC vectorizada
  CUCMC <- function(Xm, M, C) {
    Xm[M == 1] <- 1
    Xm[C == 1 & M != 1] <- 0
    Xm
  }
  
  HCAKC_specified <- function(data, K, cluster_sizes, M, C) {
    n <- nrow(data); p <- ncol(data)
    clusters <- integer(n)
    centroids <- matrix(NA_real_, nrow = K, ncol = p)
    
    start_idx <- 1
    for (k in seq_len(K)) {
      end_idx <- start_idx + cluster_sizes[k] - 1
      if (end_idx > n) stop("Índice fuera de rango en inicialización.")
      clusters[start_idx:end_idx] <- k
      centroids[k, ] <- colMeans(data[start_idx:end_idx, , drop = FALSE])
      start_idx <- end_idx + 1
    }
    
    # OPT: vectorizada en vez de sapply punto por punto
    IS_values <- IS_vectorized(data, clusters, centroids)
    
    Xmat <- matrix(0, K, K)
    for (i in seq_len(K)) {
      for (j in seq_len(K)) {
        IS_i <- IS_values[clusters == i]
        IS_j <- IS_values[clusters == j]
        if (length(IS_i) == 0 || length(IS_j) == 0) next
        min_len <- min(length(IS_i), length(IS_j))
        Xmat[i, j] <- sum(head(IS_i, min_len) * head(IS_j, min_len), na.rm = TRUE)
      }
    }
    Xmat <- CUCMC(Xmat, M, C)
    
    list(clusters = clusters, centroids = centroids)
  }
  
  result_specified <- HCAKC_specified(X, K, as.integer(target_cardinality), M, C)
  result_specified$clusters
}

# -----------------------------
# Runner: devuelve fila
# -----------------------------
run_clustering_row <- function(dataset, target_cardinality, dataset_name) {
  data <- prepare_data(dataset)
  if (is.null(data)) return(NULL)
  
  X <- data$X
  y <- data$y
  
  if (is.null(target_cardinality) || all(is.na(target_cardinality)) || length(target_cardinality) < 2) return(NULL)
  if (sum(target_cardinality) != nrow(X)) return(NULL)
  
  start_total <- proc.time()
  
  clusters <- tryCatch(run_HCAKC_algo(X, y, as.integer(target_cardinality)), error = function(e) NULL)
  
  if (is.null(clusters)) return(NULL)
  
  y_int <- as.integer(y)
  
  total_time <- (proc.time() - start_total)[3]
  
  metrics <- compute_metrics(y, clusters, X)
  
  # Guardar hiperparámetros de esta ejecución (solo en memoria)
  .hcakc_hyperparams$runs[[dataset_name]] <<- list(
    k                   = length(target_cardinality),
    n                   = length(clusters),
    n_features          = ncol(as.matrix(X)),
    target_sizes        = as.integer(target_cardinality),
    distance_method     = "cosine",
    silhouette_metric   = "euclidean",
    constraint_matrices = "M_zero_C_zero",
    init_method         = "sequential_block"
  )
  
  data.frame(
    name = dataset_name,
    n = length(clusters),
    k = length(unique(clusters)),
    y_predict = paste(as.integer(clusters), collapse = " "),
    y_true = paste(as.integer(y_int), collapse = " "),
    target_cardinality = paste(as.integer(target_cardinality), collapse = " "),
    cardinality_pred   = paste(as.integer(table(factor(clusters,
                                                       levels = 1:length(target_cardinality)))), collapse = " "),
    Execution_Time = as.numeric(total_time),
    ARI = metrics$ARI,
    AMI = metrics$AMI,
    NMI = metrics$NMI,
    Silhouette_mean = metrics$Silhouette_mean,
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Loop principal
# OPT: pre-asignar lista, seq_len, vapply
# -----------------------------
results_list <- vector("list", nrow(odatasets_unique))

for (i in seq_len(nrow(odatasets_unique))) {
  cat("\n\n--- Executing HCAKC for dataset at position:", i, "---\n")
  
  tryCatch({
    dataset <- odatasets_unique$dataset[[i]]
    dataset_name <- odatasets_unique$name[[i]]
    target_cardinality <- odatasets_unique$class_distribution_vector[[i]]
    
    if (is.null(dataset) || is.null(target_cardinality) || all(is.na(target_cardinality))) {
      cat("Dataset o cardinalidad inválida. Skipping.\n")
      next
    }
    
    row_result <- run_clustering_row(dataset, target_cardinality, dataset_name)
    
    if (!is.null(row_result)) {
      results_list[[i]] <- row_result
      cat("Time:", row_result$Execution_Time, "seconds\n")
    } else {
      cat("Dataset skipped.\n")
    }
    
  }, error = function(e) {
    cat("Error processing dataset at position", i, ":", e$message, "\n")
  })
}

# -----------------------------
# Escritura final
# -----------------------------
results_clean <- results_list[!vapply(results_list, is.null, logical(1))]

if (length(results_clean) > 0) {
  final_df <- do.call(rbind, results_clean)
  write.table(final_df, .pred_hcakc_csv, sep = ",", row.names = FALSE, col.names = TRUE)
  cat("\nHCAKC finalizado. Archivo guardado exitosamente en:", .pred_hcakc_csv, "\n")
} else {
  cat("\nHCAKC finalizado sin resultados válidos. No se generó archivo CSV.\n")
}