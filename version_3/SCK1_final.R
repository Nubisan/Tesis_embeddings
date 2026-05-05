# ==============================================================================
# SCK1_final_optimized.R — SCK1 (optimizado)
# ==============================================================================

library(lpSolve)
library(proxy)
library(dplyr)

if (!requireNamespace("aricode", quietly = TRUE)) install.packages("aricode")
library(aricode)

# -----------------------------
# Registro interno de hiperparámetros
# Solo en memoria, no se imprime ni se escribe al CSV
# Consultar: .sck1_hyperparams$global / .sck1_hyperparams$runs[["nombre"]]
# -----------------------------
.sck1_hyperparams <- list(
  global = list(
    distance_method    = "cosine",
    ilp_direction      = "min",
    default_seed       = 123L,
    max_iterations     = 1L,
    centroid_init      = "random_sample",
    assignment_method  = "ilp_exact_equality",
    centroid_update    = "colMeans"
  ),
  runs = list()
)

# -----------------------------
# Métricas
# Acepta dist_mat opcional para evitar recalcular la distancia NxN coseno
# Si no se pasa, la calcula internamente (comportamiento original)
# -----------------------------
compute_metrics <- function(y_true, y_predict, X, dist_mat = NULL) {
  y_true_int    <- as.integer(as.factor(y_true))
  y_predict_int <- as.integer(y_predict)
  ARI_val <- tryCatch(aricode::ARI(y_true_int, y_predict_int),  error = function(e) NA_real_)
  AMI_val <- tryCatch(aricode::AMI(y_true_int, y_predict_int),  error = function(e) NA_real_)
  NMI_val <- tryCatch(aricode::NMI(y_true_int, y_predict_int),  error = function(e) NA_real_)
  Sil_val <- tryCatch({
    if (length(unique(y_predict_int)) < 2) return(NA_real_)
    if (is.null(dist_mat)) {
      dist_mat <- proxy::dist(as.matrix(X), method = "cosine")
    }
    sil <- cluster::silhouette(y_predict_int, dist_mat)
    mean(sil[, "sil_width"], na.rm = TRUE)
  }, error = function(e) NA_real_)
  list(ARI = ARI_val, AMI = AMI_val, NMI = NMI_val, Silhouette_mean = Sil_val)
}

# Salida
dir.create("predictions", showWarnings = FALSE, recursive = FALSE)
.pred_sck1_csv <- "predictions/pred_SCK1.csv"

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
  
  # Eliminar columnas constantes
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
# ILP: asignación exacta (==)
# OPT: A1 y A2 construidas por indexación vectorizada (sin bucles for)
# Produce la misma matriz numérica que el original
# -----------------------------
solve_ilp_assignment <- function(cost_mat, sizes) {
  cost_mat <- as.matrix(cost_mat)
  k <- nrow(cost_mat)
  n <- ncol(cost_mat)
  
  sizes <- as.integer(sizes)
  if (length(sizes) != k) stop("sizes length != k")
  if (sum(sizes) != n) stop("sum(sizes) != n (inconsistente)")
  if (k < 2) stop("k < 2")
  
  cost_vec <- as.vector(t(cost_mat))
  total_vars <- k * n
  
  # A1: cada punto j asignado a exactamente 1 cluster
  # Original: for (j in 1:n) A1[j, seq(j, k*n, by=n)] <- 1
  # Las posiciones por fila j son: j, j+n, j+2n, ..., j+(k-1)*n
  A1 <- matrix(0, n, total_vars)
  a1_rows <- rep(1:n, each = k)
  a1_cols <- rep(1:n, each = k) + rep(seq(0, (k - 1) * n, by = n), times = n)
  A1[cbind(a1_rows, a1_cols)] <- 1
  
  # A2: cada cluster i tiene exactamente sizes[i] puntos
  # Original: for (i in 1:k) A2[i, ((i-1)*n+1):(i*n)] <- 1
  A2 <- matrix(0, k, total_vars)
  a2_rows <- rep(1:k, each = n)
  a2_cols <- sequence(rep(n, k), from = seq(1, (k - 1) * n + 1, by = n))
  A2[cbind(a2_rows, a2_cols)] <- 1
  
  res <- lp(
    direction = "min",
    objective.in = cost_vec,
    const.mat = rbind(A1, A2),
    const.dir = c(rep("=", n), rep("=", k)),
    const.rhs = c(rep(1, n), sizes),
    all.bin = TRUE
  )
  
  if (res$status != 0) stop("ILP failed")
  
  matrix(res$solution, nrow = k, byrow = TRUE)
}

# -----------------------------
# SCK1
# OPT: max.col(t()) en lugar de apply(, 2, which.max)
# -----------------------------
sck1_iterativo <- function(data, k, initial_centroids, sizes) {
  data <- as.matrix(data)
  centroids <- as.matrix(initial_centroids)
  
  if (nrow(data) != sum(sizes)) stop("nrow(data) != sum(sizes)")
  if (nrow(centroids) != k) stop("nrow(centroids) != k")
  
  cost_mat <- t(as.matrix(proxy::dist(data, centroids, method = "cosine")))
  assign_mat <- solve_ilp_assignment(cost_mat, sizes)
  
  # OPT: max.col opera sobre filas, transpose para obtener which.max por columna
  assignments <- max.col(t(assign_mat), ties.method = "first")
  
  for (i in 1:k) {
    idx <- which(assignments == i)
    if (length(idx) > 0L) centroids[i, ] <- colMeans(data[idx, , drop = FALSE])
  }
  
  list(assignments = assignments, centroids = centroids)
}

# -----------------------------
# Runner: devuelve fila
# OPT: distancia NxN coseno calculada una sola vez y pasada a compute_metrics
# -----------------------------
run_clustering_row <- function(dataset, target_cardinality, dataset_name, seed = 123) {
  data <- prepare_data(dataset)
  if (is.null(data)) return(NULL)
  
  X <- as.matrix(data$X)
  y <- data$y
  
  if (is.null(target_cardinality) || all(is.na(target_cardinality))) return(NULL)
  target_cardinality <- as.integer(target_cardinality)
  
  k <- length(target_cardinality)
  if (k < 2) return(NULL)
  if (sum(target_cardinality) != nrow(X)) return(NULL)
  
  start_total <- proc.time()
  
  set.seed(seed)
  initial_centroids <- X[sample(1:nrow(X), k), , drop = FALSE]
  
  res <- tryCatch(
    sck1_iterativo(X, k, initial_centroids, target_cardinality),
    error = function(e) NULL
  )
  
  if (is.null(res) || is.null(res$assignments)) return(NULL)
  
  y_predict <- res$assignments
  y_int <- as.integer(y)
  
  total_time <- (proc.time() - start_total)[3]
  
  # OPT: calcular distancia NxN una sola vez y reutilizarla
  dist_full <- proxy::dist(as.matrix(X), method = "cosine")
  metrics <- compute_metrics(y, y_predict, X, dist_mat = dist_full)
  
  # Guardar hiperparámetros de esta ejecución
  .sck1_hyperparams$runs[[dataset_name]] <<- list(
    seed              = seed,
    k                 = k,
    n                 = nrow(X),
    n_features        = ncol(X),
    target_sizes      = target_cardinality,
    distance_method   = "cosine",
    max_iterations    = 1L,
    centroid_init     = "random_sample",
    assignment_method = "ilp_exact_equality"
  )
  
  data.frame(
    name = dataset_name,
    n = length(y_predict),
    k = length(unique(y_predict)),
    y_predict = paste(as.integer(y_predict), collapse = " "),
    y_true = paste(as.integer(y_int), collapse = " "),
    target_cardinality = paste(as.integer(target_cardinality), collapse = " "),
    cardinality_pred   = paste(as.integer(table(factor(y_predict,
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
# OPT: pre-asignación de lista con tamaño conocido
# -----------------------------
results_list <- vector("list", nrow(odatasets_unique))

for (i in seq_len(nrow(odatasets_unique))) {
  cat("\n\n--- Executing SCK1 for dataset at position:", i, "---\n")
  
  tryCatch({
    dataset <- odatasets_unique$dataset[[i]]
    dataset_name <- odatasets_unique$name[[i]]
    target_cardinality <- odatasets_unique$class_distribution_vector[[i]]
    
    if (is.null(dataset) || is.null(target_cardinality) || all(is.na(target_cardinality))) {
      cat("Dataset o cardinalidad inválida. Skipping.\n")
      next
    }
    
    row_result <- run_clustering_row(dataset, target_cardinality, dataset_name, seed = 123)
    
    if (!is.null(row_result)) {
      results_list[[i]] <- row_result
      cat("Time:", row_result$Execution_Time, "seconds\n")
    } else {
      cat("Dataset skipped.\n")
    }
    
  }, error = function(e) cat("Error:", e$message, "\n"))
}

# -----------------------------
# Escritura final
# OPT: vapply con tipo explícito en lugar de sapply
# -----------------------------
results_clean <- results_list[!vapply(results_list, is.null, logical(1))]

if (length(results_clean) > 0) {
  final_df <- do.call(rbind, results_clean)
  write.table(final_df, .pred_sck1_csv, sep = ",", row.names = FALSE, col.names = TRUE)
  cat("\nSCK1 finalizado. Archivo guardado exitosamente en:", .pred_sck1_csv, "\n")
} else {
  cat("\nSCK1 finalizado sin resultados válidos. No se creó el CSV.\n")
}

# Hiperparámetros accesibles en memoria:
#   .sck1_hyperparams$global    -> parámetros globales del algoritmo
#   .sck1_hyperparams$runs      -> parámetros por dataset ejecutado
