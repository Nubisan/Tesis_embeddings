# -----------------------------------------------------------------------------
# Load libraries
# -----------------------------------------------------------------------------
library(cluster)
library(proxy)
library(lpSolve)
library(dplyr)

if (!requireNamespace("aricode", quietly = TRUE)) install.packages("aricode")
library(aricode)

compute_metrics <- function(y_true, y_predict, X, dist_matrix = NULL) {
  y_true_int    <- as.integer(as.factor(y_true))
  y_predict_int <- as.integer(y_predict)
  ARI_val <- tryCatch(aricode::ARI(y_true_int, y_predict_int),  error = function(e) NA_real_)
  AMI_val <- tryCatch(aricode::AMI(y_true_int, y_predict_int),  error = function(e) NA_real_)
  NMI_val <- tryCatch(aricode::NMI(y_true_int, y_predict_int),  error = function(e) NA_real_)
  Sil_val <- tryCatch({
    if (length(unique(y_predict_int)) < 2) return(NA_real_)
    # Reutilizar la matriz de distancia completa si ya fue calculada
    if (is.null(dist_matrix)) {
      dist_matrix <- proxy::dist(as.matrix(X), method = "cosine")
    }
    sil <- cluster::silhouette(y_predict_int, dist_matrix)
    mean(sil[, "sil_width"], na.rm = TRUE)
  }, error = function(e) NA_real_)
  list(ARI = ARI_val, AMI = AMI_val, NMI = NMI_val, Silhouette_mean = Sil_val)
}

# Salida
dir.create("predictions", showWarnings = FALSE, recursive = FALSE)
.pred_csclp_csv <- "predictions/pred_CSCLP.csv"

# -----------------------------------------------------------------------------
# Helpers robustos
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# CSCLP (optimizado: matrices vectorizadas, distancia cacheada)
# -----------------------------------------------------------------------------
run_CSCLP_algo <- function(dataset, target_cardinality, dist_method = "cosine", seed = 123) {
  
  prepared <- prepare_data(dataset)
  if (is.null(prepared)) stop("Dataset inválido (no se pudo preparar X/y).")
  
  X <- prepared$X
  y <- prepared$y
  
  Xmat <- as.matrix(X)
  n <- nrow(Xmat)
  k <- length(target_cardinality)
  
  if (is.na(k) || k < 2) stop("k inválido (<2).")
  if (k >= n) stop(paste0("k=", k, " inválido para n=", n, "."))
  if (any(is.na(target_cardinality))) stop("target_cardinality contiene NA.")
  if (sum(target_cardinality) != n) stop(paste0("target_cardinality suma ", sum(target_cardinality), " pero n=", n, "."))
  if (any(target_cardinality < 1)) stop("target_cardinality tiene valores < 1 (invalido).")
  
  set.seed(seed)
  pam_result <- cluster::pam(Xmat, k = k)
  c_idx <- pam_result$id.med
  v_idx <- setdiff(1:n, c_idx)
  m <- length(v_idx)
  
  # --- Calcular distancia completa (n x n) una sola vez ---
  # Se usa para: (1) cost_mat del LP, (2) silhouette en métricas
  full_dist <- proxy::dist(Xmat, method = dist_method)
  full_dist_mat <- as.matrix(full_dist)
  
  # Extraer sub-matriz de costos: medoids (filas) x no-medoids (cols)
  cost_mat <- full_dist_mat[c_idx, v_idx, drop = FALSE]
  
  # Vector de costos para LP: transponer y aplanar
  # Orden: para cada centroide j, los costos a todos los puntos i
  f <- as.vector(t(cost_mat))
  
  # --- Restricciones vectorizadas (sin bucles) ---
  n_vars <- k * m
  
  # Restricción de punto: cada punto asignado exactamente a 1 cluster
  # A_point[i, (j-1)*m + i] = 1 para j=1..k
  row_idx_point <- rep(1:m, times = k)
  col_idx_point <- rep(seq(0, (k - 1) * m, by = m), each = m) + rep(1:m, times = k)
  A_point <- matrix(0, nrow = m, ncol = n_vars)
  A_point[cbind(row_idx_point, col_idx_point)] <- 1
  b_point <- rep(1, m)
  
  # Restricción de cardinalidad: cada cluster recibe target - 1 puntos
  # A_card[j, ((j-1)*m+1):(j*m)] = 1
  row_idx_card <- rep(1:k, each = m)
  col_idx_card <- unlist(lapply(1:k, function(j) ((j - 1) * m + 1):(j * m)))
  A_card <- matrix(0, nrow = k, ncol = n_vars)
  A_card[cbind(row_idx_card, col_idx_card)] <- 1
  b_card <- as.integer(target_cardinality) - 1L
  if (any(b_card < 0)) stop("target_cardinality[j]-1 < 0, ILP inválido.")
  
  A   <- rbind(A_point, A_card)
  rhs <- c(b_point, b_card)
  dir <- rep("==", nrow(A))
  
  start_total <- proc.time()
  
  result_lp <- lp(direction = "min",
                  objective.in = f,
                  const.mat = A,
                  const.dir = dir,
                  const.rhs = rhs,
                  all.bin = TRUE)
  
  total_time <- (proc.time() - start_total)[3]
  
  if (result_lp$status != 0) stop(paste0("ILP infeasible/failed. status=", result_lp$status))
  
  x <- result_lp$solution
  assign_mat <- matrix(x, nrow = k, byrow = TRUE)
  
  y_predict <- integer(n)
  y_predict[c_idx] <- 1:k
  
  # Asignar no-medoids: vectorizado con max.col (equivalente a which.max por columna)
  # assign_mat es k x m, queremos para cada columna i el índice de fila con max
  y_predict[v_idx] <- max.col(t(assign_mat), ties.method = "first")
  
  # Guardar hiperparámetros como atributo
  hyperparams <- list(
    seed         = seed,
    dist_method  = dist_method,
    k            = k,
    n            = n,
    lp_solver    = "lpSolve::lp",
    all_bin      = TRUE
  )
  
  list(
    y_predict    = y_predict,
    y            = y,
    X            = X,
    Execution_Time = total_time,
    full_dist    = full_dist,
    hyperparams  = hyperparams
  )
}

# -----------------------------------------------------------------------------
# Runner: devuelve fila
# -----------------------------------------------------------------------------
run_clustering_row <- function(dataset, target_cardinality, dataset_name) {
  
  result <- run_CSCLP_algo(dataset, target_cardinality)
  
  y_predict  <- result$y_predict
  y_true     <- as.integer(result$y)
  X          <- result$X
  total_time <- result$Execution_Time
  
  # Reutilizar la distancia completa ya calculada para silhouette
  metrics <- compute_metrics(y_true, y_predict, X, dist_matrix = result$full_dist)
  
  data.frame(
    name = dataset_name,
    n = length(y_predict),
    k = length(unique(y_predict)),
    y_predict = paste(as.integer(y_predict), collapse = " "),
    y_true = paste(as.integer(y_true), collapse = " "),
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

# -----------------------------------------------------------------------------
# Loop principal
# -----------------------------------------------------------------------------
n_datasets <- nrow(odatasets_unique)
results_list <- vector("list", n_datasets)   # Pre-alocar lista

for (i in 1:n_datasets) {
  cat("\n\n--- Running CSCLP for dataset at position:", i, "---\n")
  
  tryCatch({
    dataset <- odatasets_unique$dataset[[i]]
    dataset_name <- odatasets_unique$name[[i]]
    target_cardinality <- odatasets_unique$class_distribution_vector[[i]]
    
    if (is.null(dataset) || is.null(target_cardinality) || all(is.na(target_cardinality)) || length(target_cardinality) < 2) {
      cat("Dataset o cardinalidad inválida. Skipping...\n")
      next
    }
    
    row_result <- run_clustering_row(dataset, target_cardinality, dataset_name)
    
    if (!is.null(row_result)) {
      results_list[[i]] <- row_result
      cat("Time (LP):", row_result$Execution_Time, "s\n")
    } else {
      cat("Dataset skipped.\n")
    }
    
  }, error = function(e) {
    cat("Error processing dataset at position", i, ":", e$message, "\n")
  })
}

# -----------------------------------------------------------------------------
# Escritura final
# -----------------------------------------------------------------------------
results_clean <- results_list[!vapply(results_list, is.null, logical(1))]

if (length(results_clean) > 0) {
  final_df <- do.call(rbind, results_clean)
  write.table(final_df, .pred_csclp_csv, sep = ",", row.names = FALSE, col.names = TRUE)
  cat("\nCSCLP finalizado. Archivo guardado exitosamente en:", .pred_csclp_csv, "\n")
} else {
  cat("\nCSCLP finalizado sin resultados válidos. No se generó archivo CSV.\n")
}