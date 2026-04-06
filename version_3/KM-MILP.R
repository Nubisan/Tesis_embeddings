# ==============================================================================
# KM-MILP.R (optimized)
# ==============================================================================

library(lpSolve)
library(proxy)
library(dplyr)

if (!requireNamespace("aricode", quietly = TRUE)) install.packages("aricode")
library(aricode)

# -----------------------------
# Default hyperparameter values
# -----------------------------
.default_seed          <- 123L
.default_max_iter      <- 10L
.default_convergence_tol <- 1e-4
.default_dist_method   <- "cosine"

compute_metrics <- function(y_true, y_predict, dist_matrix) {
  y_true_int    <- as.integer(as.factor(y_true))
  y_predict_int <- as.integer(y_predict)
  ARI_val <- tryCatch(aricode::ARI(y_true_int, y_predict_int),  error = function(e) NA_real_)
  AMI_val <- tryCatch(aricode::AMI(y_true_int, y_predict_int),  error = function(e) NA_real_)
  NMI_val <- tryCatch(aricode::NMI(y_true_int, y_predict_int),  error = function(e) NA_real_)
  Sil_val <- tryCatch({
    if (length(unique(y_predict_int)) < 2) NA_real_
    else mean(cluster::silhouette(y_predict_int, dist_matrix)[, "sil_width"], na.rm = TRUE)
  }, error = function(e) NA_real_)
  list(ARI = ARI_val, AMI = AMI_val, NMI = NMI_val, Silhouette_mean = Sil_val)
}

# Salida
dir.create("predictions", showWarnings = FALSE, recursive = FALSE)
.pred_km_milp_csv <- "predictions/pred_KM-MILP.csv"

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
# MILP: asignación con tamaño <=
# -----------------------------
solve_milp_assignment <- function(data, centroids, size_constraints,
                                  dist_method = .default_dist_method) {
  n <- nrow(data)
  k <- nrow(centroids)
  nk <- n * k
  
  if (length(size_constraints) != k) stop("size_constraints length != k")
  if (n <= 1 || k < 2) stop("n <= 1 o k < 2")
  
  # Cost vector: distance from each point to each centroid
  cost_mat <- as.matrix(proxy::dist(data, centroids, method = dist_method))
  cost_vec <- as.vector(t(cost_mat))
  
  # Constraint 1: each point assigned to exactly one cluster (vectorized)
  row_idx1 <- rep(1:n, each = k)
  col_idx1 <- seq_len(nk)
  constr1 <- matrix(0, n, nk)
  constr1[cbind(row_idx1, col_idx1)] <- 1
  
  # Constraint 2: cluster size <= size_constraints (vectorized)
  row_idx2 <- rep(1:k, times = n)
  col_idx2 <- rep(seq(0, (n - 1) * k, by = k), each = k) + rep(1:k, times = n)
  constr2 <- matrix(0, k, nk)
  constr2[cbind(row_idx2, col_idx2)] <- 1
  
  result <- lp(
    direction    = "min",
    objective.in = cost_vec,
    const.mat    = rbind(constr1, constr2),
    const.dir    = c(rep("=", n), rep("<=", k)),
    const.rhs    = c(rep(1, n), as.numeric(size_constraints)),
    all.bin      = TRUE
  )
  
  if (result$status != 0) stop("MILP failed")
  
  sol <- matrix(result$solution, nrow = n, byrow = TRUE)
  p <- max.col(sol, ties.method = "first")
  list(p = p)
}

# -----------------------------
# Loop tipo KMeans: centroides + MILP
# -----------------------------
clustering_with_size_constraints <- function(data, size_constraints,
                                             max_iter    = .default_max_iter,
                                             seed        = .default_seed,
                                             dist_method = .default_dist_method) {
  data <- as.matrix(data)
  n <- nrow(data)
  k <- length(size_constraints)
  
  if (k < 2) stop("k < 2")
  if (sum(size_constraints) < n) stop("sum(size_constraints) < n")
  
  set.seed(seed)
  init_indices <- sample.int(n, k)
  centroids <- data[init_indices, , drop = FALSE]
  p <- rep(1L, n)
  
  converged_iter <- max_iter  # default: no early convergence
  
  for (iter in seq_len(max_iter)) {
    res <- solve_milp_assignment(data, centroids, size_constraints, dist_method)
    p_new <- res$p
    
    new_centroids <- centroids
    for (j in seq_len(k)) {
      idx <- which(p_new == j)
      if (length(idx) > 0L) new_centroids[j, ] <- colMeans(data[idx, , drop = FALSE])
    }
    
    if (max(abs(centroids - new_centroids)) < .default_convergence_tol) {
      converged_iter <- iter
      p <- p_new
      break
    }
    
    centroids <- new_centroids
    p <- p_new
  }
  
  # Guardar hiperparámetros
  hyperparams <- list(
    seed              = seed,
    dist_method       = dist_method,
    max_iter          = max_iter,
    converged_at_iter = converged_iter,
    convergence_tol   = .default_convergence_tol,
    init_indices      = init_indices,
    lp_solver         = "lpSolve::lp",
    all_bin           = TRUE,
    constraint_type   = "<= (upper bound)"
  )
  
  list(p = p, hyperparams = hyperparams)
}

# -----------------------------
# Runner: devuelve fila
# -----------------------------
run_clustering_row <- function(dataset, target_cardinality, dataset_name) {
  data <- prepare_data(dataset)
  if (is.null(data)) return(NULL)
  
  X <- data$X
  y <- data$y
  
  if (is.null(target_cardinality) || all(is.na(target_cardinality))) return(NULL)
  target_cardinality <- as.integer(target_cardinality)
  
  k <- length(target_cardinality)
  if (k < 2) return(NULL)
  
  if (sum(target_cardinality) != nrow(X)) return(NULL)
  
  start_total <- proc.time()
  
  res <- tryCatch(
    clustering_with_size_constraints(X, target_cardinality),
    error = function(e) NULL
  )
  
  if (is.null(res) || is.null(res$p)) return(NULL)
  
  y_int <- as.integer(y)
  y_predict <- res$p
  hp <- res$hyperparams
  
  total_time <- (proc.time() - start_total)[3]
  
  # Compute distance matrix once for Silhouette (reuse instead of recomputing)
  dist_matrix <- proxy::dist(as.matrix(X), method = hp$dist_method)
  
  metrics <- compute_metrics(y, y_predict, dist_matrix)
  
  row_df <- data.frame(
    name = dataset_name,
    n = length(y_predict),
    k = length(unique(y_predict)),
    y_predict = paste(as.integer(y_predict), collapse = " "),
    y_true = paste(as.integer(y_int), collapse = " "),
    target_cardinality = paste(as.integer(target_cardinality), collapse = " "),
    cardinality_pred   = paste(as.integer(table(factor(y_predict,
                                                       levels = seq_len(length(target_cardinality))))), collapse = " "),
    Execution_Time = as.numeric(total_time),
    ARI = metrics$ARI,
    AMI = metrics$AMI,
    NMI = metrics$NMI,
    Silhouette_mean = metrics$Silhouette_mean,
    stringsAsFactors = FALSE
  )
  attr(row_df, "hyperparams") <- hp
  row_df
}

# -----------------------------
# Loop principal
# -----------------------------
results_list <- vector("list", nrow(odatasets_unique))

for (i in seq_len(nrow(odatasets_unique))) {
  cat("\n\n--- Executing KM-MILP for dataset at position:", i, "---\n")
  
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
    cat("Error:", e$message, "\n")
  })
}

# -----------------------------
# Escritura final
# -----------------------------
results_clean <- results_list[!vapply(results_list, is.null, FUN.VALUE = logical(1))]

if (length(results_clean) > 0) {
  final_df <- do.call(rbind, results_clean)
  write.table(final_df, .pred_km_milp_csv, sep = ",", row.names = FALSE, col.names = TRUE)
  cat("\nKM-MILP finalizado. Archivo guardado exitosamente en:", .pred_km_milp_csv, "\n")
} else {
  cat("\nKM-MILP finalizado sin resultados válidos. No se creó el CSV.\n")
}