# ==============================================================================
# PSO.R — PSO con Restricciones de Cardinalidad
# ==============================================================================

library(cluster)
library(pso)
library(proxy)
library(dplyr)

if (!requireNamespace("aricode", quietly = TRUE)) install.packages("aricode")
library(aricode)

# -----------------------------
# Registro interno de hiperparámetros (solo en memoria)
# Consultar: .pso_hyperparams$global / .pso_hyperparams$runs[["nombre"]]
# -----------------------------
.pso_hyperparams <- list(
  global = list(
    algorithm          = "PSO",
    distance_method    = "cosine",
    pso_maxit          = 20L,
    pso_swarm_size     = 40L,
    par_lower          = 0,
    par_upper          = 1,
    objective          = "intra_cluster_cosine_distance",
    assignment_method  = "order_based_partition"
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
.pred_pso_csv <- "predictions/pred_PSO.csv"

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
# PSO
# OPT: pre-calcular cumsum de cluster_sizes para asignación sin loop
# OPT: cost_function usa la cumsum pre-calculada
# -----------------------------
run_PSO <- function(X, target_cardinality) {
  
  set.seed(123)
  Xmat <- as.matrix(X)
  
  if (!all(is.finite(Xmat))) stop("X contiene NA/Inf")
  if (!all(vapply(as.data.frame(Xmat), is.numeric, logical(1)))) stop("X debe ser numérico")
  
  n <- nrow(Xmat)
  k <- length(target_cardinality)
  
  if (k < 2) stop("k < 2")
  if (sum(target_cardinality) != n) stop("sum(target_cardinality) != n (inconsistente)")
  
  D <- as.matrix(proxy::dist(Xmat, method = "cosine"))
  
  # OPT: pre-calcular los índices de corte para evitar loop en cost_function
  cs_sizes <- as.integer(target_cardinality)
  cum_sizes <- cumsum(cs_sizes)
  start_positions <- c(1L, cum_sizes[-k] + 1L)
  
  cost_function <- function(par, dist_matrix, cluster_sizes, k) {
    n <- length(par)
    cluster_assignment <- integer(n)
    
    ord <- order(par)
    # OPT: usar start/end pre-calculados
    for (i in 1:k) {
      cluster_assignment[ord[start_positions[i]:cum_sizes[i]]] <- i
    }
    
    total_distance <- 0
    for (i in 1:k) {
      ci <- which(cluster_assignment == i)
      if (length(ci) > 1) {
        cd <- dist_matrix[ci, ci, drop = FALSE]
        total_distance <- total_distance + sum(cd[lower.tri(cd, diag = FALSE)])
      }
    }
    total_distance
  }
  
  pso_result <- psoptim(
    par = runif(n),
    fn = cost_function,
    dist_matrix = D,
    cluster_sizes = cs_sizes,
    k = k,
    lower = 0,
    upper = 1,
    control = list(maxit = 20, s = 40)
  )
  
  label_pred <- integer(n)
  ord <- order(pso_result$par)
  # OPT: usar start/end pre-calculados
  for (i in seq_len(k)) {
    label_pred[ord[start_positions[i]:cum_sizes[i]]] <- i
  }
  
  list(y_predict = label_pred)
}

# -----------------------------
# Runner: devuelve fila
# -----------------------------
run_clustering_row <- function(dataset, target_cardinality, dataset_name) {
  
  set.seed(123)
  data <- prepare_data(dataset)
  if (is.null(data)) return(NULL)
  
  X <- data$X
  y <- data$y
  
  if (is.null(target_cardinality) || all(is.na(target_cardinality))) return(NULL)
  target_cardinality <- as.integer(target_cardinality)
  
  if (sum(target_cardinality) != nrow(X)) return(NULL)
  if (length(target_cardinality) < 2) return(NULL)
  
  start_total <- proc.time()
  
  results <- tryCatch(run_PSO(X, target_cardinality), error = function(e) NULL)
  
  if (is.null(results) || is.null(results$y_predict)) return(NULL)
  
  y_predict <- results$y_predict
  y_int <- as.integer(y)
  
  total_time <- (proc.time() - start_total)[3]
  
  metrics <- compute_metrics(y, y_predict, X)
  
  # Guardar hiperparámetros de esta ejecución (solo en memoria)
  .pso_hyperparams$runs[[dataset_name]] <<- list(
    k                = length(target_cardinality),
    n                = length(y_predict),
    n_features       = ncol(as.matrix(X)),
    target_sizes     = target_cardinality,
    pso_maxit        = 20L,
    pso_swarm_size   = 40L,
    par_lower        = 0,
    par_upper        = 1,
    distance_method  = "cosine",
    assignment_method = "order_based_partition"
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
# OPT: pre-asignar lista, seq_len, vapply
# -----------------------------
results_list <- vector("list", nrow(odatasets_unique))

for (i in seq_len(nrow(odatasets_unique))) {
  cat("\n\n--- Executing PSO for dataset at position:", i, "---\n")
  
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
results_clean <- results_list[!vapply(results_list, is.null, logical(1))]

if (length(results_clean) > 0) {
  final_df <- do.call(rbind, results_clean)
  write.table(final_df, .pred_pso_csv, sep = ",", row.names = FALSE, col.names = TRUE)
  cat("\nPSO finalizado. Archivo guardado exitosamente en:", .pred_pso_csv, "\n")
} else {
  cat("\nPSO finalizado sin resultados válidos. No se creó el CSV.\n")
}
