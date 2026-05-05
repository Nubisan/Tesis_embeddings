library(cluster)
library(proxy)

if (!requireNamespace("aricode", quietly = TRUE)) install.packages("aricode")
library(aricode)

# -----------------------------
# Registro interno de hiperparámetros 
# Consultar: .aco_hyperparams$global / .aco_hyperparams$runs[["nombre"]]
# -----------------------------
.aco_hyperparams <- list(
  global = list(
    algorithm          = "ACO",
    distance_method    = "cosine",
    n_ants             = 50L,
    max_iterations     = 20L,
    penalty_weight     = 100,
    perturbation_rate  = 0.1,
    kmeans_nstart      = 5L,
    kmeans_iter_max    = 30L,
    adjust_cardinality_max_iter = 1000L,
    objective          = "silhouette_cosine_minus_penalty"
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
.pred_aco_csv <- "predictions/pred_ACO.csv"

# -----------------------------
# Helpers: X numérico robusto
# -----------------------------
make_numeric_X <- function(X) {
  X <- as.data.frame(X)
  X_num <- lapply(X, function(col) {
    if (is.numeric(col) || is.integer(col)) return(as.numeric(col))
    if (is.logical(col))                    return(as.numeric(col))
    if (is.factor(col))                     return(as.numeric(col))
    if (is.character(col))                  return(as.numeric(as.factor(col)))
    return(NULL)
  })
  keep  <- !vapply(X_num, is.null, logical(1))
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
# CORRECCIÓN 1: calculate_centroids
# -----------------------------
calculate_centroids <- function(X, cluster_assignment, k) {
  X <- as.data.frame(X)
  centroids_df <- aggregate(X, by = list(cluster = cluster_assignment), FUN = mean)
  centroids_df <- centroids_df[order(centroids_df$cluster), , drop = FALSE]
  centroids    <- as.matrix(centroids_df[, -1, drop = FALSE])
  
  if (nrow(centroids) < k) {
    tmp        <- matrix(NA_real_, nrow = k, ncol = ncol(centroids))
    tmp[as.integer(centroids_df$cluster), ] <- centroids
    centroids  <- tmp
  }
  if (any(is.na(centroids))) {
    global_mu <- colMeans(as.matrix(X), na.rm = TRUE)
    for (r in seq_len(nrow(centroids))) {
      if (any(is.na(centroids[r, ]))) centroids[r, ] <- global_mu
    }
  }
  centroids
}

# -----------------------------
# CORRECCIÓN 2: adjust_cardinality
# OPT: tabulate en vez de sapply para contar tamaños
# -----------------------------
adjust_cardinality <- function(cluster_assignment, X, centroids, target_cardinality) {
  max_iterations <- 1000
  iteration      <- 0
  k              <- length(target_cardinality)
  # OPT: tabulate es O(n) en C, reemplaza sapply(1:k, function(c) sum(...))
  cluster_sizes  <- tabulate(cluster_assignment, nbins = k)
  
  while (any(cluster_sizes > target_cardinality) && iteration < max_iterations) {
    iteration <- iteration + 1
    
    for (j in which(cluster_sizes > target_cardinality)) {
      idx <- which(cluster_assignment == j)
      if (length(idx) == 0L) next
      
      element   <- idx[1]
      distances <- colSums((t(centroids) - as.numeric(X[element, ]))^2)
      
      available_clusters <- which(cluster_sizes < target_cardinality)
      valid_clusters     <- setdiff(available_clusters, j)
      
      if (length(valid_clusters) == 0L) {
        cluster_assignment[element] <- -1L
        cluster_sizes[j]            <- cluster_sizes[j] - 1L
      } else {
        chosen                       <- valid_clusters[which.min(distances[valid_clusters])]
        cluster_assignment[element]  <- chosen
        cluster_sizes[j]             <- cluster_sizes[j] - 1L
        cluster_sizes[chosen]        <- cluster_sizes[chosen] + 1L
      }
    }
  }
  
  unassigned_idx <- which(cluster_assignment == -1L)
  if (length(unassigned_idx) > 0L) {
    for (element in unassigned_idx) {
      available_clusters <- which(cluster_sizes < target_cardinality)
      if (length(available_clusters) > 0L) {
        chosen                      <- available_clusters[which.min(cluster_sizes[available_clusters])]
        cluster_assignment[element] <- chosen
        cluster_sizes[chosen]       <- cluster_sizes[chosen] + 1L
      } else {
        cluster_assignment[element] <- 1L
        cluster_sizes[1]            <- cluster_sizes[1] + 1L
      }
    }
  }
  
  if (iteration >= max_iterations) {
    warning("adjust_cardinality alcanzó el límite de iteraciones.")
  }
  
  cluster_assignment
}

# -----------------------------
# ACO principal
# OPT: dist coseno NxN se calcula una sola vez y se pasa a evaluate_solution
# OPT: penalty vectorizada con tabulate
# OPT: perturbación vectorizada sin loop punto por punto
# -----------------------------
run_ACO <- function(X, target_cardinality) {
  set.seed(123)
  Xmat <- as.matrix(X)
  n    <- nrow(Xmat)
  k    <- length(target_cardinality)
  
  if (is.na(k) || k < 2)  stop("k inválido (<2).")
  if (k >= n)              stop(paste0("k=", k, " inválido para n=", n, "."))
  
  n_ants         <- 50
  max_iterations <- 20
  penalty_weight <- 100
  
  # OPT: calcular dist coseno NxN una sola vez (antes: se recalculaba en cada evaluate_solution)
  d_cosine <- proxy::dist(Xmat, method = "cosine")
  
  # Función objetivo: silhouette coseno - penalización por violación
  # OPT: recibe d_cosine pre-calculada
  evaluate_solution <- function(cluster_assignment) {
    if (length(unique(cluster_assignment)) < 2) return(-Inf)
    ss <- tryCatch(cluster::silhouette(cluster_assignment, d_cosine), error = function(e) NULL)
    if (is.null(ss)) return(-Inf)
    # OPT: tabulate + sum vectorizado reemplaza for-loop de penalty
    current_counts <- tabulate(cluster_assignment, nbins = k)
    penalty <- penalty_weight * sum(abs(current_counts - target_cardinality))
    mean(ss[, "sil_width"]) - penalty
  }
  
  # Solución inicial: kmeans + ajuste heurístico simple
  generate_initial_solution <- function() {
    km                 <- kmeans(Xmat, centers = k, nstart = 5, iter.max = 30)
    cluster_assignment <- km$cluster
    for (j in seq_len(k)) {
      while (sum(cluster_assignment == j) > target_cardinality[j]) {
        idx                            <- which(cluster_assignment == j)
        cluster_assignment[sample(idx, 1)] <- sample(setdiff(1:k, j), 1)
      }
    }
    cluster_assignment
  }
  
  # Perturbación: cambia ~10% de los puntos aleatoriamente
  # Loop original preservado: runif(1) por punto mantiene la secuencia exacta del RNG
  perturb_solution <- function(cluster_assignment) {
    new_ca <- cluster_assignment
    for (j in 1:n) {
      if (runif(1) < 0.1) new_ca[j] <- sample(1:k, 1)
    }
    new_ca
  }
  
  best_score             <- -Inf
  best_cluster_assignment <- NULL
  
  for (iteration in seq_len(max_iterations)) {
    ants <- vector("list", n_ants)
    for (i in seq_len(n_ants)) {
      ca       <- generate_initial_solution()
      ca       <- perturb_solution(ca)
      ants[[i]] <- ca
    }
    for (i in seq_len(n_ants)) {
      score <- evaluate_solution(ants[[i]])
      if (is.finite(score) && score > best_score) {
        best_score              <- score
        best_cluster_assignment <- ants[[i]]
      }
    }
  }
  
  if (!is.null(best_cluster_assignment)) {
    centroids               <- calculate_centroids(Xmat, best_cluster_assignment, k)
    best_cluster_assignment <- adjust_cardinality(
      best_cluster_assignment,
      Xmat,
      centroids,
      target_cardinality
    )
  }
  
  list(best_solution = best_cluster_assignment, best_score = best_score)
}

# -----------------------------
# Runner: genera una fila de resultado por dataset
# -----------------------------
run_clustering_row <- function(dataset, target_cardinality, dataset_name) {
  data <- prepare_data(dataset)
  if (is.null(data)) return(NULL)
  X <- data$X
  y <- data$y
  if (is.null(target_cardinality) || length(target_cardinality) < 2 || all(is.na(target_cardinality))) return(NULL)
  
  set.seed(123)
  start_total <- proc.time()
  
  results    <- tryCatch(run_ACO(X, as.integer(target_cardinality)), error = function(e) NULL)
  
  if (is.null(results) || is.null(results$best_solution)) return(NULL)
  
  y_predict <- results$best_solution
  
  total_time <- (proc.time() - start_total)[3]
  
  metrics <- compute_metrics(y, y_predict, X)
  
  # Guardar hiperparámetros de esta ejecución (solo en memoria)
  .aco_hyperparams$runs[[dataset_name]] <<- list(
    k                 = length(target_cardinality),
    n                 = length(y_predict),
    n_features        = ncol(as.matrix(X)),
    target_sizes      = as.integer(target_cardinality),
    n_ants            = 50L,
    max_iterations    = 20L,
    penalty_weight    = 100,
    perturbation_rate = 0.1,
    distance_method   = "cosine"
  )
  
  data.frame(
    name = dataset_name,
    n = length(y_predict),
    k = length(unique(y_predict)),
    y_predict = paste(as.integer(y_predict), collapse = " "),
    y_true = paste(as.integer(y), collapse = " "), 
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
  cat("\n\n--- Executing ACO for dataset at position:", i, "---\n")
  
  tryCatch({
    dataset            <- odatasets_unique$dataset[[i]]
    dataset_name       <- odatasets_unique$name[[i]]
    target_cardinality <- odatasets_unique$class_distribution_vector[[i]]
    
    if (is.null(target_cardinality) || all(is.na(target_cardinality)) || length(target_cardinality) < 2) {
      cat("Target cardinality inválida. Skipping.\n")
      next
    }
    
    row_result <- run_clustering_row(dataset, target_cardinality, dataset_name)
    
    if (!is.null(row_result)) {
      results_list[[i]] <- row_result
      cat("Time:", row_result$Execution_Time, "s\n")
    } else {
      cat("Dataset skipped.\n")
    }
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}

results_clean <- results_list[!vapply(results_list, is.null, logical(1))]

if (length(results_clean) > 0) {
  write.table(do.call(rbind, results_clean), .pred_aco_csv,
              sep = ",", row.names = FALSE, col.names = TRUE)
  cat("\nACO guardado en:", .pred_aco_csv, "\n")
} else {
  cat("\nACO finalizado sin resultados válidos. No se generó CSV.\n")
}