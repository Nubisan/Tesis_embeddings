library(cluster)
library(proxy)

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
# Necesaria para que adjust_cardinality pueda calcular
# distancias al momento de reasignar puntos
# -----------------------------
calculate_centroids <- function(X, cluster_assignment, k) {
  X <- as.data.frame(X)
  centroids_df <- aggregate(X, by = list(cluster = cluster_assignment), FUN = mean)
  centroids_df <- centroids_df[order(centroids_df$cluster), , drop = FALSE]
  centroids    <- as.matrix(centroids_df[, -1, drop = FALSE])
  
  # Si falta algún cluster (vacío), rellenar con la media global
  if (nrow(centroids) < k) {
    tmp        <- matrix(NA_real_, nrow = k, ncol = ncol(centroids))
    tmp[as.integer(centroids_df$cluster), ] <- centroids
    centroids  <- tmp
  }
  if (any(is.na(centroids))) {
    global_mu <- colMeans(as.matrix(X), na.rm = TRUE)
    for (r in 1:nrow(centroids)) {
      if (any(is.na(centroids[r, ]))) centroids[r, ] <- global_mu
    }
  }
  centroids
}

# -----------------------------
# CORRECCIÓN 2: adjust_cardinality
# Garantiza que la solución final cumpla exactamente
# las cardinalidades objetivo. Esta función existía en
# BAT, CSCLP, KMEDOIDS, etc., pero faltaba en ACO.
# -----------------------------
adjust_cardinality <- function(cluster_assignment, X, centroids, target_cardinality) {
  max_iterations <- 1000
  iteration      <- 0
  k              <- length(target_cardinality)
  cluster_sizes  <- sapply(1:k, function(c) sum(cluster_assignment == c, na.rm = TRUE))
  
  while (any(cluster_sizes > target_cardinality) && iteration < max_iterations) {
    iteration <- iteration + 1
    
    for (j in which(cluster_sizes > target_cardinality)) {
      idx <- which(cluster_assignment == j)
      if (length(idx) == 0) next
      
      element   <- idx[1]
      # Distancia del punto a todos los centroides (vectorizado)
      distances <- colSums((t(centroids) - as.numeric(X[element, ]))^2)
      
      available_clusters <- which(cluster_sizes < target_cardinality)
      valid_clusters     <- setdiff(available_clusters, j)
      
      if (length(valid_clusters) == 0) {
        # Sin espacio disponible: marcar como no asignado temporalmente
        cluster_assignment[element] <- -1
        cluster_sizes[j]            <- cluster_sizes[j] - 1
      } else {
        # Reasignar al cluster válido más cercano
        chosen                       <- valid_clusters[which.min(distances[valid_clusters])]
        cluster_assignment[element]  <- chosen
        cluster_sizes[j]             <- cluster_sizes[j] - 1
        cluster_sizes[chosen]        <- cluster_sizes[chosen] + 1
      }
    }
  }
  
  # Reasignar puntos marcados como -1
  unassigned_idx <- which(cluster_assignment == -1)
  if (length(unassigned_idx) > 0) {
    for (element in unassigned_idx) {
      available_clusters <- which(cluster_sizes < target_cardinality)
      if (length(available_clusters) > 0) {
        chosen                      <- available_clusters[which.min(cluster_sizes[available_clusters])]
        cluster_assignment[element] <- chosen
        cluster_sizes[chosen]       <- cluster_sizes[chosen] + 1
      } else {
        # Caso extremo: asignar al cluster 1
        cluster_assignment[element] <- 1
        cluster_sizes[1]            <- cluster_sizes[1] + 1
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
# -----------------------------
run_ACO <- function(X, target_cardinality) {
  Xmat <- as.matrix(X)
  n    <- nrow(Xmat)
  k    <- length(target_cardinality)
  
  if (is.na(k) || k < 2)  stop("k inválido (<2).")
  if (k >= n)              stop(paste0("k=", k, " inválido para n=", n, "."))
  
  n_ants         <- 50
  max_iterations <- 20
  penalty_weight <- 100
  
  # Función objetivo: silhouette coseno - penalización por violación
  evaluate_solution <- function(cluster_assignment) {
    if (length(unique(cluster_assignment)) < 2) return(-Inf)
    d  <- proxy::dist(Xmat, method = "cosine")
    ss <- tryCatch(cluster::silhouette(cluster_assignment, d), error = function(e) NULL)
    if (is.null(ss)) return(-Inf)
    penalty <- 0
    for (j in 1:k) {
      diff <- abs(sum(cluster_assignment == j) - target_cardinality[j])
      if (diff > 0) penalty <- penalty + diff * penalty_weight
    }
    mean(ss[, "sil_width"]) - penalty
  }
  
  # Solución inicial: kmeans + ajuste heurístico simple
  generate_initial_solution <- function() {
    km                 <- kmeans(Xmat, centers = k, nstart = 5, iter.max = 30)
    cluster_assignment <- km$cluster
    for (j in 1:k) {
      while (sum(cluster_assignment == j) > target_cardinality[j]) {
        idx                            <- which(cluster_assignment == j)
        cluster_assignment[sample(idx, 1)] <- sample(setdiff(1:k, j), 1)
      }
    }
    cluster_assignment
  }
  
  # Perturbación: cambia ~10% de los puntos aleatoriamente
  perturb_solution <- function(cluster_assignment) {
    new_ca <- cluster_assignment
    for (j in 1:n) {
      if (runif(1) < 0.1) new_ca[j] <- sample(1:k, 1)
    }
    new_ca
  }
  
  best_score             <- -Inf
  best_cluster_assignment <- NULL
  
  for (iteration in 1:max_iterations) {
    ants <- vector("list", n_ants)
    for (i in 1:n_ants) {
      ca       <- generate_initial_solution()
      ca       <- perturb_solution(ca)
      ants[[i]] <- ca
    }
    for (i in 1:n_ants) {
      score <- evaluate_solution(ants[[i]])
      if (is.finite(score) && score > best_score) {
        best_score              <- score
        best_cluster_assignment <- ants[[i]]
      }
    }
  }
  
  # --------------------------------------------------
  # CORRECCIÓN PRINCIPAL:
  # Aplicar adjust_cardinality sobre la mejor solución
  # encontrada. Sin este paso, ACO puede violar las
  # restricciones de tamaño en hasta el 47% de los
  # datasets (34/72 observado en los resultados).
  # --------------------------------------------------
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
  
  start_algo <- Sys.time()
  results    <- tryCatch(run_ACO(X, as.integer(target_cardinality)), error = function(e) NULL)
  end_algo   <- Sys.time()
  
  if (is.null(results) || is.null(results$best_solution)) return(NULL)
  
  y_predict <- results$best_solution
  
  data.frame(
    name               = dataset_name,
    n                  = length(y_predict),
    k                  = length(unique(y_predict)),
    y_predict          = paste(as.integer(y_predict),         collapse = " "),
    y_true             = paste(as.integer(y),                  collapse = " "),
    target_cardinality = paste(as.integer(target_cardinality), collapse = " "),
    Execution_Time     = as.numeric(difftime(end_algo, start_algo, units = "secs")),
    stringsAsFactors   = FALSE
  )
}

# -----------------------------
# Loop principal
# -----------------------------
results_list <- list()

for (i in 1:nrow(odatasets_unique)) {
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

results_clean <- results_list[!sapply(results_list, is.null)]

if (length(results_clean) > 0) {
  write.table(do.call(rbind, results_clean), .pred_aco_csv,
              sep = ",", row.names = FALSE, col.names = TRUE)
  cat("\nACO guardado en:", .pred_aco_csv, "\n")
} else {
  cat("\nACO finalizado sin resultados válidos. No se generó CSV.\n")
}