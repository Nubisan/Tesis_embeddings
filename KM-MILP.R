library(lpSolve)
library(proxy)
library(dplyr)

dir.create("predictions", showWarnings = FALSE, recursive = FALSE)
.pred_milp_emb_csv <- "predictions/pred_KM-MILP_emb.csv"

# ==============================================================================
# BLOQUE 1 — PREPARACIÓN DE DATOS (NUEVA)
# ==============================================================================
prepare_embeddings <- function(emb_matrix, y_labels) {
  X <- as.matrix(emb_matrix)

  if (!is.numeric(X))       stop("emb_matrix debe ser numérica.")
  if (any(!is.finite(X)))   stop("emb_matrix contiene NA o Inf.")
  if (length(y_labels) != nrow(X)) stop("y_labels debe tener misma longitud que filas.")

  norms  <- sqrt(rowSums(X^2))
  X_norm <- X / pmax(norms, 1e-10)

  list(X = X_norm, y = as.factor(y_labels))
}

# ------------------------------------------------------------------------------
# FUNCIÓN AUXILIAR NUEVA: renormalize_centroids()
# ------------------------------------------------------------------------------
renormalize_centroids <- function(centroids) {
  norms <- sqrt(rowSums(centroids^2))
  centroids / pmax(norms, 1e-10)
}

# ==============================================================================
# BLOQUE 2 — MILP 
# ==============================================================================

solve_milp_assignment <- function(data, centroids, size_constraints) {
  n <- nrow(data)
  k <- nrow(centroids)

  # Calcular matriz de costos n×k: distancia coseno de cada punto a cada centroide
  # Esta es la operación central — n×k en lugar de n×n
  # Con n=10,000 y k=5: 50,000 valores vs 100,000,000 — 2000 veces menos memoria
  cost_matrix <- as.matrix(proxy::dist(data, centroids, method = "cosine"))
  cost_vec    <- as.vector(t(cost_matrix))   # aplanar fila por fila

  # Variables: x_ij ∈ {0,1}, tamaño n×k
  f.obj <- cost_vec

  # Restricción 1: cada punto asignado a exactamente 1 cluster
  constr1 <- matrix(0, n, n * k)
  for (i in 1:n) {
    constr1[i, ((i - 1) * k + 1):(i * k)] <- 1
  }

  # Restricción 2: cada cluster con tamaño <= size_constraints[j]
  constr2 <- matrix(0, k, n * k)
  for (j in 1:k) {
    constr2[j, seq(j, n * k, by = k)] <- 1
  }

  f.con <- rbind(constr1, constr2)
  f.dir <- c(rep("=", n), rep("<=", k))   # <= para MILP (diferente a SCK1)
  f.rhs <- c(rep(1, n), size_constraints)

  result <- lp("min", f.obj, f.con, f.dir, f.rhs, all.bin = TRUE)

  if (result$status != 0) stop("MILP no encontró solución.")

  # Extraer asignaciones: p[i] = cluster del punto i
  x_opt <- matrix(result$solution, nrow = n, byrow = TRUE)
  p     <- apply(x_opt, 1, which.max)

  list(p = p)
}

# ==============================================================================
# BLOQUE 3 — KM-MILP PRINCIPAL (cambios en inicialización y centroides)
# ==============================================================================

clustering_with_size_constraints <- function(data, size_constraints,
                                             max_iter = 100, tol = 1e-6) {
  data <- as.matrix(data)
  n    <- nrow(data)
  d    <- ncol(data)
  k    <- length(size_constraints)

  if (sum(size_constraints) != n) stop("sum(size_constraints) != n")
  if (k < 2) stop("k debe ser >= 2")

  # CAMBIO 1: inicialización con kmeans en lugar de random_indices
  km_init   <- kmeans(data, centers = k, nstart = 3, iter.max = 30)
  centroids <- km_init$centers

  # CAMBIO 2: renormalizar centroides iniciales
  centroids <- renormalize_centroids(centroids)

  converged  <- FALSE
  iteration  <- 0
  p          <- NULL

  while (!converged && iteration < max_iter) {
    iteration <- iteration + 1

    # Asignación MILP — el núcleo del algoritmo
    milp_result <- solve_milp_assignment(data, centroids, size_constraints)
    p           <- milp_result$p

    # Actualizar centroides: media de puntos asignados a cada cluster
    new_centroids <- matrix(0, nrow = k, ncol = d)
    for (j in 1:k) {
      cluster_pts <- data[p == j, , drop = FALSE]
      if (nrow(cluster_pts) > 0) {
        new_centroids[j, ] <- colMeans(cluster_pts)
      } else {
        new_centroids[j, ] <- centroids[j, ]   # mantener si cluster vacío
      }
    }

    # CAMBIO 3: renormalizar después de recalcular
    new_centroids <- renormalize_centroids(new_centroids)

    # Verificar convergencia
    if (max(abs(centroids - new_centroids)) < tol) {
      converged <- TRUE
    }

    centroids <- new_centroids
  }

  cat(sprintf("  KM-MILP convergió en %d iteraciones\n", iteration))
  list(p = p, centroids = centroids)
}

# ==============================================================================
# BLOQUE 4 — RUNNER
# ==============================================================================

run_clustering_row <- function(emb_matrix, y_labels, target_cardinality,
                               dataset_name, seed = 123) {

  data <- prepare_embeddings(emb_matrix, y_labels)
  if (is.null(data)) return(NULL)

  X <- data$X
  y <- data$y

  target_cardinality <- as.integer(target_cardinality)
  k <- length(target_cardinality)

  if (k < 2)                              return(NULL)
  if (sum(target_cardinality) != nrow(X)) return(NULL)

  set.seed(seed)

  start_algo <- Sys.time()
  res <- tryCatch(
    clustering_with_size_constraints(X, target_cardinality),
    error = function(e) {
      cat("Error en clustering_with_size_constraints:", e$message, "\n")
      NULL
    }
  )
  end_algo <- Sys.time()

  if (is.null(res) || is.null(res$p)) return(NULL)

  data.frame(
    name               = dataset_name,
    n                  = length(res$p),
    k                  = length(unique(res$p)),
    y_predict          = paste(as.integer(res$p),              collapse = " "),
    y_true             = paste(as.integer(y),                  collapse = " "),
    target_cardinality = paste(as.integer(target_cardinality), collapse = " "),
    Execution_Time     = as.numeric(difftime(end_algo, start_algo, units = "secs")),
    stringsAsFactors   = FALSE
  )
}

# ==============================================================================
# BLOQUE 5 — LOOP PRINCIPAL
# ==============================================================================

results_list <- list()

for (i in 1:nrow(odatasets_unique)) {
  cat("\n\n--- Executing KM-MILP-EMB for dataset at position:", i, "---\n")

  tryCatch({
    dataset            <- odatasets_unique$dataset[[i]]
    dataset_name       <- odatasets_unique$name[[i]]
    target_cardinality <- odatasets_unique$class_distribution_vector[[i]]

    if (is.null(dataset) || is.null(target_cardinality) ||
        all(is.na(target_cardinality))) {
      cat("Dataset o cardinalidad inválida. Skipping.\n")
      next
    }

    # Extraer X e y — reemplazar estas líneas con la carga de embeddings reales
    ds <- as.data.frame(dataset)
    if ("class" %in% colnames(ds)) {
      y_labels  <- ds$class
      emb_input <- ds[, setdiff(colnames(ds), "class"), drop = FALSE]
    } else {
      y_labels  <- ds[[ncol(ds)]]
      emb_input <- ds[, -ncol(ds), drop = FALSE]
    }
    emb_input <- as.matrix(emb_input[, sapply(emb_input, is.numeric), drop = FALSE])

    row_result <- run_clustering_row(
      emb_matrix         = emb_input,
      y_labels           = y_labels,
      target_cardinality = target_cardinality,
      dataset_name       = dataset_name,
      seed               = 123
    )

    if (!is.null(row_result)) {
      results_list[[i]] <- row_result
      cat("Time:", row_result$Execution_Time, "s\n")
    } else {
      cat("Dataset skipped.\n")
    }

  }, error = function(e) cat("Error:", e$message, "\n"))
}

# Escritura final
results_clean <- results_list[!sapply(results_list, is.null)]
if (length(results_clean) > 0) {
  write.table(do.call(rbind, results_clean), .pred_milp_emb_csv,
              sep = ",", row.names = FALSE, col.names = TRUE)
  cat("\nKM-MILP-EMB guardado en:", .pred_milp_emb_csv, "\n")
} else {
  cat("\nKM-MILP-EMB finalizado sin resultados válidos. No se creó el CSV.\n")
}
