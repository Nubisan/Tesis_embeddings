library(lpSolve)
library(proxy)
library(dplyr)

dir.create("predictions", showWarnings = FALSE, recursive = FALSE)
.pred_sck1_emb_csv <- "predictions/pred_SCK1_emb.csv"

# ==============================================================================
# BLOQUE 1 — PREPARACIÓN DE DATOS
# ==============================================================================

prepare_embeddings <- function(emb_matrix, y_labels) {
  X <- as.matrix(emb_matrix)

  # Verificaciones básicas
  if (!is.numeric(X))       stop("emb_matrix debe ser numérica.")
  if (any(!is.finite(X)))   stop("emb_matrix contiene NA o Inf.")
  if (length(y_labels) != nrow(X)) stop("y_labels debe tener la misma longitud que filas en emb_matrix.")

  # Normalizar cada fila a norma unitaria
  # rowSums(X^2) calcula la suma de cuadrados de cada fila → norma al cuadrado
  # sqrt(...) → norma euclidiana de cada vector
  # pmax(..., 1e-10) evita división por cero si algún vector es todo ceros
  norms <- sqrt(rowSums(X^2))
  X_norm <- X / pmax(norms, 1e-10)

  list(X = X_norm, y = as.factor(y_labels))
}

# ==============================================================================
# BLOQUE 2 — ILP 
# ==============================================================================

solve_ilp_assignment <- function(cost_mat, sizes) {
  cost_mat <- as.matrix(cost_mat)
  k <- nrow(cost_mat)   # número de clusters
  n <- ncol(cost_mat)   # número de puntos
  sizes <- as.integer(sizes)

  if (length(sizes) != k)  stop("sizes length != k")
  if (sum(sizes) != n)     stop("sum(sizes) != n — las cardinalidades no suman n")
  if (k < 2)               stop("k debe ser >= 2")

  # Aplanar la matriz de costos en un vector (fila por fila)
  # Si cost_mat es [[c11,c12],[c21,c22]], cost_vec = [c11,c12,c21,c22]
  cost_vec <- as.vector(t(cost_mat))

  # Restricción A1: cada punto asignado a exactamente 1 cluster
  A1 <- matrix(0, n, k * n)
  for (j in 1:n) {
    A1[j, seq(j, k * n, by = n)] <- 1
  }

  # Restricción A2: cada cluster con exactamente sizes[i] puntos
  A2 <- matrix(0, k, k * n)
  for (i in 1:k) {
    A2[i, ((i - 1) * n + 1):(i * n)] <- 1
  }

  # Resolver con lpSolve — todas las variables binarias (0 o 1)
  res <- lp(
    direction    = "min",
    objective.in = cost_vec,
    const.mat    = rbind(A1, A2),
    const.dir    = c(rep("=", n), rep("=", k)),
    const.rhs    = c(rep(1, n), sizes),
    all.bin      = TRUE
  )

  if (res$status != 0) stop("ILP no encontró solución factible.")

  # Reformatear solución: matrix k×n donde assignment[i,j]=1 si punto j → cluster i
  matrix(res$solution, nrow = k, byrow = TRUE)
}

# ==============================================================================
# BLOQUE 3 — SCK1 PRINCIPAL
# ==============================================================================

sck1_iterativo <- function(data, k, sizes, max_iter = 5) {
  data <- as.matrix(data)
  n    <- nrow(data)
  d    <- ncol(data)

  if (n != sum(sizes)) stop("nrow(data) != sum(sizes)")
  if (k < 2)           stop("k debe ser >= 2")

  # CAMBIO: inicialización con kmeans en lugar de índices aleatorios
  km_init      <- kmeans(data, centers = k, nstart = 3, iter.max = 20)
  centroids    <- km_init$centers

  assignments  <- NULL

  for (iter in 1:max_iter) {
    # Calcular distancias de cada punto a cada centroide
    cost_mat    <- t(as.matrix(proxy::dist(data, centroids, method = "cosine")))

    # Resolver ILP para asignación exacta
    assign_mat  <- solve_ilp_assignment(cost_mat, sizes)

    # Convertir matriz de asignación binaria a vector de etiquetas
    assignments <- apply(assign_mat, 2, which.max)

    # Recalcular centroides como media de los puntos asignados
    new_centroids <- matrix(0, nrow = k, ncol = d)
    for (i in 1:k) {
      pts <- data[assignments == i, , drop = FALSE]
      if (nrow(pts) > 0) new_centroids[i, ] <- colMeans(pts)
      else               new_centroids[i, ] <- centroids[i, ]
    }

    # Verificar convergencia: si los centroides casi no se movieron, parar
    if (!is.null(assignments) &&
        max(abs(centroids - new_centroids)) < 1e-6) break

    centroids <- new_centroids
  }

  list(assignments = assignments, centroids = centroids)
}

# ==============================================================================
# BLOQUE 4 — RUNNER 
# ==============================================================================

run_clustering_row <- function(emb_matrix, y_labels, target_cardinality,
                               dataset_name, seed = 123) {

  # Preparar embeddings normalizados
  data <- prepare_embeddings(emb_matrix, y_labels)
  if (is.null(data)) return(NULL)

  X <- data$X   # matriz n×d normalizada
  y <- data$y   # etiquetas reales (factor)

  # Validaciones de cardinalidad
  target_cardinality <- as.integer(target_cardinality)
  k <- length(target_cardinality)
  if (k < 2)                            return(NULL)
  if (sum(target_cardinality) != nrow(X)) return(NULL)

  set.seed(seed)

  start_algo <- Sys.time()
  res <- tryCatch(
    sck1_iterativo(X, k, target_cardinality),
    error = function(e) {
      cat("Error en sck1_iterativo:", e$message, "\n")
      NULL
    }
  )
  end_algo <- Sys.time()

  if (is.null(res) || is.null(res$assignments)) return(NULL)

  data.frame(
    name               = dataset_name,
    n                  = length(res$assignments),
    k                  = length(unique(res$assignments)),
    y_predict          = paste(as.integer(res$assignments),  collapse = " "),
    y_true             = paste(as.integer(y),                collapse = " "),
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
  cat("\n\n--- Executing SCK1-EMB for dataset at position:", i, "---\n")

  tryCatch({
    dataset            <- odatasets_unique$dataset[[i]]
    dataset_name       <- odatasets_unique$name[[i]]
    target_cardinality <- odatasets_unique$class_distribution_vector[[i]]

    if (is.null(dataset) || is.null(target_cardinality) ||
        all(is.na(target_cardinality))) {
      cat("Dataset o cardinalidad inválida. Skipping.\n")
      next
    }

    # Extraer X e y del dataset 
    # Si en el futuro usas embeddings reales, reemplaza estas dos líneas
    # por la carga de tu matriz de embeddings
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
  write.table(do.call(rbind, results_clean), .pred_sck1_emb_csv,
              sep = ",", row.names = FALSE, col.names = TRUE)
  cat("\nSCK1-EMB guardado en:", .pred_sck1_emb_csv, "\n")
} else {
  cat("\nSCK1-EMB finalizado sin resultados válidos.No se creó el CSV.\n")
}
