#' Creates a set of n_runs embeddings using either BH t-SNE or UMAP
#' @description Using either BH t-SNE or UMAP, creates a list of n_runs embeddings and returns the set of embeddings as a list of dataframes
#' @param dr_input Numeric matrix or data frame for DR (rows = observations)
#' @param n_runs Integer >= 1; number of embeddings to create
#' @param method "umap" or "tsne"
#' @param n_neighbors Positive numeric (UMAP)
#' @param min_dist Non-negative numeric (UMAP)
#' @param perplexity Positive numeric (t-SNE)
#' @param theta numeric value between (0,1) (t-SNE)
#' @return A list of length <= n_runs containing DR embeddings (each a matrix).
createEmb <- function(dr_input, n_runs = 100, method = c("umap","tsne"),
                      n_neighbors = 15, min_dist = 0.1,
                      perplexity = 30, theta = 0.5) {
  # Input validation
  if (!(is.matrix(dr_input) || is.data.frame(dr_input))) {
    stop("`dr_input` must be a matrix or data frame of numeric values.", call. = FALSE)
  }
  dr_input <- as.matrix(dr_input)
  if (!is.numeric(dr_input)) stop("`dr_input` must be numeric.", call. = FALSE)

  if (!(is.numeric(n_runs) && length(n_runs)==1 && n_runs>=1 && n_runs==as.integer(n_runs))) {
    stop("`n_runs` must be a single integer >= 1.", call. = FALSE)
  }
  method <- match.arg(method)

  # Dependency checks
  for (pkg in c("uwot","Rtsne","future.apply")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
    }
  }

  # Plan management
  old_plan <- future::plan()
  future::plan(future::multisession, workers = parallel::detectCores() - 1)
  on.exit(future::plan(old_plan), add = TRUE)

  # Seeds and chunking
  seeds  <- seq(100, length.out = n_runs)
  chunks <- if (n_runs >= 100) parallel::detectCores()*4 else parallel::detectCores()*2

  # Parallel embedding with error handling
  results <- future.apply::future_lapply(
    seq_len(n_runs), function(i) {
      set.seed(seeds[i])
      tryCatch({
        if (method == "umap") {
          uwot::umap(dr_input, n_neighbors = n_neighbors,
                     min_dist = min_dist, n_threads = 1, ret_model = FALSE)
        } else {
          Rtsne::Rtsne(dr_input, perplexity = perplexity,
                       theta = theta, pca = FALSE,
                       check_duplicates = FALSE,
                       max_iter = 500, num_threads = 1)$Y
        }
      }, error = function(e) {
        warning(sprintf("%s run %d failed: %s", toupper(method), i, e$message),
                call. = FALSE)
        NULL
      })
    }, future.seed = TRUE, future.chunk.size = chunks
  )
  gc()

  # Post-process results
  failed <- sum(vapply(results, is.null, logical(1)))
  if (failed > 0) {
    warning(sprintf("%d/%d embeddings failed and were omitted.", failed, n_runs),
            call. = FALSE)
  }
  embeddings <- results[!vapply(results, is.null, logical(1))]

  if (length(embeddings) == 0) {
    stop("All embedding runs failed; no results to return.", call. = FALSE)
  }

  return(embeddings)
}



#' Compares a list of DR embeddings and returns stability statistics
#' @description Compares a list of input embeddings and aligns them pairwise using procrustes, then computes a Kendall's Tau correlation between each pairwise alignment. Returns the mean of means correlation, density of mean correlations, range of correlations and 95% CI of means.
#' @param emb_list A list of embeddings created by createEmb
#' @importFrom stats quantile
#' @importFrom stats density
#' @importFrom utils combn
#' @importFrom pcaPP cor.fk
#' @importFrom vegan procrustes
#' @return Mean of mean correlation, density of mean correlation per embedding, range of mean correlation per embedding, 95% CI of mean correlation per embedding
compareEmb <- function(emb_list) {
  # Input checks
  utils::globalVariables("n")

  if (!is.list(emb_list) || length(emb_list) < 2) {
    stop("`emb_list` must be a list of at least two embeddings.", call. = FALSE)
  }
  dims <- unique(lapply(emb_list, dim))
  if (length(dims) != 1) {
    stop("All embeddings must have identical dimensions.", call. = FALSE)
  }

  # Dependency checks
  for (pkg in c("vegan", "pcaPP", "future.apply", "stats")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required.", pkg), call. = FALSE)
    }
  }
  # Setup
  matList <- lapply(emb_list, as.matrix)
  pairs   <- combn(length(matList), 2, simplify = FALSE)
  future::plan(future::multisession)

  # Align and compare embeddings
  cor_values <- future.apply::future_sapply(pairs, function(idx_pair) {
    i <- idx_pair[1]; j <- idx_pair[2]
    proc <- tryCatch(vegan::procrustes(matList[[i]], matList[[j]], symmetric = TRUE),
                     error = function(e) {
                       warning(sprintf("Procrustes failed for (%d, %d): %s", i, j, e$message), call. = FALSE)
                       return(NULL)
                     })
    if (is.null(proc)) return(NA_real_)
    vec1 <- as.numeric(proc$X); vec2 <- as.numeric(proc$Yrot)
    tryCatch(pcaPP::cor.fk(vec1, vec2),
             error = function(e) {
               warning(sprintf("Correlation failed for (%d, %d): %s", i, j, e$message), call. = FALSE)
               NA_real_
             })
  }, future.seed = TRUE)
  future::plan(NULL)

  # Use the computed cor_values to fill in the correlation matrix.
  for (k in seq_along(pairs)) {
    i <- pairs[[k]][1]
    j <- pairs[[k]][2]
    corr_mat[i, j] <- cor_values[k]
    corr_mat[j, i] <- cor_values[k]
  }

  # Compute the mean correlation per embedding (row means excluding the diagonal).
  mean_per_embedding <- sapply(1:n, function(i) {
    mean(corr_mat[i, -i], na.rm = TRUE)
  })

  # Compute overall mean correlation.
  cor_mean <- mean(mean_per_embedding)

  # Find CI based on quantile
  ci <- quantile(x = mean_per_embedding, probs = c(0.025, 0.975), na.rm = TRUE)

  range_vals <- range(mean_per_embedding)

  # Generate density of means
  d <- density(mean_per_embedding)

  cat(sprintf("Mean of Mean Correlation per Embedding: %.4f\n", cor_mean))
  cat(sprintf("95 Percent CI of Mean Correlation per Embedding: [%.4f, %.4f]\n", ci[1], ci[2]))
  cat(sprintf("Correlation Range: Min = %.4f Max = %.4f\n", range_vals[1], range_vals[2]))

  plot(d, type = "l", main = "Density of Mean Correlations per Embedding", sub = sprintf("%f embeddings", n),
       xlab = "Mean Correlation Per Embedding", ylab = "Density")


  # Return the results.
  return(list(all_pairwise_correlations = corr_mat,
              mean_per_embedding = mean_per_embedding,
              mean = cor_mean,
              ci = ci,
              range = range_vals))
}
