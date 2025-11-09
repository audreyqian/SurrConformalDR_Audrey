#' ============================================================================
#' Complete Implementation: Surrogate-assisted Clustered Conformal Inference
#' with Enhanced Cluster Selection Methods (BIC, AIC, Gap Statistic)
#' and GMM Clustering Option
#' ============================================================================

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS FOR OPTIMAL CLUSTER NUMBER DETERMINATION
# -----------------------------------------------------------------------------

#' Compute BIC for k-means clustering
#' @param kmeans_obj A kmeans object
#' @param data The data matrix used for clustering
#' @return BIC value (lower is better)
compute_kmeans_bic <- function(kmeans_obj, data) {
  n <- nrow(data)
  k <- length(unique(kmeans_obj$cluster))
  d <- ncol(data)
  
  # Within-cluster sum of squares
  wcss <- kmeans_obj$tot.withinss
  
  # Number of parameters: k cluster centers, each with d dimensions
  num_params <- k * d
  
  # BIC formula: n * log(WCSS/n) + num_params * log(n)
  bic <- n * log(wcss / n) + num_params * log(n)
  
  return(bic)
}

#' Compute AIC for k-means clustering
#' @param kmeans_obj A kmeans object
#' @param data The data matrix used for clustering
#' @return AIC value (lower is better)
compute_kmeans_aic <- function(kmeans_obj, data) {
  n <- nrow(data)
  k <- length(unique(kmeans_obj$cluster))
  d <- ncol(data)
  
  # Within-cluster sum of squares
  wcss <- kmeans_obj$tot.withinss
  
  # Number of parameters
  num_params <- k * d
  
  # AIC formula: n * log(WCSS/n) + 2 * num_params
  aic <- n * log(wcss / n) + 2 * num_params
  
  return(aic)
}

#' Compute BIC for GMM clustering
#' @param gmm_obj A Mclust object from mclust package
#' @return BIC value (higher is better for Mclust)
compute_gmm_bic <- function(gmm_obj) {
  # Mclust already computes BIC (note: higher is better in Mclust)
  return(gmm_obj$bic)
}

#' Compute AIC for GMM clustering
#' @param gmm_obj A Mclust object from mclust package
#' @return AIC value (lower is better)
compute_gmm_aic <- function(gmm_obj) {
  # AIC = 2k - 2ln(L)
  # where k is number of parameters and L is likelihood
  num_params <- gmm_obj$df
  loglik <- gmm_obj$loglik
  
  aic <- 2 * num_params - 2 * loglik
  
  return(aic)
}

#' Compute Gap Statistic for k-means clustering
#' @param data The data matrix
#' @param k Number of clusters to test
#' @param B Number of bootstrap samples
#' @param nstart Number of random starts for k-means
#' @return List with gap statistic and standard error
compute_gap_statistic <- function(data, k, B = 50, nstart = 10) {
  n <- nrow(data)
  d <- ncol(data)
  
  # Observed within-cluster dispersion
  kmeans_obs <- kmeans(data, centers = k, nstart = nstart)
  W_obs <- log(kmeans_obs$tot.withinss)
  
  # Generate reference distributions (uniform over range of data)
  W_refs <- numeric(B)
  for (b in 1:B) {
    # Generate uniform random data matching the range of original data
    ref_data <- matrix(0, nrow = n, ncol = d)
    for (j in 1:d) {
      ref_data[, j] <- runif(n, min = min(data[, j]), max = max(data[, j]))
    }
    
    # Cluster reference data
    kmeans_ref <- kmeans(ref_data, centers = k, nstart = nstart)
    W_refs[b] <- log(kmeans_ref$tot.withinss)
  }
  
  # Gap statistic
  gap <- mean(W_refs) - W_obs
  
  # Standard deviation and standard error
  sdk <- sqrt(mean((W_refs - mean(W_refs))^2))
  sk <- sdk * sqrt(1 + 1/B)
  
  return(list(gap = gap, se = sk))
}

#' Perform clustering using specified algorithm
#' @param embeddings The embedding matrix
#' @param k Number of clusters
#' @param method Clustering algorithm: 'kmeans' or 'gmm'
#' @param nstart Number of random starts (for k-means)
#' @param modelNames Model names for GMM (default: all models)
#' @return List with cluster assignments and clustering object
perform_clustering <- function(embeddings, k, 
                              method = c("kmeans", "gmm"),
                              nstart = 10,
                              modelNames = NULL) {
  
  method <- match.arg(method)
  
  if (method == "kmeans") {
    # K-means clustering
    clust_obj <- kmeans(embeddings, centers = k, nstart = nstart)
    return(list(
      cluster = clust_obj$cluster,
      object = clust_obj,
      method = "kmeans"
    ))
    
  } else if (method == "gmm") {
    # Gaussian Mixture Model clustering
    if (!requireNamespace("mclust", quietly = TRUE)) {
      stop("Package 'mclust' is required for GMM clustering. Install it with: install.packages('mclust')")
    }
    
    # If modelNames not specified, use default (allows all covariance structures)
    gmm_obj <- mclust::Mclust(embeddings, G = k, modelNames = modelNames, verbose = FALSE)
    
    return(list(
      cluster = gmm_obj$classification,
      object = gmm_obj,
      method = "gmm"
    ))
  }
}

#' Determine optimal number of clusters using specified method
#' @param embeddings The embedding matrix (classes x features)
#' @param labelsRY Class labels for each individual
#' @param classtable Table of class sizes
#' @param minsize Minimum cluster size
#' @param method Method for determining clusters: 'auto', 'bic', 'aic', 'gap'
#' @param max_k Maximum number of clusters to consider
#' @param nstart Number of random starts for k-means
#' @param gap_B Number of bootstrap samples for gap statistic
#' @param clustering_method Clustering algorithm: 'kmeans' or 'gmm'
#' @param gmm_modelNames Model names for GMM (NULL = automatic selection)
#' @return Optimal number of clusters
determine_optimal_clusters <- function(embeddings, 
                                      labelsRY, 
                                      classtable, 
                                      minsize,
                                      method = 'auto',
                                      max_k = NULL,
                                      nstart = 10,
                                      gap_B = 50,
                                      clustering_method = "kmeans",
                                      gmm_modelNames = NULL) {
  
  numlabel <- nrow(embeddings)
  
  # Set maximum k if not specified
  if (is.null(max_k)) {
    # Maximum k is limited by minimum cluster size requirement
    # Estimate based on total samples and minimum size
    total_samples <- sum(classtable)
    max_k <- min(numlabel, floor(total_samples / minsize))
  }
  
  # Ensure max_k is at least 2
  max_k <- max(2, min(max_k, numlabel))
  
  if (method == 'auto') {
    # Original auto method
    numcluster <- sum(classtable >= minsize) + any(classtable < minsize)
    
  } else if (method == 'bic') {
    # BIC-based selection
    bic_values <- numeric(max_k - 1)
    
    if (clustering_method == "kmeans") {
      for (k in 2:max_k) {
        clust_result <- perform_clustering(embeddings, k, method = "kmeans", nstart = nstart)
        bic_values[k - 1] <- compute_kmeans_bic(clust_result$object, embeddings)
      }
      # Select k with minimum BIC (lower is better)
      numcluster <- which.min(bic_values) + 1
      
    } else if (clustering_method == "gmm") {
      for (k in 2:max_k) {
        clust_result <- perform_clustering(embeddings, k, method = "gmm", 
                                          modelNames = gmm_modelNames)
        bic_values[k - 1] <- compute_gmm_bic(clust_result$object)
      }
      # Select k with maximum BIC (higher is better for Mclust)
      numcluster <- which.max(bic_values) + 1
    }
    
    cat("BIC values for k=2 to", max_k, "(", clustering_method, "):\n")
    print(data.frame(k = 2:max_k, BIC = bic_values))
    cat("Optimal k by BIC:", numcluster, "\n\n")
    
  } else if (method == 'aic') {
    # AIC-based selection
    aic_values <- numeric(max_k - 1)
    
    if (clustering_method == "kmeans") {
      for (k in 2:max_k) {
        clust_result <- perform_clustering(embeddings, k, method = "kmeans", nstart = nstart)
        aic_values[k - 1] <- compute_kmeans_aic(clust_result$object, embeddings)
      }
    } else if (clustering_method == "gmm") {
      for (k in 2:max_k) {
        clust_result <- perform_clustering(embeddings, k, method = "gmm",
                                          modelNames = gmm_modelNames)
        aic_values[k - 1] <- compute_gmm_aic(clust_result$object)
      }
    }
    
    # Select k with minimum AIC (lower is better)
    numcluster <- which.min(aic_values) + 1
    
    cat("AIC values for k=2 to", max_k, "(", clustering_method, "):\n")
    print(data.frame(k = 2:max_k, AIC = aic_values))
    cat("Optimal k by AIC:", numcluster, "\n\n")
    
  } else if (method == 'gap') {
    # Gap statistic only works with k-means
    if (clustering_method == "gmm") {
      warning("Gap statistic is only implemented for k-means. Switching to k-means for cluster selection.")
      clustering_method <- "kmeans"
    }
    
    gap_results <- data.frame(
      k = 2:max_k,
      gap = numeric(max_k - 1),
      se = numeric(max_k - 1)
    )
    
    for (k in 2:max_k) {
      gap_stat <- compute_gap_statistic(embeddings, k, B = gap_B, nstart = nstart)
      gap_results$gap[k - 1] <- gap_stat$gap
      gap_results$se[k - 1] <- gap_stat$se
    }
    
    # Find first k where Gap(k) >= Gap(k+1) - se(k+1)
    # This is the "1-SE" rule from Tibshirani et al.
    numcluster <- 2
    for (i in 1:(nrow(gap_results) - 1)) {
      if (gap_results$gap[i] >= gap_results$gap[i + 1] - gap_results$se[i + 1]) {
        numcluster <- gap_results$k[i]
        break
      }
    }
    
    cat("Gap statistic results:\n")
    print(gap_results)
    cat("Optimal k by Gap statistic:", numcluster, "\n\n")
    
  } else {
    stop("Unknown method. Choose from: 'auto', 'bic', 'aic', 'gap'")
  }
  
  # Ensure numcluster is at least 1 and at most numlabel
  numcluster <- max(1, min(numcluster, numlabel))
  
  # Validate that clusters will have sufficient size
  # This is a rough check - actual cluster sizes may vary
  avg_cluster_size <- sum(classtable) / numcluster
  if (avg_cluster_size < minsize) {
    warning(paste("Optimal k =", numcluster, "may result in clusters smaller than minsize =", 
                  minsize, ". Consider using a smaller k or reducing minsize."))
  }
  
  return(numcluster)
}


# -----------------------------------------------------------------------------
# MAIN FUNCTION: Surrogate-assisted Clustered Conformal Inference
# -----------------------------------------------------------------------------

#' Surrogate-assisted Clustered Conformal Inference for Efficient Individual 
#' Causal Effect Estimation
#'
#' Enhanced version with multiple cluster selection methods (BIC, AIC, Gap Statistic)
#' and GMM clustering option
#'
#' @param df Data frame containing the variables
#' @param train.idx Training data indices
#' @param eval.idx Evaluation data indices
#' @param numcluster Number of clusters or method: 'auto', 'bic', 'aic', 'gap', or numeric value
#' @param minsize Minimum size of each cluster if numcluster = 'auto'
#' @param qvec Selected quantiles for the non-conformity score embeddings
#' @param cluster_max_k Maximum number of clusters to consider for optimization methods
#' @param cluster_nstart Number of random starts for k-means
#' @param gap_B Number of bootstrap samples for gap statistic
#' @param clustering_method Clustering algorithm: 'kmeans' or 'gmm' (Gaussian Mixture Model)
#' @param gmm_modelNames GMM covariance structure: NULL (auto), 'EII', 'VII', 'EEI', 'VEI', 'EVI', 'VVI', 'EEE', 'EVE', 'VEE', 'VVE', 'EEV', 'VEV', 'EVV', 'VVV'
#' @param outcome.type Type of outcome: "Continuous" or "Categorical"
#' @param SL.library SuperLearner library for fitting models
#' @param alphaCI Significance level (default 0.05)
#' @param nested Whether to use nested conformal prediction
#'
#' @return Data frame with prediction intervals/sets for evaluation samples
#' 
#' @details
#' Clustering methods:
#' * kmeans: Traditional k-means (spherical clusters, hard assignments)
#' * gmm: Gaussian Mixture Model (elliptical clusters, soft assignments, more flexible)
#' 
#' GMM modelNames (covariance structures):
#' * NULL: Automatically select best model
#' * 'EII': Spherical, equal volume
#' * 'VII': Spherical, unequal volume
#' * 'EEE': Ellipsoidal, equal volume/shape/orientation
#' * 'VVV': Ellipsoidal, varying volume/shape/orientation (most flexible)
#'
#' @examples
#' # Use k-means with BIC
#' result <- SurrClusterConformalDR(
#'   df = mydata,
#'   train.idx = train_indices,
#'   eval.idx = eval_indices,
#'   numcluster = 'bic',
#'   clustering_method = 'kmeans',
#'   outcome.type = "Continuous"
#' )
#'
#' # Use GMM with automatic model selection
#' result <- SurrClusterConformalDR(
#'   df = mydata,
#'   train.idx = train_indices,
#'   eval.idx = eval_indices,
#'   numcluster = 'bic',
#'   clustering_method = 'gmm',
#'   gmm_modelNames = NULL,
#'   outcome.type = "Continuous"
#' )
#'
#' # Use GMM with flexible covariance structure
#' result <- SurrClusterConformalDR(
#'   df = mydata,
#'   train.idx = train_indices,
#'   eval.idx = eval_indices,
#'   numcluster = 'aic',
#'   clustering_method = 'gmm',
#'   gmm_modelNames = 'VVV',  # Most flexible
#'   outcome.type = "Continuous"
#' )
#'
#' @export
SurrClusterConformalDR <- function(df,
                                  train.idx, eval.idx,
                                  numcluster = 'auto',
                                  minsize = 150,
                                  qvec = c(0.5, 0.6, 0.7, 0.8, 0.9),
                                  cluster_max_k = NULL,
                                  cluster_nstart = 10,
                                  gap_B = 50,
                                  clustering_method = c("kmeans", "gmm"),
                                  gmm_modelNames = NULL,
                                  outcome.type = c("Continuous", "Categorical"),
                                  SL.library = c("SL.glm"),
                                  alphaCI = 0.05,
                                  nested = TRUE) {
  
  # Match arguments
  clustering_method <- match.arg(clustering_method)
  outcome.type <- match.arg(outcome.type)
  
  # Begin estimation
  N <- nrow(df)
  
  ## -------------------
  # Create non-conformity score
  ## -------------------
  
  if (outcome.type == "Continuous") {
    # CQR by quantile regression function
    qrf.Y1.obj <- grf::quantile_forest(
      X = df[train.idx, ] %>% filter(A == 1 & D == 1) %>%
        dplyr::select(grep("^([XS])", colnames(df), value = TRUE)),
      Y = df[train.idx, ] %>% filter(A == 1 & D == 1) %>% pull(Y),
      quantiles = c(alphaCI / 2, 1 - alphaCI / 2),
      honesty = FALSE
    )
    
    # Change for SuperLearner model
    df_train_A1 <- data.frame(
      tau = predict(qrf.Y1.obj,
                    newdata = df[train.idx, grep("^([XS])", colnames(df), value = TRUE)])$predictions,
      df[train.idx, grep("^([XA])", colnames(df), value = TRUE)]
    ) %>%
      filter(A == 1) %>%
      dplyr::select(-A)
    
    m.X1q1.obj <- SuperLearner(
      Y = df_train_A1$tau.1,
      X = df_train_A1[, grep("([X])", colnames(df_train_A1), value = TRUE), drop = FALSE],
      SL.library = SL.library
    )
    
    m.X1q2.obj <- SuperLearner(
      Y = df_train_A1$tau.2,
      X = df_train_A1[, grep("([X])", colnames(df_train_A1), value = TRUE), drop = FALSE],
      SL.library = SL.library
    )
    
    q.Y1 <- cbind(
      predict(m.X1q1.obj, newdata = df[, grep("([X])", colnames(df), value = TRUE), drop = FALSE])$pred,
      predict(m.X1q2.obj, newdata = df[, grep("([X])", colnames(df), value = TRUE), drop = FALSE])$pred
    )
    
    qrf.Y0.obj <- grf::quantile_forest(
      X = df[train.idx, ] %>% filter(A == 0 & D == 1) %>%
        dplyr::select(grep("^([XS])", colnames(df), value = TRUE)),
      Y = df[train.idx, ] %>% filter(A == 0 & D == 1) %>% pull(Y),
      quantiles = c(alphaCI / 2, 1 - alphaCI / 2),
      honesty = FALSE
    )
    
    # Change for SuperLearner model
    df_train_A0 <- data.frame(
      tau = predict(qrf.Y0.obj,
                    newdata = df[train.idx, grep("^([XS])", colnames(df), value = TRUE)])$predictions,
      df[train.idx, grep("^([XA])", colnames(df), value = TRUE)]
    ) %>%
      filter(A == 0) %>%
      dplyr::select(-A)
    
    m.X0q1.obj <- SuperLearner(
      Y = df_train_A0$tau.1,
      X = df_train_A0[, grep("([X])", colnames(df_train_A0), value = TRUE), drop = FALSE],
      SL.library = SL.library
    )
    
    m.X0q2.obj <- SuperLearner(
      Y = df_train_A0$tau.2,
      X = df_train_A0[, grep("([X])", colnames(df_train_A0), value = TRUE), drop = FALSE],
      SL.library = SL.library
    )
    
    q.Y0 <- cbind(
      predict(m.X0q1.obj, newdata = df[, grep("([X])", colnames(df), value = TRUE), drop = FALSE])$pred,
      predict(m.X0q2.obj, newdata = df[, grep("([X])", colnames(df), value = TRUE), drop = FALSE])$pred
    )
    
    df$R.wS <- with(df, pmax(q.Y1[, 1] - Y, Y - q.Y1[, 2]) * A) +
      with(df, pmax(q.Y0[, 1] - Y, Y - q.Y0[, 2]) * (1 - A))
  }
  
  if (outcome.type == "Categorical") {
    # wS
    objY1.XS <- nnet::multinom(
      paste("Y~", paste(grep("^([XS])", colnames(df), value = TRUE), collapse = "+")),
      data = df[train.idx, ],
      subset = A == 1 & D == 1, trace = FALSE
    )
    probY1.XS <- predict(objY1.XS,
                         newdata = df[, grep("^([XS])", colnames(df), value = TRUE)],
                         type = "probs")
    
    objY0.XS <- nnet::multinom(
      paste("Y~", paste(grep("^([XS])", colnames(df), value = TRUE), collapse = "+")),
      data = df[train.idx, ],
      subset = A == 0 & D == 1, trace = FALSE
    )
    probY0.XS <- predict(objY0.XS,
                         newdata = df[, grep("^([XS])", colnames(df), value = TRUE)],
                         type = "probs")
    
    # If a single vector is returned, bind it to provide probabilities for both classes
    if (is.null(dim(probY0.XS))) {
      probY0.XS <- cbind(probY0.XS, 1 - probY0.XS)
      colnames(probY0.XS) <- objY0.XS$lev
    }
    if (is.null(dim(probY1.XS))) {
      probY1.XS <- cbind(probY1.XS, 1 - probY1.XS)
      colnames(probY1.XS) <- objY1.XS$lev
    }
    
    # Check for missing factor level and impute with zero
    if (ncol(probY1.XS) < length(unique(df$Y))) {
      level.prob <- colnames(probY1.XS)
      level.NA <- setdiff(unique(df$Y), as.numeric(level.prob))
      probY1.XS <- cbind(probY1.XS, matrix(0, nrow = nrow(probY1.XS), ncol = length(level.NA)))
      colnames(probY1.XS) <- c(level.prob, level.NA)
    }
    
    if (ncol(probY0.XS) < length(unique(df$Y))) {
      level.prob <- colnames(probY0.XS)
      level.NA <- setdiff(unique(df$Y), as.numeric(level.prob))
      probY0.XS <- cbind(probY0.XS, matrix(0, nrow = nrow(probY0.XS), ncol = length(level.NA)))
      colnames(probY0.XS) <- c(level.prob, level.NA)
    }
    
    # Compute the non-conformity score of each outcome
    R.XSY1Mat.XS <- apply(probY1.XS, 1, function(x) {
      1 - cumsum(sort(x))[rank(x)]
    }) %>% t()
    R.XSY0Mat.XS <- apply(probY0.XS, 1, function(x) {
      1 - cumsum(sort(x))[rank(x)]
    }) %>% t()
    
    # Choose the observed one based on the outcomes
    df$R.wS <- with(df, mapply(function(row_index, col_index) R.XSY1Mat.XS[row_index, col_index],
                               row_index = 1:nrow(R.XSY1Mat.XS),
                               col_index = match(df$Y, as.numeric(colnames(probY1.XS)))) * A) +
      with(df, mapply(function(row_index, col_index) R.XSY0Mat.XS[row_index, col_index],
                      row_index = 1:nrow(R.XSY0Mat.XS),
                      col_index = match(df$Y, as.numeric(colnames(probY0.XS)))) * (1 - A))
  }
  
  ## -------------------------
  # Model fitting for nuisance functions
  ## -------------------------
  
  df.train <- df[train.idx, ]
  df.eval <- df[eval.idx, ]
  
  model.obj.wS <- list()
  
  # Model training for D|X, A
  psD.obj <- model.obj.wS$psD.obj <- SuperLearner(
    Y = df.train$D,
    X = df.train[, grep("([XA])", colnames(df.train), value = TRUE), drop = FALSE],
    SL.library = SL.library,
    family = binomial()
  )
  
  # Model for A|X
  psA.obj <- model.obj.wS$psA.obj <- SuperLearner(
    Y = as.numeric(df.train$A),
    X = df.train[, grep("([X])", colnames(df.train), value = TRUE), drop = FALSE],
    SL.library = SL.library,
    family = binomial()
  )
  
  ## Further split df.train into two folds, I_11 and I_12
  ## One for computing the initial estimate, the other for training the outcome model
  folds.train.idx <- caret::createFolds(1:length(train.idx), k = 2)
  train1.idx <- folds.train.idx[[1]]
  train2.idx <- folds.train.idx[[2]]
  df.train1 <- df.train[train1.idx, ]
  df.train2 <- df.train[train2.idx, ]
  
  # Initial estimate for theta by IPW (with first split data I_11)
  thetaA1.init <- initialize.theta(df.train1, psA.obj, psD.obj,
                                   target.A = 1, alphaCI = alphaCI)
  thetaA0.init <- initialize.theta(df.train1, psA.obj, psD.obj,
                                   target.A = 0, alphaCI = alphaCI)
  
  # Fit the model on the second split data I_12
  ## Model training with S
  model.obj.wS$mA1.obj <- SuperLearner(
    Y = as.numeric(df.train2$R.wS <= thetaA1.init)[df.train2$D == 1],
    X = df.train2[df.train2$D == 1, grep("^([XS])", colnames(df.train2), value = TRUE), drop = FALSE],
    SL.library = SL.library,
    family = binomial()
  )
  
  model.obj.wS$mA0.obj <- SuperLearner(
    Y = as.numeric(df.train2$R.wS <= thetaA0.init)[df.train2$D == 1],
    X = df.train2[df.train2$D == 1, grep("^([XS])", colnames(df.train2), value = TRUE), drop = FALSE],
    SL.library = SL.library,
    family = binomial()
  )
  
  ## -------------------------
  # ENHANCED CLUSTERING SECTION WITH GMM OPTION
  ## -------------------------
  
  # The class label for the evaluation data
  df.eval <- df.eval %>%
    mutate(class = factor(R) %>% as.numeric())
  
  # Merge rare class to a single class
  ## Merge class with size less than minsize
  classtable <- table(df.eval$class)
  rareclass <- names(classtable)[classtable < minsize]
  df.eval$class <- sapply(df.eval$class %>% as.character,
                          function(cls){
                            ifelse(cls %in% rareclass, 'Merged', cls)
                          }) %>% factor() %>% as.numeric()
  labelsRY <- df.eval$class
  numlabel <- length(unique(labelsRY))
  
  # Quantile embeddings for the training data
  embeddings  <- matrix(0, nrow = numlabel, ncol = length(qvec))
  scores <- df.eval$R.wS
  for (i in 1:numlabel) {
    if(any(labelsRY == i, na.rm = TRUE)){
      class_i_scores <- scores[which(labelsRY == i)]
      embeddings[i, ] <- quantile(class_i_scores, q = qvec, na.rm = TRUE)
    }
  }
  # A robust check for NA
  embeddings[is.na(embeddings)] <- 0
  
  # REVISED: Determine number of clusters using specified method
  numbers_only <- function(x) !grepl("\\D", x)
  
  if(numbers_only(numcluster)){
    # User specified numeric value
    numcluster <- as.numeric(numcluster)
  } else {
    # Use optimization method (auto, bic, aic, or gap)
    numcluster <- determine_optimal_clusters(
      embeddings = embeddings,
      labelsRY = labelsRY,
      classtable = classtable,
      minsize = minsize,
      method = numcluster,
      max_k = cluster_max_k,
      nstart = cluster_nstart,
      gap_B = gap_B,
      clustering_method = clustering_method,
      gmm_modelNames = gmm_modelNames
    )
  }
  
  # Perform clustering using selected method
  if(numcluster < nrow(embeddings)){
    clustering_result <- perform_clustering(
      embeddings = embeddings,
      k = numcluster,
      method = clustering_method,
      nstart = cluster_nstart,
      modelNames = gmm_modelNames
    )
    
    # Remap the cluster to each individual from the evaluation set
    df.eval$cluster <- clustering_result$cluster[df.eval$class]
    
    # Store clustering object for potential diagnostics
    clustering_obj <- clustering_result$object
    
  } else {
    # If numcluster >= nrow(embeddings), keep original class structure
    clustering_result <- list(cluster = 1:nrow(embeddings))
    df.eval$cluster <- clustering_result$cluster[df.eval$class]
    clustering_obj <- NULL
  }
  
  ## -------------------------
  # Construct the estimated cutoff values for the observed outcomes
  ## -------------------------
  
  ## EIF-based with surrogates
  cutoff.Y1.wS_obs <- by(df.eval, list(cluster = df.eval$cluster),
                         function(dfCluster){
                           tryCatch(
                             uniroot(
                               f = function(theta) {
                                 eif_theta(theta,
                                           df.eval = dfCluster,
                                           model.obj = model.obj.wS,
                                           target.A = 1, wS = TRUE,
                                           alphaCI = alphaCI,
                                           counterfactual = FALSE)
                               },
                               interval = quantile(dfCluster$R.wS, c(0.01, 0.99), na.rm = TRUE),
                               extendInt = "yes",
                               maxiter = 100
                             )$root,
                             error = function(e) {
                               quantile(dfCluster$R.wS, 1 - alphaCI, na.rm = TRUE)
                             }
                           )
                         }) %>% array2DF() %>%
    arrange(cluster) %>% pull(Value)
  
  cutoff.Y0.wS_obs <- by(df.eval, list(cluster = df.eval$cluster),
                         function(dfCluster){
                           tryCatch(
                             uniroot(
                               f = function(theta) {
                                 eif_theta(theta,
                                           df.eval = dfCluster,
                                           model.obj = model.obj.wS,
                                           target.A = 0, wS = TRUE,
                                           alphaCI = alphaCI,
                                           counterfactual = FALSE)
                               },
                               interval = quantile(dfCluster$R.wS, c(0.01, 0.99), na.rm = TRUE),
                               extendInt = "yes",
                               maxiter = 100
                             )$root,
                             error = function(e) {
                               quantile(dfCluster$R.wS, 1 - alphaCI, na.rm = TRUE)
                             }
                           )
                         }) %>% array2DF() %>%
    arrange(cluster) %>% pull(Value)
  
  ## -------------------------
  # Construct the prediction sets
  ## -------------------------
  
  if (outcome.type == "Continuous") {
    # Conformal inference on the observed outcomes
    # wS
    cls <- clustering_result$cluster[df.eval$class]
    lower.Y1_A1.wS <- c(predict(m.X1q1.obj,
                                newdata = df.eval[, grep("([XA])", colnames(df.eval), value = TRUE)])$pred) -
      cutoff.Y1.wS_obs[cls]
    
    upper.Y1_A1.wS <- c(predict(m.X1q2.obj,
                                newdata = df.eval[, grep("([XA])", colnames(df.eval), value = TRUE)])$pred) +
      cutoff.Y1.wS_obs[cls]
    
    lower.Y0_A0.wS <- c(predict(m.X0q1.obj,
                                newdata = df.eval[, grep("([XA])", colnames(df.eval), value = TRUE)])$pred) -
      cutoff.Y0.wS_obs[cls]
    
    upper.Y0_A0.wS <- c(predict(m.X0q2.obj,
                                newdata = df.eval[, grep("([XA])", colnames(df.eval), value = TRUE)])$pred) +
      cutoff.Y0.wS_obs[cls]
    
    ## For ITE
    # Construct the estimated cutoff values for the counterfactual outcomes
    cutoff.Y1.wS <- by(df.eval, list(cluster = df.eval$cluster),
                       function(dfCluster){
                         tryCatch(
                           uniroot(
                             f = function(theta) {
                               eif_theta(theta,
                                         df.eval = dfCluster,
                                         model.obj = model.obj.wS,
                                         target.A = 1, wS = TRUE,
                                         alphaCI = alphaCI)
                             },
                             interval = quantile(dfCluster$R.wS, c(0.01, 0.99), na.rm = TRUE),
                             extendInt = "yes",
                             maxiter = 100
                           )$root,
                           error = function(e) {
                             quantile(dfCluster$R.wS, 1 - alphaCI, na.rm = TRUE)
                           }
                         )
                       }) %>% array2DF() %>%
      arrange(cluster) %>% pull(Value)
    
    cutoff.Y0.wS <- by(df.eval, list(cluster = df.eval$cluster),
                       function(dfCluster){
                         tryCatch(
                           uniroot(
                             f = function(theta) {
                               eif_theta(theta,
                                         df.eval = dfCluster,
                                         model.obj = model.obj.wS,
                                         target.A = 0, wS = TRUE,
                                         alphaCI = alphaCI)
                             },
                             interval = quantile(dfCluster$R.wS, c(0.01, 0.99), na.rm = TRUE),
                             extendInt = "yes",
                             maxiter = 100
                           )$root,
                           error = function(e) {
                             quantile(dfCluster$R.wS, 1 - alphaCI, na.rm = TRUE)
                           }
                         )
                       }) %>% array2DF() %>%
      arrange(cluster) %>% pull(Value)
    
    # wS
    lower.Y1_A0.wS <- c(predict(m.X1q1.obj,
                                newdata = df.eval[, grep("([XA])", colnames(df.eval), value = TRUE)])$pred) -
      cutoff.Y1.wS[cls]
    
    upper.Y1_A0.wS <- c(predict(m.X1q2.obj,
                                newdata = df.eval[, grep("([XA])", colnames(df.eval), value = TRUE)])$pred) +
      cutoff.Y1.wS[cls]
    
    lower.Y0_A1.wS <- c(predict(m.X0q1.obj,
                                newdata = df.eval[, grep("([XA])", colnames(df.eval), value = TRUE)])$pred) -
      cutoff.Y0.wS[cls]
    
    upper.Y0_A1.wS <- c(predict(m.X0q2.obj,
                                newdata = df.eval[, grep("([XA])", colnames(df.eval), value = TRUE)])$pred) +
      cutoff.Y0.wS[cls]
    
    # Inference on the observed outcomes
    lower.Y.wS <- lower.Y1_A1.wS * (df.eval$A == 1) +
      lower.Y0_A0.wS * (df.eval$A == 0)
    upper.Y.wS <- upper.Y1_A1.wS * (df.eval$A == 1) +
      upper.Y0_A0.wS * (df.eval$A == 0)
    
    # Inference on tau when A = 0 + A = 1
    lower.tau.wS <- (lower.Y1_A0.wS - (df.eval$Y)) * (df.eval$A == 0) +
      ((df.eval$Y) - upper.Y0_A1.wS) * (df.eval$A == 1)
    upper.tau.wS <- (upper.Y1_A0.wS - (df.eval$Y)) * (df.eval$A == 0) +
      ((df.eval$Y) - lower.Y0_A1.wS) * (df.eval$A == 1)
    
    # Organize the results
    df.eval <- cbind(df.eval,
                     lower.Y = lower.Y.wS,
                     upper.Y = upper.Y.wS,
                     lower.tau = lower.tau.wS,
                     upper.tau = upper.tau.wS)
  }
  
  if (outcome.type == "Categorical") {
    
    R.XSY1Mat.XS.eval <- R.XSY1Mat.XS[eval.idx, ]
    R.XSY0Mat.XS.eval <- R.XSY0Mat.XS[eval.idx, ]
    
    # For A = 0, Y1
    Y1predSet.A1wS <- sapply(1:nrow(R.XSY1Mat.XS.eval), function(idx){
      x <- R.XSY1Mat.XS.eval[idx, ]
      cls <- clustering_result$cluster[df.eval$class[idx]]
      level.Y <- colnames(probY1.XS)[order(x)]
      level.Y[which(x[order(x)] <= cutoff.Y1.wS_obs[cls])]
    }, simplify = FALSE)
    
    Y0predSet.A0wS <- sapply(1:nrow(R.XSY0Mat.XS.eval), function(idx){
      x <- R.XSY0Mat.XS.eval[idx, ]
      cls <- clustering_result$cluster[df.eval$class[idx]]
      level.Y <- colnames(probY0.XS)[order(x)]
      level.Y[which(x[order(x)] <= cutoff.Y0.wS_obs[cls])]
    }, simplify = FALSE)
    
    # Organize the output for the evaluation sets
    ## For the observed one
    ## With surrogates
    df.eval$sets_observed <- Y0predSet.A0wS
    df.eval$sets_observed[df.eval$A == 1] <- Y1predSet.A1wS[df.eval$A == 1]
  }
  
  # Return the evaluation df without non-conformal scores
  return(df.eval %>% dplyr::rename(scores = R.wS))
}
