#' Surrogate-assisted Clustered Conformal Inference for Efficient Individual Causal Effect
#' Estimation
#'
#' Surrogate-assisted Clustered Conformal Inference for Efficiency and Fairness
#'
#' @param numcluster The number of cluster for the non-conformity scores. The default is 'auto'.
#' @param minsize The minimum size of each cluster if numcluster = 'auto'.
#' @param qvec The selected quantiles for the non-conformity score embeddings.
#'
#'
#' See details of other parameters in [SurrConformalDR::SurrConformalDR()].
#'
#' @return
#' when \code{outcome.type = "Continuous"}:
#' * lower.Y, upper Y: prediction regions for the observed outcomes
#' with index `eval.idx`
#' * lower.tau, upper.tau: prediction regions for the individualized treatment effects
#' with index `eval.idx`
#'
#' when \code{outcome.type = "Categorical"}:
#' * sets_observed: prediction sets for the observed outcomes
#' with index `eval.idx`
#' @export
SurrClusterConformalDR <- function(df,
                            train.idx, eval.idx,
                            numcluster = 'auto', # number of cluster
                            minsize = 150,  # minimum size of class
                            qvec = c(0.5, 0.6, 0.7, 0.8, 0.9), # quantiles for embeddings
                            outcome.type = c("Continuous", "Categorical"), # Categorical
                            SL.library = c("SL.glm"),
                            alphaCI = 0.05,
                            nested = TRUE) {
  # begin estimation
  N <- nrow(df)
  outcome.type <- match.arg(outcome.type)
  ## -------------------
  # create non-conformity score
  if (outcome.type == "Continuous") {
    # CQR by quantile regression function
    qrf.Y1.obj <- grf::quantile_forest(
      X = df[train.idx, ] %>% filter(A == 1 &
                                       D == 1) %>%
        dplyr::select(grep("^([XS])",
                           colnames(df),
                           value = TRUE
        )),
      Y = df[train.idx, ] %>% filter(A == 1 &
                                       D == 1) %>%
        pull(Y),
      quantiles = c(alphaCI / 2, 1 - alphaCI / 2),
      honesty = FALSE
    )

    # change for SuperLearner model
    df_train_A1 <- data.frame(
      tau = predict(qrf.Y1.obj,
                    newdata = df[train.idx, grep("^([XS])",
                                                 colnames(df),
                                                 value = TRUE
                    )]
      )$predictions,
      df[train.idx, grep("^([XA])",
                         colnames(df),
                         value = TRUE
      )]
    ) %>%
      filter(A == 1) %>%
      dplyr::select(-A)

    m.X1q1.obj <- SuperLearner(
      Y = df_train_A1$tau.1,
      X = df_train_A1[, grep("([X])",
                             colnames(df_train_A1),
                             value = TRUE
      ),
      drop = FALSE
      ],
      SL.library = SL.library
    )



    m.X1q2.obj <- SuperLearner(
      Y = df_train_A1$tau.2,
      X = df_train_A1[, grep("([X])",
                             colnames(df_train_A1),
                             value = TRUE
      ),
      drop = FALSE
      ],
      SL.library = SL.library
    )

#Produces out-of-sample lower/upper quantile predictions for all rows in df (still representing the A=1 world).
    q.Y1 <- cbind(
      predict(m.X1q1.obj, newdata = df[, grep("([X])", colnames(df),
                                              value = TRUE
      ),
      drop = FALSE
      ])$pred,
      predict(m.X1q2.obj, newdata = df[, grep("([X])", colnames(df),
                                              value = TRUE
      ),
      drop = FALSE
      ])$pred
    )

    qrf.Y0.obj <- grf::quantile_forest(
      X = df[train.idx, ] %>% filter(A == 0 &
                                       D == 1) %>%
        dplyr::select(grep("^([XS])",
                           colnames(df),
                           value = TRUE
        )),
      Y = df[train.idx, ] %>% filter(A == 0 &
                                       D == 1) %>%
        pull(Y),
      quantiles = c(alphaCI / 2, 1 - alphaCI / 2),
      honesty = FALSE
    )

    # change for SuperLearner model
    df_train_A0 <- data.frame(
      tau = predict(qrf.Y0.obj,
                    newdata = df[train.idx, grep("^([XS])",
                                                 colnames(df),
                                                 value = TRUE
                    )]
      )$predictions,
      df[train.idx, grep("^([XA])",
                         colnames(df),
                         value = TRUE
      )]
    ) %>%
      filter(A == 0) %>%
      dplyr::select(-A)

    m.X0q1.obj <- SuperLearner(
      Y = df_train_A0$tau.1,
      X = df_train_A0[, grep("([X])",
                             colnames(df_train_A0),
                             value = TRUE
      ),
      drop = FALSE
      ],
      SL.library = SL.library
    )



    m.X0q2.obj <- SuperLearner(
      Y = df_train_A0$tau.2,
      X = df_train_A0[, grep("([X])",
                             colnames(df_train_A0),
                             value = TRUE
      ),
      drop = FALSE
      ],
      SL.library = SL.library
    )

    q.Y0 <- cbind(
      predict(m.X0q1.obj, newdata = df[, grep("([X])", colnames(df),
                                              value = TRUE
      ),
      drop = FALSE
      ])$pred,
      predict(m.X0q2.obj, newdata = df[, grep("([X])", colnames(df),
                                              value = TRUE
      ),
      drop = FALSE
      ])$pred
    )
#Build the non-conformity score
    df$R.wS <- with(df, pmax(q.Y1[, 1] - Y, Y - q.Y1[, 2]) * A) +
      with(df, pmax(q.Y0[, 1] - Y, Y - q.Y0[, 2]) * (1 - A))
  }

  if (outcome.type == "Categorical") {
    # wS
    objY1.XS <- nnet::multinom(
      paste("Y~", paste(grep("^([XS])",
                             colnames(df),
                             value = TRUE
      ), collapse = "+")),
      data = df[train.idx, ],
      subset = A == 1 & D == 1, trace = FALSE
    )
    probY1.XS <- predict(objY1.XS,
                         newdata = df[, grep("^([XS])",
                                             colnames(df),
                                             value = TRUE
                         )],
                         type = "probs"
    )

    objY0.XS <- nnet::multinom(
      paste("Y~", paste(grep("^([XS])",
                             colnames(df),
                             value = TRUE
      ), collapse = "+")),
      data = df[train.idx, ],
      subset = A == 0 & D == 1, trace = FALSE
    )
    probY0.XS <- predict(objY0.XS,
                         newdata = df[, grep("^([XS])",
                                             colnames(df),
                                             value = TRUE
                         )],
                         type = "probs"
    )
    # If a single vector is returned, it means the probabilities for one class are being returned
    # Bind it to provide probabilities for both classes
    if (is.null(dim(probY0.XS))) {
      probY0.XS <- cbind(probY0.XS, 1 - probY0.XS)
      colnames(probY0.XS) <- objY0.XS$lev
    }
    if (is.null(dim(probY1.XS))) {
      probY1.XS <- cbind(probY1.XS, 1 - probY1.XS)
      colnames(probY1.XS) <- objY1.XS$lev
    }
    # check the for missing factor level and impute with zero
    if (ncol(probY1.XS) < length(unique(df$Y))) {
      level.prob <- colnames(probY1.XS)
      level.NA <- setdiff(
        unique(df$Y),
        as.numeric(level.prob)
      )
      probY1.XS <- cbind(
        probY1.XS,
        matrix(0,
               nrow = nrow(probY1.XS),
               ncol = length(level.NA)
        )
      )
      colnames(probY1.XS) <- c(level.prob, level.NA)
    }

    if (ncol(probY0.XS) < length(unique(df$Y))) {
      level.prob <- colnames(probY0.XS)
      level.NA <- setdiff(
        unique(df$Y),
        as.numeric(level.prob)
      )
      probY0.XS <- cbind(
        probY0.XS,
        matrix(0,
               nrow = nrow(probY0.XS),
               ncol = length(level.NA)
        )
      )
      colnames(probY0.XS) <- c(level.prob, level.NA)
    }

    # compute the non-conformity score of each outcomes
    R.XSY1Mat.XS <- apply(probY1.XS, 1, function(x) {
      1 - cumsum(sort(x))[rank(x)]
    }) %>% t()
    R.XSY0Mat.XS <- apply(probY0.XS, 1, function(x) {
      1 - cumsum(sort(x))[rank(x)]
    }) %>% t()


    # choose the observed one based on the outcomes
    df$R.wS <- with(df, mapply(function(row_index, col_index) R.XSY1Mat.XS[row_index, col_index],
                               row_index = 1:nrow(R.XSY1Mat.XS),
                               col_index = match(
                                 df$Y,
                                 as.numeric(colnames(probY1.XS))
                               )
    ) * A) +
      with(df, mapply(function(row_index, col_index) R.XSY0Mat.XS[row_index, col_index],
                      row_index = 1:nrow(R.XSY0Mat.XS),
                      col_index = match(
                        df$Y,
                        as.numeric(colnames(probY0.XS))
                      )
      ) * (1 - A))
  }

  # -------------------------
  # model fitting for nuisance functions
  df.train <- df[train.idx, ]
  df.eval <- df[eval.idx, ]

  model.obj.wS <- list()

  # model training for D|X, A
  psD.obj <-
    model.obj.wS$psD.obj <-
    SuperLearner(
      Y = df.train$D,
      X = df.train[, grep("([XA])",
                          colnames(df.train),
                          value = TRUE
      ),
      drop = FALSE
      ],
      SL.library = SL.library,
      family = binomial()
    )
  # model for A|X
  psA.obj <-
    model.obj.wS$psA.obj <-
    SuperLearner(
      Y = as.numeric(df.train$A),
      X = df.train[, grep("([X])",
                          colnames(df.train),
                          value = TRUE
      ),
      drop = FALSE
      ],
      SL.library = SL.library,
      family = binomial()
    )


  ## further split df.train into two folds, I_11 and I_12
  ## one for computing the initial estimate,
  ## the other one is for training the outcome model
  folds.train.idx <- caret::createFolds(1:length(train.idx), k = 2)
  train1.idx <- folds.train.idx[[1]]
  train2.idx <- folds.train.idx[[2]]
  df.train1 <- df.train[train1.idx, ]
  df.train2 <- df.train[train2.idx, ]

  # initial estimate for theta by IPW (with first split data I_11)
  thetaA1.init <- initialize.theta(df.train1, psA.obj, psD.obj,
                                   target.A = 1, alphaCI = alphaCI)
  thetaA0.init <- initialize.theta(df.train1, psA.obj, psD.obj,
                                   target.A = 0, alphaCI = alphaCI)
  # fit the model on the second split data I_12
  ## model training with S
  model.obj.wS$mA1.obj <- SuperLearner(
    Y = as.numeric(df.train2$R.wS <= thetaA1.init)[df.train2$D == 1],
    X = df.train2[df.train2$D == 1, grep("^([XS])",
                                         colnames(df.train2),
                                         value = TRUE
    ),
    drop = FALSE
    ],
    SL.library = SL.library,
    family = binomial()
  )

  model.obj.wS$mA0.obj <- SuperLearner(
    Y = as.numeric(df.train2$R.wS <= thetaA0.init)[df.train2$D == 1],
    X = df.train2[df.train2$D == 1, grep("^([XS])",
                                         colnames(df.train2),
                                         value = TRUE
    ),
    drop = FALSE
    ],
    SL.library = SL.library,
    family = binomial()
  )
  # the class label for the evaluation data
  df.eval <- df.eval %>%
    mutate(class = factor(R) %>% as.numeric())
  # merge rare class to a single class
  ## merge class with size less than minsize
  classtable <- table(df.eval$class)
  rareclass <- names(classtable)[classtable < minsize]
  df.eval$class <- sapply(df.eval$class %>% as.character,
                           function(cls){
    ifelse(cls %in% rareclass, 'Merged', cls)
  }) %>% factor() %>% as.numeric()
  labelsRY <- df.eval$class
  numlabel <- length(unique(labelsRY))

  # quantile embeddings for the training data
  embeddings  <- matrix(0, nrow = numlabel,
                        ncol = length(qvec))
  scores <- df.eval$R.wS
  for (i in 1:numlabel) {
    if(any(labelsRY== i, na.rm = TRUE)){
      class_i_scores <- scores[which(labelsRY == i)]
      embeddings[i, ] <- quantile(class_i_scores,
                                  q = qvec, na.rm = TRUE)
    }
  }
  # a robust check for NA
  embeddings[is.na(embeddings)] <- 0
  numbers_only <- function(x) !grepl("\\D", x)
  # the number of clusters
  if(!numbers_only(numcluster)){
    # group the goups with less minsize samples
    numcluster <- sum(classtable >= minsize) +
      any(classtable < minsize)
    # cat(numcluster)
    # numcluster <-
    #   min(min(min(classtable)), 1/alphaCI - 1)
  }else{
    numcluster <- as.numeric(numcluster)
  }
  if(numcluster < nrow(embeddings)){
    # Compute distance matrix
    dist_matrix <- dist(embeddings, method = "euclidean")
    
    # Perform hierarchical clustering
    hclust_result <- hclust(dist_matrix, method = "ward.D2")
    
    # Cut the dendrogram to get desired number of clusters
    cluster_assignments <- cutree(hclust_result, k = numcluster)
    
    # Create result object with same structure as kmeans
    kmeans_result <- list(cluster = cluster_assignments)
    
    # Remap the cluster to each individual from the evaluation set
    df.eval$cluster <- kmeans_result$cluster[df.eval$class]
  }else{
    kmeans_result <- list(cluster = 1:nrow(embeddings))
    df.eval$cluster <- kmeans_result$cluster[df.eval$class]
  }
  ## construct the estimated cutoff values for the observed outcomes
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
                                           counterfactual = FALSE
                                 )
                               },
                               interval = quantile(dfCluster$R.wS, c(0.01, 0.99),
                                                   na.rm = TRUE
                               ),
                               extendInt = "yes",
                               maxiter = 100
                             )$root,
                             error = function(e) {
                               quantile(dfCluster$R.wS, 1 - alphaCI,
                                        na.rm = TRUE
                               )
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
                       counterfactual = FALSE
             )
           },
           interval = quantile(dfCluster$R.wS, c(0.01, 0.99),
                               na.rm = TRUE
           ),
           extendInt = "yes",
           maxiter = 100
         )$root,
         error = function(e) {
           quantile(dfCluster$R.wS, 1 - alphaCI,
                    na.rm = TRUE
           )
         }
       )
     }) %>% array2DF() %>%
    arrange(cluster) %>% pull(Value)
  ## construct the prediction sets

  if (outcome.type == "Continuous") {
    # conformal inference on the observed outcomes
    # wS
    cls <- kmeans_result$cluster[df.eval$class]
    lower.Y1_A1.wS <- c(predict(m.X1q1.obj,
                                newdata = df.eval[, grep("([XA])",
                                                         colnames(df.eval),
                                                         value = TRUE
                                )]
    )$pred) -
      cutoff.Y1.wS_obs[cls]

    upper.Y1_A1.wS <- c(predict(m.X1q2.obj,
                                newdata = df.eval[, grep("([XA])",
                                                         colnames(df.eval),
                                                         value = TRUE
                                )]
    )$pred) +
      cutoff.Y1.wS_obs[cls]

    lower.Y0_A0.wS <- c(predict(m.X0q1.obj,
                                newdata = df.eval[, grep("([XA])",
                                                         colnames(df.eval),
                                                         value = TRUE
                                )]
    )$pred) -
      cutoff.Y0.wS_obs[cls]

    upper.Y0_A0.wS <- c(predict(m.X0q2.obj,
                                newdata = df.eval[, grep("([XA])",
                                                         colnames(df.eval),
                                                         value = TRUE
                                )]
    )$pred) +
      cutoff.Y0.wS_obs[cls]

    ## for ITE
    # construct the estimated cutoff values for the counterfactual outcomes
    cutoff.Y1.wS <- by(df.eval, list(cluster = df.eval$cluster),
                           function(dfCluster){
                             tryCatch(
                               uniroot(
                                 f = function(theta) {
                                   eif_theta(theta,
                                             df.eval = dfCluster,
                                             model.obj = model.obj.wS,
                                             target.A = 1, wS = TRUE,
                                             alphaCI = alphaCI
                                   )
                                 },
                                 interval = quantile(dfCluster$R.wS, c(0.01, 0.99),
                                                     na.rm = TRUE
                                 ),
                                 extendInt = "yes",
                                 maxiter = 100
                               )$root,
                               error = function(e) {
                                 quantile(dfCluster$R.wS, 1 - alphaCI,
                                          na.rm = TRUE
                                 )
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
                                             alphaCI = alphaCI
                                   )
                                 },
                                 interval = quantile(dfCluster$R.wS, c(0.01, 0.99),
                                                     na.rm = TRUE
                                 ),
                                 extendInt = "yes",
                                 maxiter = 100
                               )$root,
                               error = function(e) {
                                 quantile(dfCluster$R.wS, 1 - alphaCI,
                                          na.rm = TRUE
                                 )
                               }
                             )
                           }) %>% array2DF() %>%
      arrange(cluster) %>% pull(Value)
    # wS
    lower.Y1_A0.wS <- c(predict(m.X1q1.obj,
                                newdata = df.eval[, grep("([XA])",
                                                         colnames(df.eval),
                                                         value = TRUE
                                )]
    )$pred) -
      cutoff.Y1.wS[cls]

    upper.Y1_A0.wS <- c(predict(m.X1q2.obj,
                                newdata = df.eval[, grep("([XA])",
                                                         colnames(df.eval),
                                                         value = TRUE
                                )]
    )$pred) +
      cutoff.Y1.wS[cls]

    lower.Y0_A1.wS <- c(predict(m.X0q1.obj,
                                newdata = df.eval[, grep("([XA])",
                                                         colnames(df.eval),
                                                         value = TRUE
                                )]
    )$pred) -
      cutoff.Y0.wS[cls]

    upper.Y0_A1.wS <- c(predict(m.X0q2.obj,
                                newdata = df.eval[, grep("([XA])",
                                                         colnames(df.eval),
                                                         value = TRUE
                                )]
    )$pred) +
      cutoff.Y0.wS[cls]

    # inference on the observed outcomes
    lower.Y.wS <- lower.Y1_A1.wS * (df.eval$A == 1) +
      lower.Y0_A0.wS * (df.eval$A == 0)
    upper.Y.wS <- upper.Y1_A1.wS * (df.eval$A == 1) +
      upper.Y0_A0.wS * (df.eval$A == 0)

    # inference on tau when A = 0 + A = 1
    lower.tau.wS <- (lower.Y1_A0.wS - (df.eval$Y)) * (df.eval$A == 0) +
      ((df.eval$Y) - upper.Y0_A1.wS) * (df.eval$A == 1)
    upper.tau.wS <- (upper.Y1_A0.wS - (df.eval$Y)) * (df.eval$A == 0) +
      ((df.eval$Y) - lower.Y0_A1.wS) * (df.eval$A == 1)

    # organize the results
    df.eval <- cbind(df.eval,
                     lower.Y = lower.Y.wS,
                     upper.Y = upper.Y.wS,
                     lower.tau = lower.tau.wS,
                     upper.tau = upper.tau.wS
    )
  }
  if (outcome.type == "Categorical") {

    R.XSY1Mat.XS.eval <- R.XSY1Mat.XS[eval.idx, ]
    R.XSY0Mat.XS.eval <- R.XSY0Mat.XS[eval.idx, ]
    # for A = 0, Y1
    Y1predSet.A1wS <- sapply(1:nrow(R.XSY1Mat.XS.eval), function(idx){
      x <- R.XSY1Mat.XS.eval[idx, ]
      cls <- kmeans_result$cluster[df.eval$class[idx]]
      level.Y <- colnames(probY1.XS)[order(x)]
      level.Y[which(x[order(x)] <= cutoff.Y1.wS_obs[cls])]
    }, simplify = FALSE)

    Y0predSet.A0wS <- sapply(1:nrow(R.XSY0Mat.XS.eval), function(idx){
      x <- R.XSY0Mat.XS.eval[idx, ]
      cls <- kmeans_result$cluster[df.eval$class[idx]]
      level.Y <- colnames(probY0.XS)[order(x)]
      level.Y[which(x[order(x)] <= cutoff.Y0.wS_obs[cls])]
    }, simplify = FALSE)

    # Y1predSet.A1wS <- apply(R.XSY1Mat.XS, 1, function(x) {
    #   level.Y <- colnames(probY1.XS)[order(x)]
    #   level.Y[which(x[order(x)] <= cutoff.Y1.wS_obs)]
    # }, simplify = FALSE)
    #
    # # for A = 1, Y0
    # Y0predSet.A0wS <- apply(R.XSY0Mat.XS, 1, function(x) {
    #   level.Y <- colnames(probY0.XS)[order(x)]
    #   level.Y[which(x[order(x)] <= cutoff.Y0.wS_obs)]
    # }, simplify = FALSE)[eval.idx]

    # organize the output for the evaluation sets
    ## for the observed one
    ## with surrogates
    df.eval$sets_observed <- Y0predSet.A0wS
    df.eval$sets_observed[df.eval$A == 1] <-
      Y1predSet.A1wS[df.eval$A == 1]
  }

  # return the evaluation df without non-conformal scores
  return(df.eval %>% dplyr::rename(scores = R.wS))
}
