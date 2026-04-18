# Refrence Data Function --------------------------------------------------

#' Cereate Reference data object for LIA endmember distance and probablity estimate functions.
#'
#' @param x Reference data frame or tibble.
#' @param cols Nector of column names which represent 206Pb/204Pb, 207Pb/204Pb, 208Pb/204Pb
#' @param group Name of the column containing siotope groups as character string.
#'
#' @returns
#' ref.data object with Isotope groupings and Isotope colums of
#' 206Pb/204Pb, 207Pb/204Pb, 208Pb/204Pb
#' with shorter names.
#' @export
as.ref_data <- function(x, cols, group) {
        if (!all(cols %in% names(x))) {
                stop("column names not found")
        }
        cols <- c(grep("6", cols, value = T),
                  grep("7", cols, value = T),
                  grep("8", cols, value = T))

        x <- na.omit(x[, c(group, cols)])
        names(x) <- c("groups", "pb64", "pb74", "pb84")
        class(x) <- c("ref.data", class(x))
        x
}

# Training Function ----------------------------------------------------------------

#' Train XGBOOST Model
#'
#' Trains XGBOOST model for predictive isotope analysis using DCSCAN clustring
#' SMOTE data imputation for training data set, and XGBOOST, following the method
#' of Shnyr et.al (2026)
#'
#' @references Shnyr, E., Kuflik, T., Desai, K., & Eshel, T. (2026). Determining the origins of Phoenician silver: Exploring the potential of machine learning for lead isotope analysis. Journal of Archaeological Science, 188, 106–499. https://doi.org/10.1016/j.jas.2026.106499
#'
#' @param ref `ref.data` object created by ref_data.
#' @param .minSize Minisume number of samples in group to be used for sampeling (Default = 3).
#' @param .minPts_fac scaling factor minimum needed points from each group, ranging from 0:1. (Default = 0.1)
#' @param .eps size (radius) of the epsilon neighborhood. (Default = 0.18)
#' @param .eta Step size shrinkage used in update to prevent overfitting. After each boosting step, we can directly get the weights of new features, and eta shrinks the feature weights to make the boosting process more conservative. (Defualt = 0.1)
#' @param .max_depth Maximum depth of a tree. Increasing this value will make the model more complex and more likely to overfit. 0 indicates no limit on depth. Beware that XGBoost aggressively consumes memory when training a deep tree. "exact" tree method requires non-zero value. (Default = 6)
#' @param .nrounds Max number of boosting iterations. (Default = 100)
#' @param nthread Number of threadts for parallel processing. When choosing it, please keep thread contention and hyperthreading in mind. (Default = 4)
#' @inheritDotParams dbscan::dbscan weights borderPoints
#' @inheritDotParams xgboost::xgb.train early_stopping_rounds maximize
#'
#' @importFrom dbscan dbscan
#' @importFrom smotefamily SMOTE
#' @importFrom xgboost xgb.train xgb.DMatrix
#' @returns
#' List of xgboot.model objects.
#' @export
train_data <- function(ref,
                       .minSize = 3,
                       .minPts_fac = 0.1,
                       .eps = 0.18,
                       .eta = 0.1,
                       .max_depth = 6,
                       .nrounds = 100,
                       nthread = 4L,
                       ...) {
        if (!"ref.data" %in% class(ref)) {
                stop(
                        "ref must be of class ref.data. Use function `ref_data`
                     for preporcessing ref data first."
                )
        }
        ox <- ref
        uni_groups <- unique(ox[[1]])

        # DBSCAN ------------------------------------------------------------------
        dbscan_groups <- function(g_name) {
                group_df <- ox[ox[[1]] == g_name, ]

                # Ensure we have enough data to cluster
                if (nrow(group_df) < .minSize) {
                        return(NULL)
                }
                if(!.minPts_fac > 0 & !.minPts_fac < 1 ) {
                        stop(".minPts_fac should be bettwee 0 or 1")
                }
                minPts <- nrow(group_df) * .minPts_fac
                res <- dbscan::dbscan(group_df[, -1], minPts = minPts , eps = .eps, ...)
                cluster <- res$cluster
                group_df <- group_df[cluster > 0, ]
                cluster_s <- cluster[cluster > 0]
                if(sum(cluster_s) <= 0){
                        return(NULL)
                }
                if(length(unique(cluster_s))> 1){
                        group_df[[1]] <- paste0(group_df[[1]], "_", cluster_s)
                }
                return(group_df)

        }
        dbscan_df <- do.call(rbind, lapply(uni_groups, dbscan_groups))
        # SMOTE (Multi-class handling) --------------------------------------------
        subgroups <- unique(dbscan_df[[1]])
        apply_smote <- function(target_subgroup) {
                temp_df <- dbscan_df
                temp_df$target <- ifelse(temp_df[[1]] == target_subgroup, "p", "n")

                p_count <- sum(temp_df$target == "p")

                # SMOTE requires at least K+1 samples in the minority class
                if (p_count < 3)
                        return(NULL)

                k_val <- min(2, p_count - 1)

                smote_out <- smotefamily::SMOTE(temp_df[, 2:4],
                                                target = temp_df$target,
                                                K = k_val)

                final_df <- rbind(smote_out$orig_P, smote_out$syn_data)
                final_df$group <- target_subgroup
                final_df
        }
        smote_df_final <- do.call(rbind, lapply(subgroups, apply_smote))
        # XGBOOST Implementation --------------------------------------------------
        # Define a function to train a binary model for a specific group
        train_group_model <- function(target_group, data = smote_df_final) {
                # Create binary labels: 1 if target_group, 0 otherwise
                labels <- ifelse(data$group == target_group, 1, 0)

                # Check if we have both classes represented
                if (length(unique(labels)) < 2) {
                        warning(
                                paste(
                                        "Skipping group",
                                        target_group,
                                        "- no negative samples available."
                                )
                        )
                        return(NULL)
                }

                X_train <- as.matrix(data[, 1:3])
                dtrain <- xgboost::xgb.DMatrix(data = X_train, label = labels)

                # Binary logistic parameters
                params <- list(
                        objective = "binary:logistic",
                        eta = .eta,
                        max_depth = .max_depth,
                        nthread = nthread,
                        eval_metric = "logloss",
                        # scale_pos_weight can help if your group is much smaller than the rest
                        scale_pos_weight = sum(labels == 0) / sum(labels == 1)
                )

                model <- xgboost::xgb.train(params = params,
                                            data = dtrain,
                                            nrounds = .nrounds,
                                            ...)

                return(model)
        }
        list <- setNames(lapply(subgroups, train_group_model), subgroups)
}


# Prediction Function -----------------------------------------------------

xgboost_predict <- function(x, model_list = NULL, .n) {
        dtest <- xgboost::xgb.DMatrix(as.matrix(x))

        all_preds <- lapply(model_list, function(m) {
                if (is.null(m))
                        return(rep(NA, nrow(x)))
                predict(m, dtest)
        })
        df <- stack(as.data.frame(all_preds))
        names(df) <- c("prob", "group")
        df <- df[, c("group", "prob")]
        df$group <- gsub("\\.", " ", df$group)
        head(df[order(df$prob, decreasing = TRUE), ], .n)
}



# End member Prediction function ------------------------------------------

#' Predict Isotope Provenance
#'
#' Predicuts Pb Isotope provencnace of a sample matrix in reference ot a xgboost trained list
#'
#' @param x Matrix of liaendmembers object of pbisotope samples
#' @param model_list Model list gnerated by `train_data()`
#' @param .n Intiger for numer of observations.
#'
#' @importFrom stats na.omit predict setNames
#'
#' @returns
#' data.frame object or list of data.frames
#' @export
isoprov_predict <- function(x, model_list = NULL, .n = 6) {
     if (inherits(x, "liaendmembers")) {
          target_groups <- x[3:4]
          pred <- lapply(target_groups, \(grp) {
                  xgboost_predict(grp, model_list = model_list, .n = .n)
          })
          return(pred)
     }

     return(xgboost_predict(x, model_list = model_list, .n = .n))
}
