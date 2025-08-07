

#' Find Endmembers
#'
#' @description
#' Finds endmembers from a set of Lead Isotope Points using Principle Component
#' analysis and the Geochron slope according to the two-stage model.
#' @param x data.frame or matrix object or Pb 206/204, 207/204, 208/204 isotope ratios.
#' @param col Isotope column names.
#' @param tolerance Tolerance value for points considered to be intercepted.
#'      (Default 0.01)
#' @param ... Additional Parameters
#'
#' @returns
#' An object of class 'liaendmembers' as a list of 6
#'      \item{data}{Isotope data}
#'      \item{group1, group2}{End member groups}
#'      \item{mixing}{Mixing Group}
#'      \item{tolerance}{Totlarance value}
#'      \item{pca}{List of PCA analysis}
#' @export
#' @examples
#' # Create object with class liaendmembers
#' data("dor")
#' end_members <- endmembers(
#'         dor,
#'         colnames(dor),
#'         tolerance = 0.01
#' )
#' # Prind summary of the liaendmembers object
#' summary.liaendmembers(end_members)
endmembers <- function(x, col = NULL, tolerance = 0.01, ...) {
        requireNamespace("stats")
        if (!inherits(x, "data.frame") && !inherits(x, "matrix")) {
                stop(paste(
                        deparse(substitute(x)),
                        "is not a dataframe or a matrix"
                ))
        }

        if (is.null(col)) {
                stop("column names needed!")
        }

        x <- x[, col]

        if (inherits(x, "data.frame")) {
                x <- as.matrix(x)
        }
        if (!is.numeric(x)) {
                stop("Non-numeric values in dataframe or matrix")
        }
        row.names(x) <- 1:nrow(x)

        isotope_matrix <- x
        if (nrow(isotope_matrix) < 3) {
                stop("To few samples. Need to  be more than 3")
        }
        # Concduct PCA and check if there are only 2 end memebrs
        pca_result <- prcomp(isotope_matrix, scale = FALSE)
        pca_values <- pca_result$x
        row.names(pca_values) <- row.names(isotope_matrix)
        sum <- summary(pca_result)
        PC1_vla <- sum$importance[["Cumulative Proportion", "PC1"]]
        if (PC1_vla < 0.95) {
                warning(
                        "PC1 represents less than 95% of the Variance, There may be more than two end members."
                )
                print(sum)
        }

        # Check if PC2 and PC3 are normally distributed
        norm_test_pc1 <- shapiro.test(pca_values[, "PC1"])$p.value < 0.95
        norm_test_pc2 <- shapiro.test(pca_values[, "PC2"])$p.value > 0.95
        norm_test_pc3 <- shapiro.test(pca_values[, "PC3"])$p.value > 0.95

        if (!(norm_test_pc1 & norm_test_pc2 & norm_test_pc3)) {
                warning(
                        "PC2 or PC3 are not normally distributed. This may indicate that their variation may not be random noise."
                )
        }

        # Derive end Memebrs
        pca_ends <- pca_values[pca_values[, "PC1"] %in% c(min(pca_values[, "PC1"]), max(pca_values[, "PC1"])), ]

        isotpe_ends <- isotope_matrix[as.numeric(row.names(pca_ends)), ]

        geo_slope = 0.626208


        end_member_filter <- function(x, ...) {
                geo_intercept <- isotpe_ends[x, "pb74"] - isotpe_ends[x, "pb64"] * geo_slope
                point_itercept <- geo_slope * isotope_matrix[, "pb64"] + geo_intercept
                prob_end <- abs(isotope_matrix[, "pb74"] - point_itercept) < tolerance
                isotope_matrix[prob_end, ]
        }

        end_group1 <- end_member_filter(1)
        end_group2 <- end_member_filter(2)
        if (nrow(end_group1) < 2 || nrow(end_group2) < 2) {
                warning(
                        "End Member group has less than two points. Likely hood of the point being an endmember is low"
                )
        }
        mixing_group <- isotope_matrix[-as.integer(c(row.names(end_group1), row.names(end_group2))), ]

        endmember_list <- list(
                data = isotope_matrix,
                group1 = end_group1,
                group2 = end_group2,
                mixing = mixing_group,
                tolarance = tolerance,
                pca = pca_result
        )
        endmember_list <- structure(endmember_list, class = "liaendmembers")
        return(endmember_list)
}
