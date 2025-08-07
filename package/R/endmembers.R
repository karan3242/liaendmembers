
#' Find Endmembers
#'
#' @description
#' Finds endmembers from a set of Lead Isotope Points using Principle Component
#' analysis and the Geochron slope according to the two-stage model.
#'
#' @param pb64,pb74,pb84 Vectors of Pb Isotope rations 206/204, 207/204, 208/204
#' @param tolerance Tolerance value for points considered to be intercepted.
#'      (Default 0.01)
#' @param ... Additional Parameters
#'
#' @returns
#' An object of class 'liaendmembers' as a list of 5
#' @export
#' @examples
#' # Create object with class liaendmembers
#' data("dor_silver_hoard")
#' end_members <- endmembers(
#'         pb64 = dor_silver_hoard[[1]],
#'         pb74 = dor_silver_hoard[[2]],
#'         pb84 = dor_silver_hoard[[3]],
#'         tolerance = 0.01
#' )
#' # Prind summary of the liaendmembers object
#' summary.liaendmembers(end_members)
endmembers <- function(pb64, pb74, pb84, tolerance = 0.01, ...) {
        requireNamespace("stats")
        # Creat a matrix of isotope rations for PCA
     safe_cbind <- function(pb64, pb74, pb84, ...) {

          # 1. Check for equal lengths
          if (length(pb64) != length(pb74) || length(pb64) != length(pb84)) {
               stop("Error: All vectors must have the same length.")
          }

          # 2. Check for NA values
          if (any(is.na(pb64)) || any(is.na(pb74)) || any(is.na(pb84))) {
               stop("Error: Vectors cannot contain NA values.")
          }

          # 3. Check if they are numeric
          if (!is.numeric(pb64) || !is.numeric(pb74) || !is.numeric(pb84)) {
               stop("Error: All vectors must be numeric.")
          }

          # If all checks pass, combine the vectors
          return(cbind(pb64, pb74, pb84))
     }
     isotope_matrix <- safe_cbind(pb64, pb74, pb84)
     row.names(isotope_matrix) <- seq(1, nrow(isotope_matrix))
     # Concduct PCA and check if there are only 2 end memebrs
     pca_result <- prcomp(isotope_matrix, scale = FALSE)
     pca_values <- pca_result$x
     row.names(pca_values) <- row.names(isotope_matrix)
     sum <- summary(pca_result)
     PC1_vla <- sum$importance[["Cumulative Proportion", "PC1"]]
     if (PC1_vla < 0.95){
          warning("PC1 represents less than 95% of the Variance, There may be more than two end members.")
          print(sum)
     }

     # Check if PC2 and PC3 are normally distributed
     norm_test_pc1 <- shapiro.test(pca_values[,"PC1"])$p.value < 0.95
     norm_test_pc2 <- shapiro.test(pca_values[,"PC2"])$p.value > 0.95
     norm_test_pc3 <- shapiro.test(pca_values[,"PC3"])$p.value > 0.95

     if (!(norm_test_pc1 & norm_test_pc2 & norm_test_pc3)){
          warning("PC2 or PC3 are not normally distributed. This may indicate that their variation may not be random noise.")
     }

     # Derive end Memebrs
     pca_ends <- pca_values[pca_values[,"PC1"] %in% c(min(pca_values[,"PC1"]), max(pca_values[,"PC1"])),]

     isotpe_ends <- isotope_matrix[as.numeric(row.names(pca_ends)),]

     geo_slope = 0.626208


     end_member_filter <- function(x, ...){
          geo_intercept <- isotpe_ends[x, "pb74"] - isotpe_ends[x, "pb64"] * geo_slope
          point_itercept <- geo_slope * isotope_matrix[,"pb64"] + geo_intercept
          prob_end <- abs(isotope_matrix[,"pb74"] - point_itercept) < tolerance
          isotope_matrix[prob_end,]
     }

     end_group1 <- end_member_filter(1)
     end_group2 <- end_member_filter(2)
     if(nrow(end_group1) < 2 || nrow(end_group2) < 2){
          warning("End Member group has less than two points. Likely hood of the point being an endmember is low")
     }
     mixing_group <- isotope_matrix[-as.integer(c(row.names(end_group1), row.names(end_group2))),]

     endmember_list <- list(group1 = end_group1,
                            group2 = end_group2,
                            mixing = mixing_group,
                            tolarance = tolerance,
                            pca = pca_result)
     endmember_list <- structure(endmember_list, class = "liaendmembers")
}
