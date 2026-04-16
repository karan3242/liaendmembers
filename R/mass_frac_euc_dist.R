#' Euclidean Distance for Pb Isotope rations to ore Sources
#'
#' This function
#'
#' @param x matrix of Pb Isotopes with colums in the order of
#' 206Pb/204Pb, 207Pb/204Pb,208Pb/204Pb
#' @param ox dataframe of renfrence data contained the Pb Isotope rations
#' @param ox_col vector contained names of the Pb Isotope rations.
#' Names much include digits 6, 7, 8 for isotope ratio identification.
#' @param group groups belonging to reach isotope ration vector
#' @param .n Length of result output (Default = 10)
#' @param as.dist Option to show distances or only a vector of closest groups.
#' @param ...
#'
#' @returns data frame or charecter vector
#' @export
#'
euc_dist <- function(x,
                     ox,
                     ox_col,
                     group,
                     .n = 10,
                     as.dist = FALSE ,
                     ...) {
        if (!all(ox_col %in% names(ox))) {
                stop("column names not found")
        }
        if (nrow(ox) != length(group)) {
                stop("ox rows and groups not equal")
        }
        x1 <- as.matrix(x)
        ox1 <- ox[, ox_col]
        ox1 <- ox1[, c(grep("6", names(ox1), ),
                       grep("7", names(ox1)),
                       grep("8", names(ox1)))]
        ox1 <- as.matrix(ox[, ox_col])
        ox_group <- group

        norm_x1 <- rowSums(x1^2)
        norm_ox1 <- rowSums(ox1^2)

        dot_product <- x1 %*% t(ox1)

        dist_sq <- sweep(sweep(-2 * dot_product, 1, norm_x1, "+"), 2, norm_ox1, "+")

        dist_matrix <- t(sqrt(pmax(dist_sq, 0)))
        min_dists <- apply(dist_matrix, 1, min)
        hits <- data.frame(group = ox_group, dist = min_dists)
        hits_sorted <- hits[order(hits$dist), ]

        if (as.dist) {
                return(head(hits_sorted, .n))
        }
        sorted_unique_groups <- unique(hits_sorted$group)
        head(sorted_unique_groups, n = .n)
}
#' @rdname euc_dist
#' @param s mass-fractionation factor (Defualt = 0.001)
#' @export
mf_dist <- function(x,
                    ox,
                    ox_col,
                    group,
                    .n = 10,
                    as.dist = FALSE,
                    s = 0.001,
                    ...) {
        if (!all(ox_col %in% names(ox))) {
                stop("column names not found")
        }
        if (nrow(ox) != length(group)) {
                stop("ox rows and groups not equal")
        }
        x <- as.matrix(x)
        ox <- ox[, ox_col]
        ox <- ox[, c(grep("6", names(ox), ),
                     grep("7", names(ox)),
                     grep("8", names(ox)))]
        ox <- as.matrix(ox)

        if (nrow(ox) != length(group)) {
                stop("Groups vector not the same length as ox")
        }

        n_artifacts <- nrow(x)
        n_ores <- nrow(ox)

        # Constant Correlation Matrix R
        R <- matrix(
                c(1, 0.96, 0.94,
                  0.96, 1, 0.96,
                  0.94, 0.96, 1),
                nrow = 3,
                byrow = TRUE
        )

        # Long-form is usually better for ggplot2 or further analysis
        results_list <- lapply(1:n_artifacts, function(j) {
                x0 <- x[j, ]

                # 2. Geometry Setup for this specific artifact
                v <- x0 * c(2, 3, 4)
                n <- v / sqrt(sum(v^2))

                # Weighting Matrix W calculation (Eq 8)
                sd_diag <- diag(c(2, 3, 4) * s * x0)
                W <- sd_diag %*% R %*% sd_diag

                # Projection setup (Gram-Schmidt)
                basis1 <- if (abs(n[1]) < 0.9)
                        c(1, 0, 0)
                else
                        c(0, 1, 0)
                u1 <- basis1 - (sum(basis1 * n)) * n
                u1 <- u1 / sqrt(sum(u1^2))
                u2 <- c(n[2] * u1[3] - n[3] * u1[2], n[3] * u1[1] - n[1] *
                                u1[3], n[1] * u1[2] - n[2] * u1[1])
                P <- cbind(u1, u2)

                # 3. Project W into 2D and invert
                W_p_inv <- solve(t(P) %*% W %*% P)

                # 4. Vectorized Distance Calculation for all ores against this artifact
                delta_X <- sweep(ox, 2, x0, "-")

                # Project all delta_X rows into 2D: (N x 3) %*% (3 x 2) = (N x 2)
                dx_p <- delta_X %*% P

                # Calculate quadratic form: d_sq = rowSums((dx_p %*% W_p_inv) * dx_p)
                # This is the vectorized version of t(dx_p) %*% inv %*% dx_p
                d_sq <- rowSums((dx_p %*% W_p_inv) * dx_p)

                return(data.frame(region = group, d_sq = d_sq))
        })

        # Combine all artifact results into one data frame
        all_list <- do.call(rbind, results_list)
        all_list <- all_list[order(all_list$d_sq), ]
        if (as.dist) {
                all_list <- head(all_list, .n)
                return(all_list)
        }
        head(unique(all_list$region), .n)

}


#' @rdname euc_dist
#' @param dist_type Distance type to use, simple euclidiant or
#' mass-fractionation corrected
#' @returns list of dataframe or charecter vector
#' @export
endmember_dist <- function(df,
                           ox,
                           ox_col,
                           group,
                           dist_type = "ed",
                           .n = 10,
                           as.dist = FALSE,
                           s = 0.001,
                           ...) {
        if (!"liaendmembers" %in% class(df)) {
                stop("df is not of class 'liaendmembers'")
        }
        if (dist_type != "ed" || dist_type != "mfd"){
                stop("dist_type not defined correctly. Should be either 'ed' or 'mfd")
        }
        if (dist_type == "ed")
        list(
                "group1" = clossest_region(df$group1, ox, ox_col, group, .n, as.dist, ...),
                "group2" = clossest_region(df$group2, ox, ox_col, group, .n, as.dist, ...)
        )
        if (dist_type == "mfd")
        list(
                "group1" = mass_frac_dist_multi(df$group1, ox, ox_col, group, s, .n, as.dist),
                "group2" = mass_frac_dist_multi(df$group2, ox, ox_col, group, s, .n, as.dist)
        )
}

