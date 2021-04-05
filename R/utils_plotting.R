## It's a nice theme! 
theme_set(theme_tufte(base_size = 14, base_family = "Helvetica"))

## Tufte is nice but removes boundaries from plots. Add it back here. 
theme_update(panel.background = element_rect(fill = NA, color = "black"))

theme_update(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
)

opts <- options()  # save old options
# options(ggplot2.continuous.colour="viridis")
# options(ggplot2.continuous.fill = "viridis")
# scale_colour_discrete <- function(...) {
#   scale_colour_manual(..., values = viridis::viridis(2))
# }





################# FOREST PLOTS ###########################
forest_multi <- function(data_df, sigma_max = Inf, reorder_levels=TRUE, grid.switch = NULL) {
    if (reorder_levels) {
        cluster_levels <- data.table(data_df)[, .(Z = mean(beta)), by = Cluster][order(Z), as.character(Cluster)]
        data_df <- data_df %>% 
            dplyr::mutate(Cluster = factor(Cluster, cluster_levels)) 
    }
    
    data_df %>% 
        ## shrink sigma if it's not plotable (make a note?)
        dplyr::mutate(sigma = pmin(sigma, sigma_max)) %>% 
        ggplot(aes(Tissue, beta, color = Tissue)) + 
            geom_point() + 
            geom_errorbar(aes(ymin = beta - 1.96 * sigma, ymax = beta + 1.96 * sigma), width = 0) + 
            coord_flip() + 
#             theme_test(base_size = 14) + 
#             scale_color_tableau() + 
            facet_grid(Cluster~., space = 'free', scales = 'free', switch = grid.switch) + 
            geom_hline(yintercept = c(0), linetype = 2) + 
#             geom_hline(yintercept = c(-1, 0, 1), linetype = 2) + 
            labs(x = '', y = 'Log2 fold change') + 
            geom_hline(yintercept = c(0), linetype = 2) + 
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = 'bottom',
                strip.text.y = element_text(angle = 0, hjust = 0),
                panel.background = element_blank()
            ) + 
#             guides(color = FALSE) + 
            NULL    
    
}

forest_uni <- function(data_df, sigma_max = Inf, fdr_max=0.05, beta_min=-Inf) {    
    ## order clusters by effect size and direction
    cluster_levels <- data.table(data_df)[, .(Z = mean(beta)), by = Cluster][order(Z), as.character(Cluster)]

    data_df %>% 
        dplyr::mutate(Cluster = factor(Cluster, cluster_levels)) %>% 
        ## shrink sigma if it's not plotable (make a note?)
        dplyr::mutate(sigma = pmin(sigma, sigma_max)) %>% 
        ggplot(aes(Cluster, beta, color = fdr < fdr_max & beta > beta_min)) +     
            geom_point() + 
            geom_errorbar(aes(ymin = beta - 1.96 * sigma, ymax = beta + 1.96 * sigma), width = 0) + 
            coord_flip() + 
#             theme_test(base_size = 14) + 
            scale_color_manual(values = c('black', 'red')) + 
            geom_hline(yintercept = c(0), linetype = 2) + 
#             geom_hline(yintercept = c(-1, 0, 1), linetype = 2) + 
            labs(x = '', y = 'Log2 fold change') + 
            guides(color = FALSE) + 
            NULL    
    
}


################# FOREST PLOTS ###########################




do_scatter <- function (umap_use, meta_data, label_name, facet_var, no_guides = TRUE, 
    do_labels = TRUE, nice_names, palette_use = colors_overload, 
    pt_size = 4, point_size = 0.5, pt_shape = ".", base_size = 20, 
    do_points = TRUE, do_density = FALSE, h = 3, w = 4, 
                        alpha_fore=1, alpha_back=.3, color_back='lightgrey', 
                       nrow = 1, do_raster = FALSE) 
{
    if (do_raster) {
        geom_point_fxn <- function(...) geom_point_rast(..., width = w, height = h)
    } else {
        geom_point_fxn <- geom_point
    }
    
    plt_df <- data.frame(umap_use)[, 1:2]
    colnames(plt_df) <- c('X1', 'X2')
    plt_df <- plt_df %>% 
        cbind(meta_data) %>% 
        dplyr::sample_frac(1L)
    plt_df$given_name <- plt_df[[label_name]]
    if (!missing(nice_names)) {
        plt_df %<>% dplyr::inner_join(nice_names, by = "given_name") %>% 
            subset(nice_name != "" & !is.na(nice_name))
        plt_df[[label_name]] <- plt_df$nice_name
    }
    
    plt <- plt_df %>% ggplot(aes_string("X1", "X2", col = label_name, 
        fill = label_name)) + 
#         theme_tufte(base_size = base_size) + 
#         theme(panel.background = element_rect(fill = NA, color = "black")) + 
        guides(color = guide_legend(override.aes = list(stroke = 1, 
            alpha = 1, shape = 16, size = 4)), alpha = FALSE) + 
        scale_color_manual(values = palette_use) + scale_fill_manual(values = palette_use) + 
        theme(plot.title = element_text(hjust = 0.5)) + labs(x = "UMAP 1", 
        y = "UMAP 2")
    if (do_points) {
        ## this facets while keeping non-facet points in the background
        if (!missing(facet_var)) {
            if (!is(facet_var, 'quosure')) {
                stop('facet_var must be a quosure. e.g. quo(\'donor\')')
            }            

            plt <- plt + geom_point_fxn(
                data = dplyr::select(plt_df, -!!facet_var), 
                shape = pt_shape, size = point_size,
                color = color_back, fill = color_back, alpha = alpha_back
            ) +
                facet_wrap(vars(!!facet_var), nrow = nrow)
        }
        plt <- plt + geom_point_fxn(shape = pt_shape, size = point_size, alpha = alpha_fore)
    }
    if (do_density) 
        plt <- plt + geom_density_2d()
    if (no_guides) 
        plt <- plt + guides(col = FALSE, fill = FALSE, alpha = FALSE)
    if (do_labels) {
        plt <- plt + 
#             geom_text_repel(
#                 data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], 
#                 label.size = NA, aes_string(label = label_name), 
#                 color = "black", 
#                 size = pt_size, alpha = 1, segment.size = 0
#             ) + 
            geom_label(
                data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], 
                label.size = NA, aes_string(label = label_name, color = label_name), 
#                 color = "black", 
                fill = 'white', 
                size = pt_size, alpha = .6, segment.size = 0
            ) + 
            geom_text(
                data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], 
                label.size = NA, aes_string(label = label_name, color = label_name), 
#                 color = "black", 
                size = pt_size, alpha = 1, segment.size = 0
            ) + 
            guides(col = FALSE, fill = FALSE)
    }
    return(plt)
}


setupVals <- function(data_mat, feature, qlo, qhi) {
    .x <- data_mat[feature, , drop = FALSE] %>% as("dgTMatrix")
    cutoffs <- quantileSparse(.x, c(qlo, qhi))
    cutoffs[2] <- max(cutoffs[2], min(.x@x))
    if (qlo == 0 & qhi == 1) {
        return(.x)
    } 
    
    if (qlo > 0) {
        .x@x[.x@x < cutoffs[1]] <- cutoffs[1]
#         message(sprintf("For %s, lo = %.3f", feature, ifelse(length(.x@x) == ncol(.x), cutoffs[1], NA)))
    }
    if (qhi < 1) {
        .x@x[.x@x > cutoffs[2]] <- cutoffs[2]
#         message(sprintf("For %s, hi = %.3f", feature, cutoffs[2]))
        
    }
    return(.x)
}


quantileSparse <- function(.x, qlist) {
    ratio_zero <- 1 - (length(.x@x) / ncol(.x))
    q_nz <- which(qlist > ratio_zero)
    q_adj <- (qlist[q_nz] - ratio_zero) / (1 - ratio_zero)
    res <- rep(0, length(qlist))
    res[q_nz] <- quantile(.x@x, q_adj)
    res
}

## TODO: test is feature is present
## TODO: allow for different cutoffs, for each marker
## TODO: somehow draw canvas first, then do plotting? 
library(patchwork)
library(ggthemes)

plotFeatures <- function(data_mat, dim_df, features, nrow = 1, 
                         qlo = 0.05, qhi = 1, order_by_expression = FALSE, 
                         pt_shape = 16, pt_size = .5, no_guide = FALSE,
                         .xlim = c(NA, NA), .ylim = c(NA, NA), color_high = muted("blue")) {
    plt_df <- data.frame(dim_df[, 1:2])
    colnames(plt_df) <- c("X1", "X2")


    plt_list <- lapply(features, function(feature) {
        .x <- setupVals(data_mat, feature, qlo, qhi)
        plt_df$value <- 0
        plt_df[.x@j + 1, "value"] <- .x@x
        if (order_by_expression) {
            plt_df %<>% dplyr::arrange(value)             
        } else {
            plt_df %<>% dplyr::sample_frac(1L)
        }

        plt <- plt_df %>% 
            ggplot(aes(X1, X2, color = value)) + 
#             geom_point_rast(dpi = 300, width = 6, height = 4, size = .5, shape = pt_shape) + 
            geom_point(shape = ".") + 
            scale_color_gradient2(na.value = "lightgrey", mid = "lightgrey", midpoint = 0, high = color_high) + 
#             theme_tufte(base_size = 14, base_family = "Helvetica") + 
#             theme(panel.background = element_rect(), plot.title = element_text(hjust = .5)) +
            theme(plot.title = element_text(hjust = .5)) +
            labs(x = "UMAP 1", y = "UMAP 2", title = feature) + 
            NULL
        if (no_guide) {
            plt <- plt + 
            guides(color = FALSE) 
        }
        
        if (sum(is.na(.xlim)) < 2) 
            plt <- plt + xlim(.xlim)
        if (sum(is.na(.ylim)) < 2) 
            plt <- plt + ylim(.ylim)
        plt

    })
    if (length(plt_list) > 1) {
        Reduce(`+`, plt_list) + patchwork::plot_layout(nrow = nrow)
    } else {
        plt_list[[1]]
    }
}

