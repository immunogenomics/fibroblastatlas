colors_overload <- union(ggthemes::tableau_color_pal('Tableau 20')(20), RColorBrewer::brewer.pal(12, 'Set3'))
colors_overload <- c(colors_overload, 'black')

fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
}


do_umap <- function(
    Xmat, cache_fname=NULL, 
    .spread=0.3, .min_dist=0.05,
    .metric='euclidean', .init='laplacian',
    .a=NULL, .b=NULL,
    .n_components=2L,
    ...
) {
    umap_object <- uwot::umap(
        X = Xmat,
        n_threads = 20,
        n_neighbors = 30L,
        n_components = .n_components,
        metric = .metric,
        init = .init,
        n_epochs = NULL,
        learning_rate = 1.0,
#         min_dist = 0.3,
#         spread = 1.0,
        min_dist = .min_dist, 
        spread = .spread,        
        set_op_mix_ratio = 1.0,
        local_connectivity = 1L,
        repulsion_strength = 1,
        negative_sample_rate = 1,
        a = .a,
        b = .b,
        fast_sgd = FALSE,
        verbose = FALSE, 
#         ret_model = TRUE,
#         ret_nn = TRUE
        ret_extra = c('nn', 'fgraph', 'model'),
        ...
    ) 

    ## save object for mapping new data
    if (!is.null(cache_fname)) {
        uwot::save_uwot(umap_object, file = cache_fname)#, unload = FALSE, verbose = FALSE)
    }

    ## fxn from dist to kernel from UWOT 
    nn_idx <- umap_object$nn[[1]]$idx
    adj <- Matrix::sparseMatrix(
        i = rep(1:nrow(nn_idx), each = ncol(nn_idx)), 
        j = c(t(nn_idx)), 
        x = c(t(exp(-(pmax(umap_object$nn[[1]]$dist, .min_dist) - .min_dist)/.spread)))
    )
    
    ## return embeddings
    return(list(
        embedding=umap_object$embedding,
        adj=umap_object$fgraph + Matrix::Diagonal(n = nrow(umap_object$fgraph)),
        knnadj=adj
    ))
}


make_predictions <- function(ref_obj, labels_to_predict, query_obj) {
    if (!labels_to_predict %in% names(ref_obj$models)) {
        stop(as.character(glue('{labels_to_predict} not available to predict in this Atlas.')))
    }
    ## per-cluster predictions
    pred <- map(ref_obj$models[[labels_to_predict]], do.call, list(.xnew=t(query_obj$Z)))

    ## dot product with cluster probabilities
    map(seq_len(nrow(query_obj$R)), function(k) {
        Diagonal(x = t(query_obj$R[k, ])) %*% pred[[k]]
    }) %>% 
        purrr::reduce(`+`)

}     

normalizeData <- function(A, scaling_factor=10000, method='log') {
    if (!"dgCMatrix" %in% class(A)) 
        A <- as(A, "dgCMatrix")
    if (method == 'log') {
        if ('dgCMatrix' %in% class(A)) {
            A@x <- A@x/rep.int(Matrix::colSums(A), diff(A@p))
            A@x <- scaling_factor * A@x
            A@x <- log(1 + A@x)            
        } else {
            
        }
        
    } else {
        stop('Only logCPX transform implemented here')
    }
}

