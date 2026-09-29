#' Find Silhouette-Optimized Clusters
#'
#' This helper function assists with finding the optimal number of clusters using the silhouette method.
#'
#' @param hc An object produced by the hclust() function. Hierarchical clustered data.
#' @param dist_obj An object produced by the as.dist() function. A one-dimensional vector of values, calculated as distance between pairwise correlations
#' @param input_counts Counts to determine if silhouette k+1 results in the freeing of a subcluster
#' @param k_min Minimum number of clusters to test
#' @param k_max Maximum number of clusters to test
#'
#' @return A list of labels assigned, the optimal k, and the silhouette curve.
#'
#' @export

best_k_labels <- function(

    hc,
    dist_obj,
    input_counts = NULL,
    k_min = 2L,
    k_max = 10L

    ) {

  kk <- k_min:k_max #k values to test

  curve <- sapply(kk, function(k) { #for each k value...
    lab <- cutree(hc, k = k) #cut the tree there
    if (length(unique(lab)) < 2) return(NA_real_) #NA redundency for 1 cluster
    mean(silhouette(lab, dist = dist_obj)[, "sil_width"]) #take average silhouette score
  })
  if (all(is.na(curve))) #return error if the entire silhouette failed
    stop("No valid k in [", k_min, ",", k_max,
         "] produced >= 2 clusters.")

  best_k <- kk[which.max(curve)] #take the best k
  assign_labels <- cutree(hc, k = best_k)

  curve_d <- curve[1:length(curve)-1] - lead(curve)[1:length(curve)-1]
  best_gain <- kk[which.max(curve_d)]

  split_res <- FALSE
  if(!is.null(input_counts)){
    alt_labels <- cutree(hc, k = (best_k + 1))

    count_labels <- tabulate(assign_labels) %>% setNames(nm = unique(assign_labels))
    count_alt <- tabulate(alt_labels) %>% setNames(nm = unique(alt_labels))

    if(max(count_labels) != max(alt_labels)){


      filter_cluster <- count_alt[!count_alt %in% count_labels]

      # Filter top and find eigenvector
      initial_top <- assign_labels[assign_labels == names(count_labels[which.max(count_labels)])]
      ev_initial <- input_counts[row.names(input_counts) %in% names(initial_top),] %>%
        t() %>% cor(method = "pearson") %>%
        eigs_sym(k = min(c(10, length(initial_top)-1)), which = "LA")
      ev_initial <- ev_initial[["values"]]
      ev_initial <- ev_initial / sum(ev_initial)

      # Filter the larger split and find eigenvector
      split_top <- assign_labels[alt_labels == names(filter_cluster[which.max(filter_cluster)])]
      ev_split_top <- input_counts[row.names(input_counts) %in% names(split_top),] %>%
        t() %>% cor(method = "pearson") %>%
        eigs_sym(k = min(c(10, length(split_top)-1)), which = "LA")
      ev_split_top <- ev_split_top[["values"]]
      ev_split_top <- ev_split_top / sum(ev_split_top)

      # Filter the smaller split and find eigenvector
      split_bot <- assign_labels[alt_labels == names(filter_cluster[which.min(filter_cluster)])]
      ev_split_bot <- input_counts[row.names(input_counts) %in% names(split_bot),] %>%
        t() %>% cor(method = "pearson") %>%
        eigs_sym(k = min(c(10, length(split_bot)-1)), which = "LA")
      ev_split_bot <- ev_split_bot[["values"]]
      ev_split_bot <- ev_split_bot / sum(ev_split_bot)

      # Calculate if the split increases or decreases A2/A1
      if(
        (ev_initial[2]/ev_initial[1]) > ((ev_split_top[2]/ev_split_top[1]) + (ev_split_bot[2]/ev_split_bot[1]))
        ){
        assign_labels <- alt_labels
        best_k <- (best_k + 1)
        split_res <- TRUE
      }
    }
  }

  if(1 %in% tabulate(assign_labels)){3
    best_k <- kk[which.max(curve)+1]
    assign_labels <- cutree(hc, k = best_k)
    bump_k <- TRUE
  }else{bump_k <- FALSE}

  return(list( #return our data
    labels   = assign_labels,
    best_k   = best_k,
    gain_k = best_gain,
    curve    = data.frame(
      k       = kk,
      mean_sil = curve
    ),
    singlet = bump_k,
    split_k = split_res
  ))
}


#' Calculate Cluster Variance and Standard Deviation
#'
#' This helper function takes input data and clusters to calculate each cluster's standard deviation and variance per variable (i.e. time point)
#'
#' @param input_counts A data.frame or matrix object, row names as observations (i.e. genes) and column names as variables (i.e. time points)
#' @param input_clusters A data.frame object with two columns: observations (i.e. genes) and clusters. To allow for different naming conventions, the ORDER is important: initial_clusters column 1 == character class (observations), initial_clusters column 2 == numeric class (cluster assignments).
#'
#' @return A data.frame object containing metrics for the input clusters.
#'
#' @export

calculate_cluster_var <- function(

    input_counts,
    input_clusters

    ) {

  # For each cluster, isolate and find the values
  cluster_total <- unique(input_clusters)

  cluster_res <- map_dfr(cluster_total, \(cluster_i) {
    pickCluster <- input_clusters[input_clusters == cluster_i]

    if(length(pickCluster) > 1){
    pickGenes <- input_counts[names(pickCluster),]
    clusterVar <- apply(pickGenes, 2, var)
    clusterSD <- apply(pickGenes, 2, sd)

    ev_res <- eigen(cor(t(pickGenes),method = "pearson"), symmetric = TRUE, only.values = TRUE)$values
    prop <- ev_res / sum(ev_res)

    clusterInfo <- list( # Combine for graphing
      Var  = mean(clusterVar),
      sum_var = sum(clusterVar),
      SD = mean(clusterSD),
      Cluster = as.character(cluster_i),
      n = length(pickCluster),
      eigen = prop[1:5],
      row.names = NULL
    )
    return(clusterInfo)
    }else{

    clusterInfo <- data.frame( # Combine for graphing
      Var  = 0,
      mean_var = mean(clusterVar),
      SD = 0,
      Cluster = as.character(cluster_i),
      n = length(pickCluster),
      row.names = NULL
    )
    return(clusterInfo)
    }
  })
  return(cluster_res)

}

#' Cluster and Filter Data Using Silhouette Method
#'
#' @param expr A data.frame or matrix object, row names as observations (i.e. genes) and column names as variables (i.e. time points)
#' @param sil_threshold Threshold score, removing genes which score below this value
#' @param cor_method Which method to use for correlation. Pearson is strongly recommended, but option for Kendall or Spearman is present.
#' @param clus_method Which method to use for hierarchical clustering. Complete method allows for non-Euclidian distances, since we are taking pairwise correlations between objects.
#' @param condition_name The name of your condition, which will be paired with the output.
#' @param return_output Input one of "expr", "data", or "metric", or both by default. Will return filtered expression, calculated data, and silhouette metrics, respectively.
#' @param ... Additional parameters to be passed on to [totoInferno::best_k_labels]
#' @inheritDotParams best_k_labels k_min k_max
#'
#' @return A list object containing various data metrics, filter metrics, and returned data.
#'
#' @export

silhouette_filter <- function(

    expr,
    sil_threshold  = 0,
    cor_method     = "pearson",
    clus_method    = "complete",
    condition_name = "condition",
    return_output  = c("expr", "data", "metric"),

    ...

    ) {
  if(!"data" %in% return_output & !"metric" %in% return_output & !"expr" %in% return_output){
    stop("Unrecognized return value, please see documentation for options")
  }

  expr <- as.matrix(expr) #turn into matrix for computational efficiency
  if (is.null(rownames(expr))) stop("Input data must have gene IDs as rownames.")

  ## -- Step 1: distance and clustering ----
  d1 <- as.dist(1 - cor(t(expr), method = cor_method))
  hc1 <- hclust(d1, method = clus_method)

  ## -- Step 2: First cluster assignment (intended for silhouette curve) ----

  # The initial silhouette grouping genes inappropriately is a major source of error
  # Too little or (especially) too many groups affects how this pipeline determines "well fit" clusters

  # If the output is chopping a lot of genes at a low threshold, such as -0.5, the input data is being clustered poorly with the 'complete' method. Use a mean correlation filter to remove some noise and assist the complete clustering. Or (not recommended) use another clustering method.

  cut1 <- best_k_labels(hc1, d1, expr, ...) #take best K
  labels_p1 <- setNames(cut1[["labels"]], row.names(expr))

  n_clusters_p1 <- length(unique(labels_p1))

  ## -- Step 3: Trim silhouette assignment ----

  sil_tab <- list(
    gene       = row.names(expr),
    cluster_p1 = labels_p1
  )

  if (n_clusters_p1 >= 2) {
    sil_tab[["sil_width"]] <- silhouette(labels_p1, dist = d1)[, "sil_width"]
  } else {
    warning("[", condition_name, "] Pass ", sil_threshold, " produced < 2 clusters; ",
            "silhouette undefined.")
  }

  ## Keep genes with silhouette >= threshold.
  keep_mask     <- !is.na(sil_tab[["sil_width"]]) &
    (sil_tab[["sil_width"]] > sil_threshold) #remove genes which have a silhouette score below the threshold

  expr_keep     <- expr[keep_mask, , drop = FALSE]
  n_removed     <- abs(nrow(expr) - nrow(expr_keep))
  rm(d1)
  rm(hc1)
  rm(labels_p1)
  rm(sil_tab)

  if (nrow(expr_keep) < 2)
    stop("[", condition_name, "] Fewer than 2 genes survived the silhouette ",
         "trim (threshold = ", sil_threshold, ")")

  ## -- Step 4: Pass-2 re-clustering on remaining genes ----

  d2 <- as.dist(1 - cor(t(expr_keep), method = cor_method))
  hc2 <- hclust(d2, method = clus_method)

  cut2 <- best_k_labels(hc2, d2, expr_keep, ...)
  labels_final <- setNames(cut2$labels, row.names(expr_keep))

  ## -- Step 6: validation metrics on FINAL clusters ----
  ## mean silhouette of the final clustering
  if (length(unique(labels_final)) >= 2) {
    sil_final <- silhouette(labels_final, dist = d2)
    mean_sil_final <- mean(sil_final[, "sil_width"])
    sil_per_cluster <- tapply(sil_final[, "sil_width"], labels_final, mean)
  } else {
    mean_sil_final  <- NA_real_
    sil_per_cluster <- NA_real_
  } #metrics for reporting

  ## within-cluster SD/SE of the expression profile (per cluster + overall)
  # mean over the 6 time points of each cluster's per-time-point SD.

  find_var_sd     <- calculate_cluster_var(expr_keep, cut2$labels) #overall score for reporting


  ## -- Step 7: Create list of metrics based on how this pipeline operated ----
  # No one metric is the definitive 'grade'
  # multiple metrics give insight into data refinement and changes over threshold iterations
  res <- list()
  if("metric" %in% return_output) {
    res[["metric"]] <- list(
      mean_silhouette_final = cut2$curve[cut2$curve$k == cut2$best_k,]$mean_sil,
      mean_silhouette_second = cut2$curve[cut2$curve$k == cut2$gain_k,]$mean_sil,
      within_cluster_sd_overall = find_var_sd$SD,
      within_cluster_var_overall = find_var_sd$Var,
      genes_removed = n_removed,
      n_clusters_final = cut2$best_k,
      n_clusters_second = cut2$gain_k,
      singlet = cut2$singlet,
      n_genes_in  = nrow(expr),
      n_genes_out = nrow(expr_keep),
      split_k = cut2[["split_k"]]
    )
  }

  ## For scan plotting: silhouette, SD, and 'percent genes remaining' scale naturally to 0-1 scale

  ## -- Step 8: Output list of values plus metrics (Step 7) ----
  if("expr" %in% return_output){
    res[["expr"]] <- expr_keep
  }

  if("data" %in% return_output){
    res[["data"]] <- list(
      genes_keep     = row.names(expr_keep),
      dist_final     = d2,
      hclust_final   = hc2,
      labels_final   = labels_final,    # Named survivors
      sil_curve_p1   = cut1[["curve"]],      # Silhouette curve via silhouette optimized (curve) method
      sil_curve_final= mean_sil_final
    )
  }

  if(length(return_output) == 1){
    res <- list_c(res)
  }
  return(res)

}



#' Apply Silhouette Filter Iterations and Return Data.
#'
#' This function is similar to silhouette_filter_scan but removes many redundancies for reporting how the pipeline works through the data. This increases speed and reduces memory burden. The intended use is with parameters determined through silhouette_filter_scan, then using this function to more efficiently apply the filter and return the data.
#'
#' @param expr A data.frame or matrix object, row names as observations (i.e. genes) and column names as variables (i.e. time points)
#' @param threshold_sequence A vector or list, preferably made by seq(), which contains the thresholds to iterate through
#' @param return_output Input one of "expr", "metric", or both by default. Returns the filtered data or the cluster metrics, respectively.
#' @param verbose Prints the threshold at each step. Useful for monitoring progress.
#' @param ... Addition arguments passed to [totoInferno::silhouette_filter]
#' @inheritDotParams silhouette_filter -expr -sil_threshold -return_output
#' @inheritDotParams best_k_labels k_min k_max
#'
#' @return A list object containing various data metrics, filter metrics, and returned data.
#'
#' @export

silhouette_filter_iterate <- function(

    expr,
    threshold_sequence = seq(-0.5, 0.5, by = 0.01),
    return_output      = c("expr", "metric"),
    verbose            = FALSE,

    ...

    ) {

  expr <- as.matrix(expr)

  ## Run the iterative scan
  metric_res <- data.frame()

  for(i in seq_along(threshold_sequence)){

    n_removed <- 1

    while(n_removed > 0){ #repeat filters until the filter is not removing genes. Gives data a chance to stabilize, especially with larger jumps
      if(verbose == TRUE){print(threshold_sequence[i])}
      ## run the full pipeline at this threshold; capture failures as NA rows
      res_it <- tryCatch(

        silhouette_filter(
          expr = expr,
          sil_threshold = threshold_sequence[i],
          return_output = c("expr", "metric"),
          ...
        ),

        error = function(e) NULL

      )


      if(!is.null(res_it)){
        n_removed <- res_it[["metric"]][["genes_removed"]]
        expr <- res_it[["expr"]]

        if("metric" %in% return_output){
          metric_res <- rbind(metric_res, data.frame(
            threshold = threshold_sequence[i],
            n_remaining = res_it[["metric"]][["n_genes_out"]],
            n_removed = res_it[["metric"]][["genes_removed"]],
            silhouette = res_it[["metric"]][["mean_silhouette_final"]],
            var = mean(res_it[["metric"]][["within_cluster_var_overall"]]),
            sd = mean(res_it[["metric"]][["within_cluster_sd_overall"]])

          ))
        }

      }else{n_removed <- 1}

    }
  }

  res <- list()
  if("metric" %in% return_output){
    res[["metric"]] <- metric_res
  }
  if("expr" %in% return_output){
    res[["expr"]] <- res_it[["expr"]]
  }
  if(length(return_output) == 1){
    res <- list_c(res)
  }

  return(res)
}

#' Apply High Pass Filter and Return Filter Information
#'
#' This function uses the average correlation of observations (i.e. genes) to remove those which have the smallest correlations to the system. This increases the cluster separation when using the 'complete' hierarchical clustering method. Increased cluster separation allows for the silhouette_filter function to handle input data more predictably. This function specifically returns various metrics for use in the cor_filter_scan function to determine the optimal initial high pass filter.
#'
#' @param expr A data.frame object, row names as observations (i.e. genes) and column names as variables (i.e. time points)
#' @param cor_threshold The cutoff to remove genes with the lowest average correlation
#' @param cor_method Which method to use for correlation. Pearson is strongly recommended, but option for Kendall or Spearman is present.
#' @param clus_method Which method to use for hierarchical clustering. Complete method allows for non-Euclidian distances, since we are taking pairwise correlations between objects.
#' @param return_output Input one of "expr", "metric", or both by default. Returns the filtered data or the cluster metrics, respectively.
#' @param return_second_cluster TRUE/FALSE. Determines to pass on a second silhouette-based cluster amount, determined by the largest gain in silhouette width as cluster amount decreases. Useful for separating similar, but distinct groups.
#' @param verbose Print progress messages
#' @param ... Additional parameters to be passed on to [totoInferno::best_k_labels]
#' @inheritDotParams best_k_labels k_min k_max
#'
#' @return A data.frame object containing Standard deviation, Variance, and Silhouette curve.
#'
#' @export

cor_filter <- function(

  expr,
  cor_threshold      = 0,
  cor_method  = "pearson",
  clus_method = "complete",
  return_output = c("expr", "metric"),
  return_second_cluster = FALSE,
  verbose     = FALSE,

  ...

) {
  ## Initial processing of raw df
  mean_cor <- seq(1, nrow(expr)) %>%
    map(\(i) {
      res <- cor(expr[i,], y = t(expr[-i,]), method = cor_method) %>%
        abs() %>%
        mean()
      return(res)
    }) %>%
    list_c() %>%
    setNames(nm = row.names(expr))

  filter_mask <- mean_cor[mean_cor > cor_threshold] #remove genes below the cutoff
  expr_filter <- expr[row.names(expr) %in% names(filter_mask),]

  if(verbose == TRUE) print(paste0(length(mean_cor) - length(filter_mask), " Filtered at ", cor_threshold))

  res <- list()

  if("expr" %in% return_output){
    res[["expr"]] <- expr_filter
  }

  if("metric" %in% return_output){
    dist_res <- as.dist(1 - cor(t(expr_filter), method = cor_method))
    hclust_res <- hclust(dist_res, method = clus_method)

    sil_res <- best_k_labels(hclust_res, dist_res, expr_filter, ...)
    best_k <- sil_res[["best_k"]]

    var_res <- calculate_cluster_var(expr_filter, sil_res$labels)

    #mean_abs_cor = mean(filter_mask),
    res[["metric"]] <- list(
      threshold = cor_threshold,
      silhouette = sil_res$curve[sil_res$curve$k == sil_res$best_k,]$mean_sil,
      sd = mean(var_res[["SD"]][var_res[["SD"]] > 0]),
      var = mean(var_res[["Var"]][var_res[["Var"]] > 0]),
      remaining = length(filter_mask),
      n_cluster = best_k,
      singlet = sil_res$singlet,
      split_k = sil_res[["split_k"]]
    )
    if(return_second_cluster == TRUE){
      res[["metric"]][["n_cluster_second"]] <- sil_res$gain_k
      res[["metric"]][["silhouette_second"]] <- sil_res$curve[sil_res$curve$k == sil_res$gain_k,]$mean_sil
    }

  }

  if(length(return_output) == 1){
    res <- list_c(res)
  }
  return(res)
}

#' Apply High Pass Filter Iterations and Return Metrics
#'
#' @param input_data A data.frame object, row names as observations (i.e. genes) and column names as variables (i.e. time points)
#' @param threshold_sequence The sequence of cutoffs to iterate through.
#' @param return_output Input one of "expr", "metric", or both by default. Returns the filtered data or the cluster metrics, respectively.
#' @param iterate_data TRUE/FALSE. Determines whether to overwrite the expression data at each iterated threshold.
#' @param verbose Print progress messages
#' @param ... Additional parameters to be passed on to [totoInferno::cor_filter]
#' @inheritDotParams cor_filter -expr -cor_threshold -return_output -verbose
#' @inheritDotParams best_k_labels k_min k_max
#'
#' @export

cor_filter_iterate <- function(

  input_data,
  threshold_sequence,
  return_output = c("expr", "metric"),
  iterate_data = FALSE,
  verbose     = FALSE,

  ...

) {

  expr <- as.matrix(input_data)

  # Run iterations of scanned parameters

  res_list <- list()

  for (i in seq_along(threshold_sequence)) {
    res_it <- tryCatch(
      cor_filter(
        expr,
        cor_threshold = threshold_sequence[i],
        return_output = c("expr", "metric"),
        ...
      ),
      error = function(e) NULL
    )
    if (!is.null(res_it)){
      if(iterate_data == TRUE){expr <- res_it[["expr"]]}
      res_list[[i]] <- as.data.frame(res_it[["metric"]])
    }
  }


  res <- list()
  # First run for comparison
  if("metric" %in% return_output){
    res[["metric"]] <- list_c(res_list)
    scan_first <- cor_filter(
      input_data,
      cor_threshold      = 0,
      return_output = "metric",
      ...
    )
    scan_append <- as.data.frame(scan_first)
  }
  if("expr" %in% return_output){
    res[["expr"]] <- expr
  }

  if(length(return_output) == 1){
    res <- list_c(res)
  }

  return(res)
}
