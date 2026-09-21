#' Cluster and Filter Data Using Silhouette Scores
#'
#' @param expr A data.frame or matrix object, row names as observations (i.e. genes) and column names as variables (i.e. time points)
#' @param k_min Minimum number of clusters to test
#' @param k_max Maximum number of clusters to test
#' @param sil_threshold Threshold score, removing genes which score below this value
#' @param cor_method Which method to use for correlation. Pearson is strongly recommended, but option for Kendall or Spearman is present.
#' @param clus_method Which method to use for hierarchical clustering. Complete method allows for non-Euclidian distances, since we are taking pairwise correlations between objects.
#' @param condition_name The name of your condition, which will be paired with the output.
#'
#' @return A list object containing various data metrics, filter metrics, and returned data.
#'
#' @export

silhouette_filter <- function(expr,
                              k_min          = 2L,
                              k_max          = 10L,
                              sil_threshold  = 0,
                              cor_method     = "pearson",
                              clus_method    = "complete",
                              condition_name = "condition") {

  expr <- as.matrix(expr) #turn into matrix for computational efficiency
  if (is.null(rownames(expr))) stop("Input data must have gene IDs as rownames.")

  ## -- Step 1: distance and clustering ----
  d1 <- as.dist(1 - cor(t(expr), method = cor_method))
  hc1 <- hclust(d1, method = clus_method)

  ## -- Step 2: First cluster assignment (intended for silhouette curve) ----

  # The initial silhouette grouping genes inappropriately is a major source of error
  # Too little or (especially) too many groups affects how this pipeline determines "well fit" clusters

  # If the output is chopping a lot of genes at a low threshold, such as -0.5, the input data is being clustered poorly with the 'complete' method. Use a mean correlation filter to remove some noise and assist the complete clustering. Or (not recommended) use another clustering method.

  cut1 <- best_k_labels(hc1, d1, k_min = k_min, k_max = k_max) #take best K
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

  cut2 <- best_k_labels(hc2, d2, k_min = k_min, k_max = k_max)
  labels_final <- setNames(cut2$labels, row.names(expr_keep))

  n_clusters_final <- length(unique(labels_final))

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
  #within_sd_per_cluster <- rowMeans(sd_profiles) # sd_profiles is time x time matrix of SD, not the gene expression data
  var_clusters <- data.frame(Cluster = cut2$labels,
                             Gene = names(cut2$labels))
  find_var_sd     <- calculate_cluster_var(expr_keep, var_clusters) #overall score for reporting


  ## -- Step 7: Create list of metrics based on how this pipeline operated ----
  # No one metric is the definitive 'grade'
  # multiple metrics give insight into data refinement and changes over threshold iterations

  metrics <- list(
    mean_silhouette_final = mean_sil_final,
    silhouette_per_cluster = sil_per_cluster,
    within_cluster_sd_overall = find_var_sd$SD,
    within_cluster_var_overall = find_var_sd$Var,
    genes_removed = n_removed,
    n_clusters_final = n_clusters_final,
    n_genes_in  = nrow(expr),
    n_genes_out = nrow(expr_keep)
  )

  ## For scan plotting: silhouette, SD, and 'percent genes remaining' scale naturally to 0-1 scale

  ## -- Step 8: Output list of values plus metrics (Step 7) ----

  list(
    condition      = condition_name,
    expr_keep      = expr_keep,       # Genes remaining
    genes_keep     = row.names(expr_keep),
    dist_final     = d2,
    hclust_final   = hc2,
    labels_final   = labels_final,    # Named survivors
    sil_curve_p1   = cut1[["curve"]],      # Silhouette curve via silhouette optimized (curve) method
    sil_curve_final= mean_sil_final,  # Silhouette curve via dynamic tree method
    metrics        = metrics,
    params = list(k_min = k_min,
                  k_max = k_max,
                  sil_threshold = sil_threshold,
                  cor_method = cor_method)
  )
}


#' Find Silhouette-Optimized Clusters
#'
#' This helper function assists with finding the optimal number of clusters using the silhouette method.
#'
#' @param hc An object produced by the hclust() function. Hierarchical clustered data.
#' @param dist_obj An object produced by the as.dist() function. A one-dimensional vector of values, calculated as distance between pairwise correlations
#' @param k_min Minimum number of clusters to test
#' @param k_max Maximum number of clusters to test
#'
#' @return A list of labels assigned, the optimal k, and the silhouette curve.
#'
#' @export

best_k_labels <- function(hc,
                          dist_obj,
                          k_min = 2L,
                          k_max = 10L) {

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

  return(list( #return our data
    labels   = cutree(hc, k = best_k),
    best_k   = best_k,
    curve    = data.frame(
      k       = kk,
      mean_sil = curve
    )
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

calculate_cluster_var <- function(input_counts,
                                  input_clusters) {

  # For each cluster, isolate and find the values
  cluster_count <- unique(input_clusters)

  cluster_res <- map_dfr(cluster_count, \(cluster_i) {
    pickCluster <- input_clusters[input_clusters == cluster_i]
    pickGenes <- input_counts[names(pickCluster),]

    clusterVar <- mean(apply(pickGenes, 2, var))
    clusterSD <- mean(apply(pickGenes, 2, sd))

    clusterInfo <- data.frame( # Combine for graphing
      Var  = clusterVar,
      SD = clusterSD,
      Cluster = as.character(cluster_i),
      row.names = NULL
    )
    return(clusterInfo)
  })
  return(cluster_res)

}

#' Apply Silhouette Filter Iterations and Return Data.
#'
#' This function is similar to silhouette_filter_scan but removes many redundancies for reporting how the pipeline works through the data. This increases speed and reduces memory burden. The intended use is with parameters determined through silhouette_filter_scan, then using this function to more efficiently apply the filter and return the data.
#'
#' @param expr A data.frame or matrix object, row names as observations (i.e. genes) and column names as variables (i.e. time points)
#' @param condition_name The name of your condition, which will be paired with the output.
#' @param initial_threshold Which silhouette threshold to begin iterations at. -0.5 is a safe start. Starting too low will massively waste time, starting too high will lead to 'true negative' filters where genes were placed in incorrect clusters via the hierarchical clustering method. The pipeline becomes faster as more genes are removed, so a low initial_threshold is equal to the longest computation time and provided the least (sometimes zero) filtering.
#' @param final_threshold Which silhouette threshold to end at. The maximum value is 1, so thresholding too close to that will just remove nearly all the genes and usually create errors in doing so. 0.8-0.9 recommended, since clustering a little too far will help visualize the 'drop' in quality gained by removing genes. Lowering the top threshold generally doesn't save too much time, given that the pipeline works faster as the number of genes decreases. This requires input, given that the final threshold should be determined through use of the silhouette_filter_scan function.
#' @param initial_clusters A data.frame object with two columns: observations (i.e. genes) and clusters. To allow for different naming conventions, the ORDER is important: initial_clusters column 1 == character class (observations), initial_clusters column 2 == numeric class (cluster assignments). If the input data is particularly stubborn to prepare for this pipeline, this is the option to manually set clusters using settings which aren't available or aren't conducive to the silhouette method.
#' @param stepsize The step size for increasing the silhouette threshold during iterations. Keep in mind our full threshold range is only -1 to 1. 0.01 is recommended for allowing some specificity to small changes while still scanning through all iterations quickly enough. Going smaller generally will not yield more information, but will drastically increase the total computational cost. Small iterations will not remove the inherent behavior of cutting many genes with one step. Going larger may be useful for quickly scanning large data sets or to find a good initial threshold, but may cut too many genes at once for the final test.
#' @param cor_method Which method to use for correlation. Pearson is strongly recommended, but option for Kendall or Spearman is present.
#' @param clus_method Which method to use for hierarchical clustering. Complete method allows for non-Euclidian distances, since we are taking pairwise correlations between objects.
#' @param k_min Minimum number of clusters to test
#' @param k_max Maximum number of clusters to test
#' @param verbose Prints the threshold at each step. Useful for monitoring progress.
#'
#' @return A list object containing various data metrics, filter metrics, and returned data.
#'
#' @export

silhouette_filter_return <- function(expr,
                                     condition_name    = "condition",
                                     initial_threshold = -0.5,
                                     final_threshold   = 0.5,
                                     initial_clusters  = NULL,
                                     stepsize          = 0.01,
                                     cor_method        = "pearson",
                                     clus_method       = "complete",
                                     k_min             = 2L,
                                     k_max             = 10L,
                                     verbose           = FALSE) {
  thresholds <- seq(initial_threshold, final_threshold, by = stepsize) #create threshold sequence

  ## If initial clusters are designated
  # This code is the same as used in the silhouette_filter, but makes use of assigned clusters here
  if(length(initial_clusters) > 1){ #if initial_clusters have been inputted, use them

    dmat_initial <- as.dist(1 - cor(t(expr), method = cor_method))
    dmat_initial <- as.matrix(dmat_initial) #distance and correlation

    sil_tab <- data.frame(
      gene       = initial_clusters[,1],
      cluster_p1 = as.numeric(initial_clusters[,2]),
      sil_width  = NA_real_,
      stringsAsFactors = FALSE
    )
    if (length(unique(initial_clusters)) >= 2) {
      sil <- silhouette(as.numeric(initial_clusters[,2]), dmatrix = dmat_initial)
      sil_tab$sil_width <- sil[, "sil_width"]
    } else {
      warning("[", condition_name, "] Pass 1 produced < 2 clusters; ",
              "silhouette undefined.")
    }

    keep_mask     <- !is.na(sil_tab$sil_width) &
      (sil_tab$sil_width > initial_threshold) #remove the input cluster genes based off of initial threshold
    genes_keep    <- rownames(expr)[keep_mask]
    genes_removed <- setdiff(rownames(expr), genes_keep)
    n_removed     <- length(genes_removed)

    expr_keep <- expr[genes_keep, , drop = FALSE]
    newexpr <- expr_keep
  }
  else{ #this runs if there are not intitial clusters present
    newexpr <- expr
  }

  ## Run the iterative scan
  for(th in thresholds){
    n_removed <- 1 #set this at the start to allow for the next step

    while(n_removed > 0){ #repeat filters until the filter is not removing genes. Gives data a chance to stabilize, especially with larger jumps
        if(verbose == TRUE){print(th)}
        ## run the full pipeline at this threshold; capture failures as NA rows
        res <- tryCatch(
          silhouette_filter(newexpr, sil_threshold = th,
                            condition_name = condition_name,
                            clus_method = clus_method, cor_method = cor_method,
                            k_min = k_min, k_max = k_max),
          error = function(e) NULL
        )

        if(!is.null(res)){
          m <- res$metrics
          n_removed <- m$genes_removed
          newexpr <- newexpr %>% filter(row.names(newexpr) %in% res$genes_keep)
        }

  }}
  return(res)
}

#' Apply Silhouette Filter Iterations and Report Quality Metrics
#'
#' This function applies the silhouette_filter function across a designated sequence of iterations. The outputs are various data quality measurements to assess an appropriate start and end for iterative filtering. Currently, graphing the data metrics is the most convenient way to determine appropriate parameters for use in silhouette_filter_return. Additional methods to automatically determine the optimal final threshold are provided elsewhere.
#'
#' @param expr A data.frame or matrix object, row names as observations (i.e. genes) and column names as variables (i.e. time points)
#' @param condition_name The name of your condition, which will be paired with the output.
#' @param initial_threshold Which silhouette threshold to begin iterations at. -0.5 is a safe start. Starting too low will massively waste time, starting too high will lead to 'true negative' filters where genes were placed in incorrect clusters via the hierarchical clustering method. The pipeline becomes faster as more genes are removed, so a low initial_threshold is equal to the longest computation time and provided the least (sometimes zero) filtering.
#' @param final_threshold Which silhouette threshold to end at. The maximum value is 1, so thresholding too close to that will just remove nearly all the genes and usually create errors in doing so. 0.8-0.9 recommended, since clustering a little too far will help visualize the 'drop' in quality gained by removing genes. Lowering the top threshold generally doesn't save too much time, given that the pipeline works faster as the number of genes decreases.
#' @param initial_clusters A data.frame object with two columns: observations (i.e. genes) and clusters. To allow for different naming conventions, the ORDER is important: initial_clusters column 1 == character class (observations), initial_clusters column 2 == numeric class (cluster assignments). If the input data is particularly stubborn to prepare for this pipeline, this is the option to manually set clusters using settings which aren't available or aren't conducive to the silhouette method.
#' @param stepsize The step size for increasing the silhouette threshold during iterations. Keep in mind our full threshold range is only -1 to 1. 0.01 is recommended for allowing some specificity to small changes while still scanning through all iterations quickly enough. Going smaller generally will not yield more information, but will drastically increase the total computational cost. Small iterations will not remove the inherent behavior of cutting many genes with one step. Going larger may be useful for quickly scanning large data sets or to find a good initial threshold, but may cut too many genes at once for the final test.
#' @param cor_method Which method to use for correlation. Pearson is strongly recommended, but option for Kendall or Spearman is present.
#' @param clus_method Which method to use for hierarchical clustering. Complete method allows for non-Euclidian distances, since we are taking pairwise correlations between objects.
#' @param k_min Minimum number of clusters to test
#' @param k_max Maximum number of clusters to test
#' @param verbose Prints the threshold and genes removed at the step. Useful for monitoring progress and speed, or quickly visualizing unexpected results (i.e. too many genes removed early) and halt the function and save time.
#'
#' @return A data.frame containing data metrics as the filter threshold increases through iterations.
#'
#' @export

silhouette_filter_scan <- function(expr,
                           initial_threshold = -0.5,
                           final_threshold   = 0.9,
                           stepsize          = 0.01,
                           condition_name    = "condition",
                           initial_clusters  = NULL,
                           cor_method        = "pearson",
                           clus_method       = "complete",
                           k_min             = 2L,
                           k_max             = 10L,
                           verbose           = FALSE) {


  ## If initial clusters are designated
  # This code is the same as used in the silhouette_filter, but makes use of assigned clusters here
  if(length(initial_clusters) > 1){ #if initial_clusters have been inputted, use them

    dmat_initial <- as.dist(1 - cor(t(expr), method = cor_method))
    dmat_initial <- as.matrix(dmat_initial) #distance and correlation

    sil_tab <- data.frame(
      gene       = initial_clusters[,1],
      cluster_p1 = as.numeric(initial_clusters[,2]),
      sil_width  = NA_real_,
      stringsAsFactors = FALSE
    )
    if (length(unique(initial_clusters)) >= 2) {
      sil <- silhouette(as.numeric(initial_clusters[,2]), dmatrix = dmat_initial)
      sil_tab$sil_width <- sil[, "sil_width"]
    } else {
      warning("[", condition_name, "] Pass 1 produced < 2 clusters; ",
              "silhouette undefined.")
    }

    keep_mask     <- !is.na(sil_tab$sil_width) &
      (sil_tab$sil_width > initial_threshold) #remove the input cluster genes based off of initial threshold
    genes_keep    <- rownames(expr)[keep_mask]
    genes_removed <- setdiff(rownames(expr), genes_keep)
    n_removed     <- length(genes_removed)

    expr_keep <- expr[genes_keep, , drop = FALSE]
    newexpr <- expr_keep
  }
  else{ #this runs if there are not intitial clusters present
    newexpr <- expr
  }

  thresholds <- seq(initial_threshold, final_threshold, by = stepsize)
  scan <- data.frame(
    threshold                 = NA_integer_,
    n_gene_in                 = NA_integer_,
    genes_removed             = NA_integer_,
    n_genes_out               = NA_integer_,
    mean_silhouette_final     = NA_real_,
    within_cluster_sd_overall = NA_real_,
    within_cluster_var_overall = NA_real_,
    n_clusters_filtered       = NA_integer_,
    n_clusters_final          = NA_integer_,
    between_cluster_separation= NA_integer_,
    stringsAsFactors = FALSE
  )
  n_removed <- 1
  for(th in thresholds){
    while(n_removed > 0){
      ## run the full pipeline at this threshold; capture failures as NA rows
      res <- tryCatch(
        silhouette_filter(newexpr, sil_threshold = th,
                          condition_name = condition_name, k_min = k_min, k_max = k_max,
                          clus_method = clus_method, cor_method = cor_method),
        error = function(e) NULL
      )
      if (is.null(res)) { #return NA during iterations which created an error
        rbind(scan, data.frame(
          threshold                 = th,
          n_gene_in                 = NA_integer_,
          genes_removed             = NA_integer_,
          n_genes_out               = NA_integer_,
          mean_silhouette_final     = NA_real_,
          within_cluster_sd_overall = NA_real_,
          within_cluster_var_overall = NA_real_,
          n_clusters_filtered       = NA_integer_,
          n_clusters_final          = NA_integer_,
          between_cluster_separation= NA_integer_,
          stringsAsFactors = FALSE
        ))
          if(verbose == TRUE){print(paste0("Threshold: ",th,", Error returned"))}
        n_removed <- 0
      }else{ #return metrics in iterations which succeeded
        m <- res$metrics
        newexpr <- newexpr %>% filter(row.names(newexpr) %in% res$genes_keep)
        scan <- rbind(scan,
                      data.frame(
                        threshold                  = th,
                        n_gene_in                  = m$n_genes_in,
                        genes_removed              = m$genes_removed,
                        n_genes_out                = m$n_genes_out,
                        mean_silhouette_final      = m$mean_silhouette_final,
                        within_cluster_sd_overall  = m$within_cluster_sd_overall,
                        within_cluster_var_overall = m$within_cluster_var_overall,
                        n_clusters_filtered        = m$n_clusters_pass1,
                        n_clusters_final           = m$n_clusters_final,
                        between_cluster_separation = m$between_cluster_separation,
                        stringsAsFactors = FALSE
                      ))
        if(verbose == TRUE){print(paste0("Threshold: ",th,", Removed: ",m$genes_removed, " observations"))}
        n_removed <- m$genes_removed
      }}
    n_removed <- 1
  }

  return(scan)
}

#' Apply High Pass Filter to Remove Isolated Observations
#'
#' This function uses the average correlation of observations (i.e. genes) to remove those which have the smallest correlations to the system. This increases the cluster separation when using the 'complete' hierarchical clustering method. Increased cluster separation allows for the silhouette_filter function to handle input data more predictably. Intended for use after finding an appropriate cutoff using the cor_filter_scan function.
#'
#' @param input_data A data.frame object, row names as observations (i.e. genes) and column names as variables (i.e. time points)
#' @param cor_method Which method to use for correlation. Pearson is strongly recommended, but option for Kendall or Spearman is present.
#' @param clus_method Which method to use for hierarchical clustering. Complete method allows for non-Euclidian distances, since we are taking pairwise correlations between objects.
#' @param cutoff The cutoff to remove genes with the lowest average correlation
#' @param k_min Minimum number of clusters to test
#' @param k_max Maximum number of clusters to test
#' @param verbose Print progress messages
#'
#' @return A list object containing two data.frames: The filtered version of the input data.frame, and a data.frame containing the genes (observations) and cluster assignments after filtering.
#'
#' @export

cor_filter_return <- function(input_data,
                              cor_method  = "pearson",
                              clus_method = "complete",
                              cutoff      = 0,
                              k_min       = 2L,
                              k_max       = 10L,
                              verbose     = FALSE) {
  ## Initial processing of raw df
  inputCor <- abs(cor(t(input_data), method = cor_method))
  meanCor <- apply(inputCor, 1, mean) #rowwise mean of correlations

  meanFilter <- meanCor[meanCor >= cutoff] #remove genes below the cutoff

  if(verbose == TRUE) print(paste0(length(meanCor) - length(meanFilter), " Filtered at ", cutoff))

  filterData <- input_data %>% filter(row.names(input_data) %in% names(meanFilter))
  filterCor <- cor(t(filterData), method = cor_method)

  filterDist <- as.dist(1 - filterCor)
  filterHclust <- hclust(filterDist, method = clus_method)

  # Report silhouette on Curve method, since this is what we plug in to the filter
  filterCurve <- best_k_labels(filterHclust, filterDist, k_min = k_min, k_max = k_max)
  filterLabels <- data.frame(Gene = names(filterCurve$labels), Cluster = filterCurve$labels)

  return(list(filterData = filterData, filterLabels = filterLabels))
}

#' Apply High Pass Filter and Return Filter Information
#'
#' This function uses the average correlation of observations (i.e. genes) to remove those which have the smallest correlations to the system. This increases the cluster separation when using the 'complete' hierarchical clustering method. Increased cluster separation allows for the silhouette_filter function to handle input data more predictably. This function specifically returns various metrics for use in the cor_filter_scan function to determine the optimal initial high pass filter.
#'
#' @param input_data A data.frame object, row names as observations (i.e. genes) and column names as variables (i.e. time points)
#' @param cor_method Which method to use for correlation. Pearson is strongly recommended, but option for Kendall or Spearman is present.
#' @param clus_method Which method to use for hierarchical clustering. Complete method allows for non-Euclidian distances, since we are taking pairwise correlations between objects.
#' @param cutoff The cutoff to remove genes with the lowest average correlation
#' @param k_min Minimum number of clusters to test
#' @param k_max Maximum number of clusters to test
#' @param verbose Print progress messages
#'
#' @return A data.frame object containing Standard deviation, Variance, and Silhouette curve.
#'
#' @export

cor_filter <- function(input_data,
                       cor_method  = "pearson",
                       clus_method = "complete",
                       cutoff      = 0,
                       k_min       = 2L,
                       k_max       = 10L,
                       verbose     = FALSE) {
  # Initial processing of raw df
  inputCor <- abs(cor(t(input_data), method = cor_method))
  meanCor <- apply(inputCor, 1, mean)

  meanFilter <- meanCor[meanCor >= cutoff]
  if(verbose == TRUE){print(paste0(length(meanCor) - length(meanFilter), " Filtered at ", cutoff))}

  filterData <- input_data %>% filter(row.names(input_data) %in% names(meanFilter))
  filterCor <- cor(t(filterData), method = cor_method)

  filterDist <- as.dist(1 - filterCor)
  filterHclust <- hclust(filterDist, method = clus_method)

  # Report silhouette on Curve method, since this is what we plug in to the filter
  filterCurve <- best_k_labels(filterHclust, filterDist, k_min = k_min, k_max = k_max)

  # Report Var and SD
  clusterCurve <- data.frame(Cluster = filterCurve$labels, Gene = names(filterCurve$labels))
  clusterVar <- calculate_cluster_var(filterData, clusterCurve)
  clusterVar <- clusterVar %>% filter(.data$Group == "Mean")

  metric_data <- data.frame(
    Curve = filterCurve$curve[as.numeric(filterCurve$best_k)-1,-3],
    Curve.SD = mean(clusterVar$SD),
    Curve.Var = mean(clusterVar$Var),
    n_remaining = nrow(filterCor)
  )

  return(metric_data)
}

#' Apply High Pass Filter Iterations and Return Metrics
#'
#' @param input_data A data.frame object, row names as observations (i.e. genes) and column names as variables (i.e. time points)
#' @param cutoff_sequence The sequence of cutoffs to iterate through.
#' @param cor_method Which method to use for correlation. Pearson is strongly recommended, but option for Kendall or Spearman is present.
#' @param clus_method Which method to use for hierarchical clustering. Complete method allows for non-Euclidian distances, since we are taking pairwise correlations between objects.
#' @param k_min Minimum number of clusters to test
#' @param k_max Maximum number of clusters to test
#' @param verbose Print progress messages
#'
#' @export

cor_filter_scan <- function(input_data,
                            cutoff_sequence,
                            cor_method  = "pearson",
                            clus_method = "complete",
                            k_min       = 2L,
                            k_max       = 10L,
                            verbose     = FALSE) {
  scan_append <- data.frame()
  # First run for comparison
  scan_first <- cor_filter(input_data,
                           cutoff      = 0,
                           cor_method  = cor_method,
                           clus_method = clus_method,
                           k_min       = k_min,
                           k_max       = k_max)
  scan_append <- rbind(scan_append, data.frame(
    Cutoff = 0,
    scan_first
  ))

  # Run iterations of scanned parameters
  for(ct in cutoff_sequence){ # along each cutoff
    scan_cutoff <- tryCatch(cor_filter(input_data,
                                       cutoff = ct,
                                       cor_method = cor_method,
                                       clus_method = clus_method),
                            error = function(e) NULL)
    if(is.null(scan_cutoff)){
      if(verbose == TRUE){print(paste0(ct, " Cutoff returns error"))}

      scan_append <- rbind(scan_append, data.frame(
        Cutoff = ct,
        Curve.k = NA_integer_,
        Curve.mean_sil = NA_integer_,
        Curve.SD = NA_integer_,
        Curve.Var = NA_integer_,
        n_remaining = NA_integer_
      ))
    }
    else{
      if(verbose == TRUE){print(paste0(ct, " Cutoff, Removed: ",
                                       as.numeric(nrow(input_data)) - as.numeric(scan_cutoff$n_remaining)))}
      scan_append <- rbind(scan_append, data.frame(
        Cutoff = ct,
        scan_cutoff
      ))
    }
  }
  return(scan_append)
}

#' Find Metrics on Best Cluster
#'
#' @param input_data A data.frame object, row.names as observations, column names as variables
#' @param cor_method String for which method to pass on to cor(), i.e. "pearson".
#' @param clus_method String for which method to pass on to hclust(), i.e. "complete".
#' @param k_min Minimum number of clusters to test
#' @param k_max Maximum number of clusters to test
#'
#' @export

cluster_metric <- function(input_data,
                           cor_method  = "pearson",
                           clus_method = "complete",
                           k_min       = 2L,
                           k_max       = 10L) {
  input_mat <- as.matrix(input_data)

  cor_res <- cor(t(input_mat), method = cor_method)
  cor_hclust <- hclust(as.dist(1 - cor_res), method = clus_method)

}
