#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end

#' @importFrom rlang .data
#' @importFrom dplyr count filter setdiff
#' @importFrom cluster silhouette
#' @importFrom stats cor as.dist hclust cutree var sd
#' @importFrom utils globalVariables
NULL

utils::globalVariables(c("Cluster", "Group"))
