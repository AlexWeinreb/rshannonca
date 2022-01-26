#' Title
#'
#' @param seu_object Seurat object
#' @param assay Assay
#' @param npcs Desired dimensionality of the reduction.
#' @param n_iter Number of iterations of SCA. More iterations usually strengthens signal, stabilizing around 3-5
#' @param nbhd_size Size of neighborhoods used to assess the local expression of a gene. Should be smaller than the smallest subpopulation.
#' @param model Model used to test for local enrichment of genes, used to compute information scores. One of [“wilcoxon”,”binomial”,”ttest”], default “wilcoxon” (recommended).
#' @param metric Metric used to compute k-nearest neighbor graphs for SCA score computation. Default “euclidean”. See sklearn.neighbors.DistanceMetric for list of choices.
#' @param verbose If True, print progress.
#' @param n_tests Effective number of independent genes per cell, use for FWER multiple testing correction. Set to “auto” to automatically determine by bootstrapping.
#' @param reduction.name Name to store dimensional reduction under in the Seurat object.
#' @param reduction.key Dimensional reduction key, specifies the string before the number for the dimension names.
#' @param ... Additional arguments passed to the chosen scorer.
#'
#' @return Returns a Seurat object containing an SCA representation.
#' @export
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- RunSCA(pbmc_small)
#' DimPlot(pbmc_small, reduction = "sca", dims = 1:2)
RunSCA <- function (seu_object, assay = NULL, npcs = 50L, n_iter = 5L,
                    nbhd_size = 15L, model = c("wilcoxon", "binomial", "ttest"),
                    metric = "euclidean", verbose = FALSE, n_tests = "auto",
                    reduction.name = "sca", reduction.key = "SC_", ...)
{

  assay <- if(!is.null(assay)) assay else Seurat::DefaultAssay(object = seu_object)
  assay_object <- Seurat::GetAssay(object = seu_object, assay = assay)
  object <- Seurat:::PrepDR(object = assay_object, verbose = verbose)

  npcs <- min(npcs, nrow(x = object)-1)

  model <- match.arg(model)


  sca <- reticulate::import("shannonca")

  sca.results <- sca$reduce(t(object),
                            n_comps = as.integer(npcs),
                            iters = as.integer(n_iter),
                            nbhd_size = as.integer(nbhd_size),
                            model = model,
                            metric = metric,
                            verbose = verbose,
                            n_tests = n_tests,
                            kwargs = ...)



  rownames(x = sca.results) <- rownames(x = t(object))
  colnames(x = sca.results) <- paste0(reduction.key,
                                      1:npcs)

  reduction.data <- Seurat::CreateDimReducObject(embeddings = sca.results,
                                         assay = assay,
                                         key = reduction.key,
                                         global = TRUE)


  seu_object[[reduction.name]] <- reduction.data
  seu_object <- Seurat::LogSeuratCommand(object = seu_object)
  return(seu_object)
}
