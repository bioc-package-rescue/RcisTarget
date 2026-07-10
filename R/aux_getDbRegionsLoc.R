#' @title Gets the region location for the given database IDs.
#' @description
#' Gets the region location based on the region ID.
#' Only needed for drosophila\bold{drosophila} (fruit fly) regions.
#' 
#' For \bold{human/mouse} the region locations are stored in a separate object 
#' (i.e. \code{}). This function is not needed.
#' @param featherFilePath Path to the rankings database
#' @param spltChr Character(s) used to split the prefix from the region location.
#' The default is used for current Drosophila versions. Use NULL to skip.
#' @param indexCol Column name including the motif/feature IDs. By default it will use the first column.
#' @return The region locations in a GRanges object, with the original region ID as name.
#'
#' See the package vignette for examples and more details:
#' \code{vignette("RcisTarget")}
#'
#' @examples
#' library(arrow)
#' # Create a dummy data frame with drosophila-style region column names
#' df <- data.frame(
#'   features = c("motif1", "motif2"),
#'   `motif1__chr3L:1000-2000` = c(1, 2),
#'   `motif2__chr3L:3000-4000` = c(3, 4),
#'   check.names = FALSE
#' )
#' 
#' # Write to a temporary feather file
#' tmpFile <- tempfile(fileext = ".feather")
#' write_feather(df, tmpFile)
#' 
#' # Get the region locations
#' getDbRegionsLoc(tmpFile)
#' 
#' # Clean up
#' file.remove(tmpFile)
#' @rdname getDbRegionsLoc
#' @importFrom GenomeInfoDb sortSeqlevels seqlevels 
#' @importFrom GenomicRanges GRanges
#' @export
getDbRegionsLoc <- function(featherFilePath, spltChr="__", indexCol=NULL)
{
  featherFilePath <- path.expand(featherFilePath)
  
  dbCols <- getColumnNames(featherFilePath)
  indexCol <- .getIndexCol(dbCols, indexCol=indexCol, verbose=FALSE)
  
  dbCols <- setdiff(dbCols, indexCol); length(dbCols)
  dbCols <- setNames(dbCols, dbCols)
  
  if(!is.null(spltChr)) dbCols <- sapply(strsplit(dbCols, spltChr), function(x) x[[2]])
  dbRegions <- GenomicRanges::GRanges(unname(dbCols), name=names(dbCols))
  dbRegions <- GenomeInfoDb::sortSeqlevels(dbRegions)
  return(dbRegions)
}

