#' rcTrain
#'
#' Train a rcissus model
#'
#' @param x input filenames
#' @return model is output
#' @author Avants BB
#' @seealso \code{\link[ANTsRCore]{ripmmarcPop}} \url{https://antsx.github.io/ANTsRCore/reference/ripmmarcPop.html}
#' @examples
#' pop = getANTsRData( "population" ) # list of example images
#' popmasks = list( )
#' for ( i in 1:length( pop ) )
#'   popmasks[[ i ]] = getMask( pop[[ i ]] )
#' rp = ripmmarcPop( pop, popmasks, patchRadius=3,
#'   meanCenter = TRUE, patchSamples=1000 )
#' nv = 15
#' rippedTest <- ripmmarc( pop[[3]], popmasks[[3]], patchRadius = 3,
#'   evecBasis = rp$basisMat[1:nv,], patchVarEx = nv, meanCenter = TRUE,
#'   canonicalFrame = rp$canonicalFrame, regressProjections = TRUE )
#' mm = makeImage( popmasks[[3]], rippedTest$evecCoeffs[,1] )
#'
#' @export rcTrain
#' @import ANTsR
#' @import methods
#' @import h2o
rcTrain <- function( x ) {
  return( NA )
}
