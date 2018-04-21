#' rcBasis
#'
#' Learn a rcissus image basis
#'
#' @param x input filenames or image list
#' @param patchRadius patch radius, integer value
#' @return basis is output
#' @author Avants BB
#' @seealso \code{\link[ANTsRCore]{ripmmarcPop}} \url{https://antsx.github.io/ANTsRCore/reference/ripmmarcPop.html}
#' @examples
#'
#' pop = getANTsRData( "population" ) # list of example images
#' imgBases = rcBasis( pop )
#'
#' @export rcBasis
#' @importFrom utils write.csv
#' @import ANTsR
#' @importFrom ANTsRCore antsImageRead
#' @importFrom ANTsRCore antsImageWrite
#' @importFrom ANTsRCore imageFileNames2ImageList
#' @importFrom ANTsRCore getMask
#' @importFrom ANTsRCore ripmmarc
#' @importFrom ANTsRCore ripmmarcPop
#' @importFrom ANTsRCore randomMask
#' @import methods
#' @import h2o
rcBasis <- function( x, patchRadius = 3  ) {
  maskType = 'dense'
  if ( ! ( maskType %in% c("sparse","dense") ) ) stop("pass dense or sparse as maskType" )
  if ( length( x ) < 1 ) stop( "pass in a non-zero length input" )

  if ( class( x[[1]] )[1] == 'character' )
    isFilename = TRUE else isFilename = FALSE

  # build a training dataset by sampling each image and binding together
  # the output of the matrices
  nsam = 1000
  popmasks = list()
  for ( i in 1:length( x ) ) {
    if ( isFilename ) img = ANTsRCore::antsImageRead( x[ i ] ) else img = x[[ i ]]
    if ( maskType == 'dense' )
      popmasks[[ i ]] = ANTsRCore::getMask( img, cleanup = 0  )
    if ( maskType == 'sparse' ) {
      temp = ANTsRCore::getMask( img, cleanup = 0 )
      popmasks[[ i ]] = ANTsRCore::randomMask( temp, nsam, perLabel = TRUE ) * temp
      }
    }
  if ( isFilename ) {
    rp = ANTsRCore::ripmmarcPop( ANTsRCore::imageFileNames2ImageList( x ),
      popmasks, patchRadius=patchRadius, meanCenter = FALSE, patchSamples=nsam )
    } else {
    rp = ANTsRCore::ripmmarcPop( x,
      popmasks, patchRadius=patchRadius, meanCenter = FALSE, patchSamples=nsam )
    }
  return( rp )
}




#' rcTrainingMatrix
#'
#' Build a rcissus training matrix
#'
#' @param y input filenames or image list defining ground truth
#' @param x input filenames or image list defining one training modality
#' @param masks input filenames or image list defining training modality masks
#' @param rcb input image basis from \code{\link{rcBasis}}
#' @param patchRadius patch radius, integer value
#' @param nsamples number of samples per image
#' @param seeds random seeds for reproducibility across training modalities
#' @return training matrix and vector (ground truth) is output in list
#' @author Avants BB
#' @seealso \code{\link[ANTsRCore]{ripmmarcPop}} \url{https://antsx.github.io/ANTsRCore/reference/ripmmarcPop.html}
#'
#' @export rcTrainingMatrix
rcTrainingMatrix <- function( y, x, masks, rcb, patchRadius = 3, nsamples = 1000, seeds ) {
  if ( length( x ) < 1 ) stop( "pass in a non-zero length input" )
  if ( length( x ) != length( y ) ) stop( "length of y should equal length of x" )

  if ( class( x[[1]] )[1] == 'character' )
    isFilename = TRUE else isFilename = FALSE

  if ( missing( seeds ) )
    seeds = sample( 1:10000 )[1:length(x)]

  if ( length( x ) != length( seeds ) ) stop( "length of seeds should equal length of x" )

  if (
    ( class( y[[1]] )[1] == 'character' & !isFilename ) |
    ( class( y[[1]] )[1] != 'character' &  isFilename )  )
    stop(" both x and y should be either filenames or image lists")

  if (
    ( class( masks[[1]] )[1] == 'character' & !isFilename ) |
    ( class( masks[[1]] )[1] != 'character' &  isFilename )  )
    stop(" both x and masks should be either filenames or image lists")

  # build a training dataset by sampling each image and binding together
  # the output of the matrices
  yvec = matrix( )
  xmat = matrix( )
  for ( i in 1:length( x ) ) {
    if ( isFilename ) yimg = ANTsRCore::antsImageRead( y[ i ] ) else yimg = y[[ i ]]
    if ( isFilename ) img = ANTsRCore::antsImageRead( x[ i ] ) else img = x[[ i ]]
    if ( isFilename ) msk = ANTsRCore::antsImageRead( masks[ i ] ) else msk = masks[[ i ]]
    ripped = ripmmarc( img, msk, patchRadius = patchRadius, meanCenter = FALSE,
        patchSamples = nsamples,
        evecBasis = rcb$basisMat, patchVarEx = nrow(rcb$basisMat),
        canonicalFrame = rcb$canonicalFrame, regressProjections = TRUE,
        rotationInvariant = FALSE )
    if ( i == 1 ) {
      xmat = ripped$evecCoeffs
      yvec = matrix( yimg[  msk == 1 ], ncol = 1 )
      spatialMat = imageDomainToSpatialMatrix( msk, msk )
      } else {
      xmat = rbind( xmat, ripped$evecCoeffs )
      yvec = rbind( yvec, matrix( yimg[  msk == 1 ], ncol = 1 ) )
      spatialMat = rbind( spatialMat, imageDomainToSpatialMatrix( msk, msk ) )
      }
    }

  return( list( x = xmat, y = yvec, position = spatialMat ) )

}




#' rcTestingMatrix
#'
#' Build a rcissus testing matrix
#'
#' @param x input filenames or image list defining one testing modality
#' @param masks input filenames or image list defining testing modality masks
#' @param rcb input image basis from \code{\link{rcBasis}}
#' @param patchRadius patch radius, integer value
#' @param nsamples number of samples per image
#' @param seeds random seeds for reproducibility across testing modalities
#' @return testing matrix is output in list
#' @author Avants BB
#'
#' @export rcTestingMatrix
rcTestingMatrix <- function( x, masks, rcb, patchRadius = 3, nsamples = 1000, seeds ) {

  if ( length( x ) < 1 ) stop( "pass in a non-zero length input" )

  if ( class( x[[1]] )[1] == 'character' )
    isFilename = TRUE else isFilename = FALSE

  if ( missing( seeds ) )
    seeds = sample( 1:10000 )[1:length(x)]

  if ( length( x ) != length( seeds ) ) stop( "length of seeds should equal length of x" )

  if (
    ( class( masks[[1]] )[1] == 'character' & !isFilename ) |
    ( class( masks[[1]] )[1] != 'character' &  isFilename )  )
    stop(" both x and masks should be either filenames or image lists")

  # build a training dataset by sampling each image and binding together
  # the output of the matrices
  xmat = matrix( )
  for ( i in 1:length( x ) ) {
    if ( isFilename ) img = ANTsRCore::antsImageRead( x[ i ] ) else img = x[[ i ]]
    if ( isFilename ) msk = ANTsRCore::antsImageRead( masks[ i ] ) else msk = masks[[ i ]]
    ripped = ripmmarc( img, msk, patchRadius = patchRadius, meanCenter = FALSE,
        patchSamples = nsamples,
        evecBasis = rcb$basisMat, patchVarEx = nrow(rcb$basisMat),
        canonicalFrame = rcb$canonicalFrame, regressProjections = TRUE,
        rotationInvariant = FALSE )
    if ( i == 1 ) {
      xmat = ripped$evecCoeffs
      spatialMat = imageDomainToSpatialMatrix( msk, msk )
      } else {
      xmat = rbind( xmat, ripped$evecCoeffs )
      spatialMat = rbind( spatialMat, imageDomainToSpatialMatrix( msk, msk ) )
      }
    }

  return( list( x = xmat, position = spatialMat ) )

}





#' rcTrain
#'
#' Train a rcissus model
#'
#' @param y outcome vector
#' @param trainingDf input training data
#' @return model is output
#' @author Avants BB
#' @importFrom stats lm
#' @seealso \code{\link[ANTsRCore]{ripmmarcPop}} \url{https://antsx.github.io/ANTsRCore/reference/ripmmarcPop.html}
#'
#' @export rcTrain
rcTrain <- function( y, trainingDf ) {
  trainingDf$y = y
  if ( !usePkg("h2o") ) {
    mdl = lm( y ~ . , data = trainingDf )
    return( mdl )
  } else {
    localH2O <- h2o::h2o.init( nthreads = -1 , max_mem_size = "100G" ) # all cores
    h2o::h2o.init()
    tempath = tempfile( pattern = "h2ofile", tmpdir = tempdir(), fileext = ".csv")
    write.csv( trainingDf, tempath, row.names = FALSE )
    train.h2o <- h2o::h2o.importFile( tempath )
    mdl = h2o::h2o.deeplearning( y = "y",  training_frame = train.h2o, epochs=200 )
    return( mdl )
  }
}




#' rcPredict
#'
#' Predict from a rcissus model
#'
#' @param mdl input trained model
#' @param testingDf input testing data
#' @return prediction is output
#' @author Avants BB
#' @importFrom stats lm
#'
#' @export rcPredict
rcPredict <- function( mdl, testingDf ) {
  if ( !usePkg("h2o") ) {
    return( predict( mdl ) )
  } else {
    tempath = tempfile( pattern = "h2otestfile", tmpdir = tempdir(), fileext = ".csv")
    write.csv( testingDf, tempath, row.names = FALSE )
    h2otest <- h2o.importFile( tempath )
    return( as.data.frame( h2o.predict( mdl, h2otest ) )[,1] )
  }
}
