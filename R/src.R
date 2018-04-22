#' rcBasis
#'
#' Learn a rcissus image basis
#'
#' @param x input filenames or image list
#' @param patchRadius patch radius, integer value
#' @param meanCenter boolean
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
#' @importFrom ANTsRCore usePkg
#' @importFrom ANTsRCore randomMask
#' @importFrom stats predict
#' @importFrom magrittr %>%
#' @import methods
#' @import h2o
rcBasis <- function( x, patchRadius = 3, meanCenter = FALSE  ) {
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
#' @param meanCenter boolean
#' @param nsamples number of samples per image
#' @param seeds random seeds for reproducibility across training modalities
#' @return training matrix and vector (ground truth) is output in list
#' @author Avants BB
#' @seealso \code{\link[ANTsRCore]{ripmmarcPop}} \url{https://antsx.github.io/ANTsRCore/reference/ripmmarcPop.html}
#'
#' @export rcTrainingMatrix
rcTrainingMatrix <- function( y, x, masks, rcb, patchRadius = 3,
  meanCenter = FALSE, nsamples = 1000, seeds ) {
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
#' @param meanCenter boolean
#' @param nsamples number of samples per image
#' @param seeds random seeds for reproducibility across testing modalities
#' @return testing matrix is output in list
#' @author Avants BB
#'
#' @export rcTestingMatrix
rcTestingMatrix <- function( x, masks, rcb, patchRadius = 3, meanCenter = FALSE,
  nsamples = 1000, seeds ) {

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
#' Train a rcissus model with \code{lm} or deep learning from \code{h2o}
#'
#' @param y outcome vector
#' @param trainingDf input training data
#' @param hidden Hidden layer sizes (e.g. \code{c(100,100)}). Defaults to \code{c(200,200)}.
#' @param classification boolean
#' @param mdlMethod string either lm, h2o or kerasdnn
#' @param epochs number of epochs (integer) over which to train
#' @param max_mem sets maximum allowable memory for h2o deep learning
#' @return model is output
#' @author Avants BB
#' @importFrom stats lm
#' @seealso \code{\link[ANTsRCore]{ripmmarcPop}} \url{https://antsx.github.io/ANTsRCore/reference/ripmmarcPop.html}
#'
#' @export rcTrain
rcTrain <- function( y,
  trainingDf,
  hidden = c( 200, 200 ),
  classification = FALSE,
  mdlMethod = 'h2o',
#  nfolds = 5,
  epochs = 200,
  max_mem = "100G" ) {
  if ( mdlMethod == 'lm' ) {
    trainingDf$y = y
    if ( classification ) trainingDf$y = factor( paste0( "class_", as.character( y ) ) )
    mdl = lm( y ~ . , data = trainingDf )
    return( mdl )
  }
  if ( mdlMethod == 'h2o' ) {
    trainingDf$y = y
    if ( classification ) trainingDf$y = factor( paste0( "class_", as.character( y ) ) )
    localH2O <- h2o::h2o.init( nthreads = -1 , max_mem_size = max_mem ) # all cores
    h2o::h2o.init()
    tempath = tempfile( pattern = "h2ofile", tmpdir = tempdir(), fileext = ".csv")
    write.csv( trainingDf, tempath, row.names = FALSE )
    train.h2o <- h2o::h2o.importFile( tempath )
    mdl = h2o::h2o.deeplearning( y = "y",
      training_frame = train.h2o,
#      activation = "Rectifier",    ## default
      hidden = hidden,       ## default: 2 hidden layers with 200 neurons each
      epochs = epochs )
    return( mdl )
  }
  if ( mdlMethod == 'kerasdnn' & ANTsRCore::usePkg( "keras" ) ) {
    mod <- keras::keras_model_sequential() %>%
      keras::layer_dense( units = hidden[1], activation = 'relu',
        input_shape = ncol( trainingDf ) )
    for ( k in 2:length( hidden ) ) {
      mod <- keras::layer_dense( mod, units = hidden[k], activation = 'relu',
        input_shape = hidden[k-1] )
      }
    if ( classification ) {
      mod <- keras::layer_dense( mod, units = ncol( y ), activation = 'softmax' ) %>%
       keras::compile( loss = 'categorical_crossentropy',
         optimizer = keras::optimizer_rmsprop(), metrics = c('accuracy') )
      }
    if ( ! classification ) {
      mod <- keras::layer_dense( mod, units = ncol( y )  ) %>%
      keras::compile( loss = 'mean_squared_error',
        optimizer = keras::optimizer_rmsprop() )
      }
    # now compile
    btch = round( nrow( trainingDf ) / 10 )
    keras::fit( mod,
      data.matrix( trainingDf ), data.matrix( y ), # batch = btch,
        epochs = epochs, verbose = 1, validation_split = 0.1 )
    return( mod )
  }
  return( NA )
}




#' rcPredict
#'
#' Predict from a rcissus model
#'
#' @param mdl input trained model
#' @param testingDf input testing data
#' @param classification boolean
#' @param mdlMethod string either lm, h2o or kerasdnn
#' @return prediction is output
#' @author Avants BB
#' @importFrom stats lm
#'
#' @export rcPredict
rcPredict <- function( mdl,
  testingDf,
  classification = FALSE,
  mdlMethod = 'h2o' ) {
  if ( mdlMethod == 'lm' ) {
    return( predict( mdl, newdata = testingDf ) )
    }
  if ( mdlMethod == 'kerasdnn' ) {
    if ( classification ) {
      outdf = keras::predict_classes( mdl, x = data.matrix( testingDf ) )
      print( dim( outdf ) )
      return( outdf )
    }
    if ( !classification ) {
      outdf = predict( mdl, x = data.matrix( testingDf ) )
      return( outdf )
    }
  }
  if ( mdlMethod == 'h2o' ) {
    tempath = tempfile( pattern = "h2otestfile", tmpdir = tempdir(), fileext = ".csv")
    write.csv( testingDf, tempath, row.names = FALSE )
    h2otest <- h2o.importFile( tempath )
    if ( classification ) {
      outdf = as.data.frame( h2o.predict( mdl, h2otest ) )
    } else {
      outdf = as.data.frame( h2o.predict( mdl, h2otest )[,1] )
    }
    return( outdf )
  }
  return( NA )
}
