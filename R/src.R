#' rcBasis
#'
#' Learn a rcissus image basis
#'
#' @param x input filenames or image list
#' @param patchRadius patch radius, integer value
#' @param meanCenter boolean
#' @param nsamples number of samples
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
#' @importFrom ANTsR writeNormalizedPopulationData
#' @importFrom ANTsR readNormalizedPopulationData
#' @importFrom stats predict
#' @importFrom knitr knit
#' @importFrom magrittr %>%
#' @import methods
#' @import h2o
rcBasis <- function( x, patchRadius = 3, meanCenter = FALSE, nsamples = 1000  ) {
  maskType = 'sparse'
  if ( ! ( maskType %in% c("sparse","dense") ) ) stop("pass dense or sparse as maskType" )
  if ( length( x ) < 1 ) stop( "pass in a non-zero length input" )

  if ( class( x[[1]] )[1] == 'character' )
    isFilename = TRUE else isFilename = FALSE

  # build a training dataset by sampling each image and binding together
  # the output of the matrices
  popmasks = list()
  for ( i in 1:length( x ) ) {
    if ( isFilename ) img = ANTsRCore::antsImageRead( x[ i ] ) else img = x[[ i ]]
    if ( maskType == 'dense' )
      popmasks[[ i ]] = ANTsRCore::getMask( img, cleanup = 0  )
    if ( maskType == 'sparse' ) {
      temp = ANTsRCore::getMask( img, cleanup = 0 )
      popmasks[[ i ]] = ANTsRCore::randomMask( temp, nsamples, perLabel = TRUE ) * temp
      }
    }
  if ( isFilename ) {
    rp = ANTsRCore::ripmmarcPop( ANTsRCore::imageFileNames2ImageList( x ),
      popmasks, patchRadius=patchRadius, meanCenter = FALSE, patchSamples=nsamples )
    } else {
    rp = ANTsRCore::ripmmarcPop( x,
      popmasks, patchRadius=patchRadius, meanCenter = FALSE, patchSamples=nsamples )
    }
  return( rp )
}




#' rcTrainingMatrix
#'
#' Build a rcissus training matrix
#'
#' @param x input filenames or image list defining one training modality
#' @param y input filenames or image list defining ground truth
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
rcTrainingMatrix <- function( x, y, masks, rcb, patchRadius = 3,
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
  if ( mdlMethod == 'kerasdnn'  ) {
    myact = 'relu'
    dropRate = 0.01
    mod <- keras::keras_model_sequential() %>%
      keras::layer_dense(
        units = hidden[1],
        activation = myact,
        input_shape = ncol( trainingDf ),
        kernel_initializer = keras::initializer_random_normal() ) %>%
#        activity_regularizer = keras::regularizer_l1_l2(l1 = 0.01, l2 = 0.01) )  %>%
        keras::layer_dropout(rate = dropRate ) #  %>% keras::layer_batch_normalization()
    for ( k in 2:length( hidden ) ) {
      mod <- keras::layer_dense( mod,
        units = hidden[k],
        activation = myact,
        input_shape = hidden[k-1] ) %>%
#        activity_regularizer = keras::regularizer_l1_l2(l1 = 0.01, l2 = 0.01) ) %>%
        keras::layer_dropout(rate = dropRate )  #  %>%  keras::layer_batch_normalization()
      }
    if ( classification ) {
      losswmx = max( table( y ) )
      losswtbl = ( as.numeric( table(y )) )
      lossw = ( losswtbl / losswmx )^(-1)
      lossw = list( lossw / sum( lossw ) )
      y = keras::to_categorical( y )
      mod <- keras::layer_dense( mod, units = ncol( y ), activation = 'softmax' )
      keras::compile( mod, loss = 'categorical_crossentropy', # loss_weights = lossw,
         optimizer = keras::optimizer_adam() )
      }
    if ( ! classification ) {
      mod <- keras::layer_dense( mod, units = ncol( y )  ) %>%
      keras::compile( loss = 'mean_squared_error',
        optimizer = keras::optimizer_adam() )
      }
    # now compile
    btch = round( nrow( trainingDf ) / 10 )
    keras::fit( mod,
      data.matrix( trainingDf ), data.matrix( y ), # batch = btch,
        epochs = epochs, verbose = 1, validation_split = 0.2 )
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
  mdlMethod ) {
  if ( missing( mdlMethod ) ) {
    if ( class( mdl )[1] == 'H2ORegressionModel'  ) {
      mdlMethod = 'h2o'
      classification = FALSE
    }
    if ( class( mdl )[1] == 'H2OMultinomialModel'  ) {
      mdlMethod = 'h2o'
      classification = TRUE
    }
    if ( class( mdl )[1] == 'keras.models.Sequential'  ) {
      mdlMethod = 'kerasdnn'
    }
    if ( class( mdl )[1] == 'lm'  ) {
      mdlMethod = 'lm'
      classification = FALSE
    }
  }
  if ( mdlMethod == 'lm' ) {
    return( predict( mdl, newdata = testingDf ) )
    }
  if ( mdlMethod == 'kerasdnn' ) {
    if ( classification ) {
      mypred = keras::predict_classes( mdl, x = data.matrix( testingDf ) )
      myprobs = keras::predict_proba( mdl, x = data.matrix( testingDf ) )
      outdf = data.frame( cbind( mypred, myprobs ) )
      colnames( outdf ) = c( "predict", paste0( "class_", 0:(ncol(myprobs)-1) ) )
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
    h2otest <- h2o::h2o.importFile( tempath )
    if ( classification ) {
      outdf = as.data.frame( h2o::h2o.predict( mdl, h2otest ) )
    } else {
      outdf = as.data.frame( h2o::h2o.predict( mdl, h2otest )[,1] )
    }
    return( outdf )
  }
  return( NA )
}






#' rcTrainTranslation
#'
#' Training for rcissus translation between image modalities
#'
#' @param x input filenames or image list defining source modality
#' @param y input filenames or image list defining target modality
#' @param masks input filenames or image list defining source modality masks
#' @param nBases number of eigenvector bases
#' @param patchRadius patch radius, integer value
#' @param meanCenter boolean
#' @param nsamples number of samples per image
#' @param seeds random seeds for reproducibility across testing modalities
#' @param hidden Hidden layer sizes (e.g. \code{c(100,100)}). Defaults to \code{c(200,200)}.
#' @param classification boolean
#' @param dataScaling scaling by which to divide input data - results may be sensitive to this value
#' @param mdlMethod string either lm, h2o or kerasdnn
#' @param epochs number of epochs (integer) over which to train
#' @return list containing the image bases and the trained model and other relevant parameters.
#' @author Avants BB
#'
#' @export rcTrainTranslation
rcTrainTranslation <- function(
  x,
  y,
  masks,
  nBases = 6,
  patchRadius = 3,
  meanCenter = FALSE,
  nsamples = 1000,
  seeds,
  hidden = c( 200, 200 ),
  classification = FALSE,
  dataScaling = 1000,
  mdlMethod = 'h2o',
  epochs = 10 ) {


  ################################### critical step - build the basis set from both features
  trnBas = rcBasis( x, patchRadius = patchRadius, meanCenter = meanCenter )
  trnBas$basisMat = trnBas$basisMat[ 1:nBases,  ] # select basis vectors
  if ( missing( seeds ) )
    seeds = c( 1:length( x ) )
  ################################### project features to the basis
  trnMat1 = rcTrainingMatrix( x, y, masks, trnBas, seeds = seeds,
    patchRadius = patchRadius, meanCenter = meanCenter   )
  ################################### train/test below
  traindf = data.frame( trnMat1$x, trnMat1$position ) / dataScaling
  traindf = data.frame( trnMat1$x ) / dataScaling
  trn = rcTrain( trnMat1$y, traindf, mdlMethod = mdlMethod,
    epochs = epochs, hidden = hidden,
    classification = classification )

  # out of sample prediction ....
  dmdl = list(
      deepNet = trn,
      meanCenter = meanCenter,
      patchRadius = patchRadius,
      basis = trnBas,
      dataScaling = dataScaling )

  return( dmdl )
}



#' rcTranslate
#'
#' Rcissus translation between image modalities
#'
#' @param x input single image from source modality
#' @param rcmdl the trained model from \code{rcTrainTranslation}
#' @param mask optional input single image mask
#' @param classification boolean
#' @return the translated image.
#' @author Avants BB
#'
#' @export rcTranslate
rcTranslate <- function(
  x,
  rcmdl,
  mask,
  classification = FALSE  ) {

  if ( missing( mask ) ) {
    mask = getMask( x )
    }

  basisRepresentation = rcTestingMatrix(
    list( x ),
    list( mask ),
    rcmdl$basis,
    seeds = 1,
    patchRadius = rcmdl$patchRadius,
    meanCenter = rcmdl$meanCenter  )

  ################################### apply model below
  testdfx = data.frame( basisRepresentation$x, basisRepresentation$position ) / rcmdl$dataScaling
  testdfx = data.frame( basisRepresentation$x ) / rcmdl$dataScaling
  prd = rcPredict( rcmdl$deepNet, testdfx, classification = classification )
  output = matrixToImages( t(data.matrix(prd)), mask )
  return( output )

}










#' writeRcissus
#'
#' Write an rcissus basis set
#'
#' @param basis the basis image to write
#' @param patchRadius the radius associated with the basis
#' @param meanCenter boolean
#' @param dataScaling optional scaling value(s)
#' @param deepNet optional network
#' @param directoryname for output
#' @return succeed or fail
#' @author Avants BB
#' @examples
#'
#' \dontrun{
#' pop = getANTsRData( "population" ) # list of example images
#' imgBases = rcBasis( pop )
#' writeRcissus( imgBases, 3, FALSE, tempfile() )
#' }
#'
#' @export writeRcissus
writeRcissus <- function(
  basis,
  patchRadius,
  meanCenter,
  dataScaling = 1,
  deepNet = NA,
  directoryname ) {

dd = data.frame(
  patchRadius = rep( patchRadius, nrow( basis$basisMat ) ),
  meanCenter =  rep( meanCenter,  nrow( basis$basisMat ) ) )
mask = basis$canonicalFrame * 0 + 1
mask[ basis$canonicalFrame == 0 ] = 0
ANTsR::writeNormalizedPopulationData(
  dd,
  basis$basisMat,
  mask,
  rep( TRUE, nrow( basis$basisMat ) ),
  directoryname
  )

write.csv( dataScaling,
  paste0( directoryname, "/dataScaling.csv" ), row.names=F )
h5fn = paste0( directoryname, "/deepNet.h5" )
if ( !is.na( deepNet ) )
  keras::save_model_hdf5( deepNet, h5fn )


}



#' readRcissus
#'
#' Read an rcissus basis set
#'
#' @param directoryname to read
#' @return rcissus list of objects
#' @author Avants BB
#' @importFrom ANTsRCore makeImage
#' @importFrom utils read.csv
#' @importFrom keras load_model_hdf5
#' @importFrom keras load_model_hdf5
#' @examples
#' \dontrun{
#' pop = getANTsRData( "population" ) # list of example images
#' imgBases = rcBasis( pop )
#' tfn = tempfile()
#' writeRcissus( imgBases, 3, FALSE, directoryname = tfn )
#' imgB2 = readRcissus( tfn )
#' }
#' @export readRcissus
readRcissus <- function( directoryname ) {

temp = ANTsR::readNormalizedPopulationData( directoryname )
canframe = ANTsRCore::makeImage( temp$imageMask, temp$imageMat[1,] )
dataScaling = utils::read.csv( paste0( directoryname, "/dataScaling.csv" ) )
deepNet = NA
h5fn = paste0( directoryname, "/deepNet.h5" )
if ( file.exists( h5fn  ) )
  deepNet = keras::load_model_hdf5( h5fn )
outbasis = list(
  canonicalFrame = canframe,
  basisMat = temp$imageMat,
  mask = temp$imageMask,
  patchRadius = temp$demographics$patchRadius[1],
  meanCenter = temp$demographics$meanCenter[1],
  deepNet = deepNet,
  dataScaling = dataScaling
  )
return( outbasis )
}
