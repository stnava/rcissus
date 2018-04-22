# rcissus

Rcissus provides a fast approach to patch-based deep learning segmentation or
image translation (mapping intensities between imaging modalities).

rcissus uses an eigenvector representation of image patches in order to allow
rapid and generalizable training from small $n$ datasets.  the representation is
based on [RIPMMARC](https://www.ncbi.nlm.nih.gov/pubmed/25449745) and [RIPMMARC-POP](https://ww5.aievolution.com/hbm1701/index.cfm?do=abs.viewAbs&abs=2175).

the model training instances are determined by the number of patches, not the
number of image examples, making the approach relevant to much smaller datasets
potentially reducing the need for augmentation and high-performance GPUs. these
models run fairly efficiently on CPU architecture.

the keras-based models are much more flexible than those provided by h2o and
will likely serve as the foundation for further development.

```
devtools::install_github("stnava/rcissus")
```

an example predicting raw intensity from gradient and laplacian images.  The
`rcissus` vignette may be a better place to start than here.

```
mth = 'kerasdnn' # may require more tuning than h2o but can yield better results
mth = 'h2o'
nBases = 5
# step 1 - collect data
library( rcissus )
popfns = getANTsRData('show')[1:5]
testfn = getANTsRData('show')[6]
# below, we use image lists but one can alternatively replaces lists
# with vectors of filenames - just do so consistently
popGT = list( )  # ground truth
popGR = list( )  # gradient
popLA = list( )  # laplacian
masks = list( )  # sample masks
nsam = 5000 # samples
myPR = 3
mc = FALSE
for ( i in 1:length( popfns ) ) {
  popGT[[ i ]] = antsImageRead( getANTsRData( popfns[ i ] ) )
  popGR[[ i ]] = iMath( popGT[[ i ]], "Grad", 1 )
  popLA[[ i ]] = iMath( popGT[[ i ]], "Laplacian", 1 )
  masks[[ i ]] = randomMask( thresholdImage( popGT[[ i ]], 1 , 255  ) , nsam )
  }
trnBas = rcBasis( lappend( popGR, popLA  ), patchRadius = myPR, meanCenter = mc )
trnBas$basisMat = trnBas$basisMat[ 1:nBases,  ] # select basis vectors
myseeds = c( 1:length( popGT ) )
trnMat1 = rcTrainingMatrix( popGT, popGR, masks, trnBas, seeds = myseeds, patchRadius = myPR, meanCenter = mc   )
trnMat2 = rcTrainingMatrix( popGT, popLA, masks, trnBas, seeds = myseeds, patchRadius = myPR, meanCenter = mc   )
print( table( trnMat2$y==trnMat1$y ) )

popGTtest = list( )  # ground truth
popGRtest = list( )  # gradient
popLAtest = list( )  # laplacian
maskstest = list( )  # sample masks
for ( i in 1:length( testfn ) ) {
  popGTtest[[ i ]] = antsImageRead( getANTsRData( testfn[ i ] ) )
  popGRtest[[ i ]] = iMath( popGTtest[[ i ]], "Grad", 1 )
  popLAtest[[ i ]] = iMath( popGTtest[[ i ]], "Laplacian", 1 )
  maskstest[[ i ]] = getMask( popGTtest[[ i ]] ) # NOTE: dense prediction!
  }
testMat1 = rcTestingMatrix( popGRtest, maskstest, trnBas, seeds = 1, patchRadius = myPR, meanCenter = mc  )
testMat2 = rcTestingMatrix( popLAtest, maskstest, trnBas, seeds = 1, patchRadius = myPR, meanCenter = mc  )


traindf = data.frame( trnMat1$x, trnMat2$x, trnMat1$position ) / 1000
testdf = data.frame( testMat1$x, testMat2$x, testMat1$position ) / 1000
# if h2o works on your machine, use deep learning
trn = rcTrain( trnMat1$y, traindf, mdlMethod = mth, epochs = 100 )
prd = rcPredict( trn, testdf, mdlMethod = mth )
mm = makeImage( maskstest[[1]], as.numeric( prd[,1] ) )
plot( mm )

```

The same approach can be used for segmentation.

```
# step 1 - collect data
# below, we use image lists but one can alternatively replaces lists
# with vectors of filenames - just do so consistently
popGT = list( )  # ground truth
popGR = list( )  # gradient
masks = list( )  # sample masks
nsam = 5000 # samples
myPR = 4
for ( i in 1:length( popfns ) ) {
  popGT[[ i ]] = antsImageRead( getANTsRData( popfns[ i ] ) )
  popGR[[ i ]] = thresholdImage( popGT[[ i ]], "Otsu", 3 )
  masks[[ i ]] = randomMask( thresholdImage( popGT[[ i ]], 1, 255  ) , nsam )
  }
trnBas = rcBasis( popGT, patchRadius = myPR )
trnBas$basisMat = trnBas$basisMat[ 1:12,  ] # select 15 basis vectors
myseeds = c( 1:length( popGT ) )
trnMat1 = rcTrainingMatrix( popGR, popGT, masks, trnBas, seeds = myseeds, patchRadius = myPR  )

# step 2 - build similar data for prediction - but use a dense mask
popGTtest = list( )  # ground truth
popGRtest = list( )  # gradient
maskstest = list( )  # sample masks
for ( i in 1:length( testfn ) ) {
  popGTtest[[ i ]] = antsImageRead( getANTsRData( testfn[ i ] ) )
  popGRtest[[ i ]] = thresholdImage( popGTtest[[ i ]], "Otsu", 3 )
  maskstest[[ i ]] = getMask( popGTtest[[ i ]] ) # NOTE: dense prediction!
  }
testMat1 = rcTestingMatrix( popGTtest, maskstest, trnBas, seeds = 1, patchRadius = myPR )

# step 3 - now implement the training and testing
# - note that we scale the data by trivial means
traindf = ( data.frame( trnMat1$x, position=trnMat1$position ) ) / 10000
testdf = ( data.frame( testMat1$x, position=testMat1$position ) ) / 10000
# if h2o works on your machine, use deep learning
trn = rcTrain( trnMat1$y, traindf, classification = TRUE, mdlMethod=mth, epochs=33 )
prd = rcPredict( trn, testdf, classification = TRUE, mdlMethod=mth )
prdnumerics = as.numeric(prd$predict)  - min( as.numeric(prd$predict) )
mm = makeImage( maskstest[[1]], prdnumerics  )
plot( popGTtest[[ 1 ]], mm, window.overay = range( mm ), alpha = 0.6 )
sel = maskstest[[ 1 ]] == 1
gtvals = popGRtest[[ 1 ]][ sel ]
print( table( gtvals == prdnumerics  )  )
############ the ground truth ############
plot( popGTtest[[ 1 ]], popGRtest[[ 1 ]], window.overay = range( mm ), alpha = 0.6 )
mm = makeImage( maskstest[[1]], prd[ , "class_2" ]  )
plot( mm  )
```


Note: one may want to customize the deep learning components used here, which
are simply the defaults provided in h2o.  One can inspect the code in order to
get started on this, as well as the h2o documentation.  In order to optimize
performance for a specific problem, one will need to define a validation objective
and take advantage of the many parameter search strategies available.
