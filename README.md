# rcissus

fast approach to mapping between image representations or modalities

```
devtools::install_github("stnava/rcissus")
```

an example predicting raw intensity from gradient and laplacian images

```
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
for ( i in 1:length( popfns ) ) {
  popGT[[ i ]] = antsImageRead( getANTsRData( popfns[ i ] ) )
  popGR[[ i ]] = iMath( popGT[[ i ]], "Grad", 1 )
  popLA[[ i ]] = iMath( popGT[[ i ]], "Laplacian", 1 )
  masks[[ i ]] = randomMask( thresholdImage( popGT[[ i ]], 1 , 255  ) , nsam )
  }
trnBas = rcBasis( lappend( popGR, popLA  ), patchRadius = myPR )
trnBas$basisMat = trnBas$basisMat[ 1:8,  ] # select 15 basis vectors
myseeds = c( 1:length( popGT ) )
trnMat1 = rcTrainingMatrix( popGT, popGR, masks, trnBas, seeds = myseeds, patchRadius = myPR  )
trnMat2 = rcTrainingMatrix( popGT, popLA, masks, trnBas, seeds = myseeds, patchRadius = myPR  )
print( table( trnMat2$y==trnMat1$y ) ) # should be equivalent


# step 2 - build similar data for prediction - but use a dense mask
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
testMat1 = rcTestingMatrix( popGRtest, maskstest, trnBas, seeds = 1, patchRadius = myPR )
testMat2 = rcTestingMatrix( popLAtest, maskstest, trnBas, seeds = 1, patchRadius = myPR )

# step 3 - now implement the training and testing
traindf = data.frame( trnMat1$x, trnMat2$x, trnMat1$position )
testdf = data.frame( testMat1$x, testMat2$x, testMat1$position )
# if h2o works on your machine, use deep learning
trn = rcTrain( trnMat1$y, traindf )
prd = rcPredict( trn, testdf )
mm = makeImage( maskstest[[1]], prd  )
plot( mm )


# if h2o does not work, you can use other models
library( randomForest )
library( e1071 )  # alternatively could use h2o or tensorflow
# model interactions with position
# trn = svm(  trnMat1$y ~  ( . ) * X1 * X2 , data = traindf )
# prd = predict( trn, newdata = testdf )
# mm = makeImage( maskstest[[1]], prd  )
# plot( mm )



```



The same approach can be used for segmentation.



```


# step 1 - collect data
library( rcissus )
popfns = getANTsRData('show')[1:5]
testfn = getANTsRData('show')[6]
# below, we use image lists but one can alternatively replaces lists
# with vectors of filenames - just do so consistently
popGT = list( )  # ground truth
popGR = list( )  # gradient
masks = list( )  # sample masks
nsam = 5000 # samples
myPR = 3
for ( i in 1:length( popfns ) ) {
  popGT[[ i ]] = antsImageRead( getANTsRData( popfns[ i ] ) )
  popGR[[ i ]] = thresholdImage( popGT[[ i ]], "Otsu", 3 )
  masks[[ i ]] = randomMask( thresholdImage( popGT[[ i ]], 1, 255  ) , nsam )
  }
trnBas = rcBasis( popGT, patchRadius = myPR )
trnBas$basisMat = trnBas$basisMat[ 1:15,  ] # select 15 basis vectors
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
traindf = data.frame( trnMat1$x, position=trnMat1$position )
testdf = data.frame( testMat1$x, position=testMat1$position )
# if h2o works on your machine, use deep learning
trn = rcTrain( trnMat1$y, traindf, classification = TRUE )
prd = rcPredict( trn, testdf )
mm = makeImage( maskstest[[1]], prd  )
plot( mm )
```


Note: one may want to customize the deep learning components used here, which
are simply the defaults provided in h2o.  One can inspect the code in order to
get started on this, as well as the h2o documentation.
