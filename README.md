[![Build Status](https://travis-ci.org/stnava/rcissus.png?branch=master)](https://travis-ci.org/stnava/rcissus)

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

install rcissus via `devtools` in `R`:

```
devtools::install_github("stnava/rcissus")
```

primarily, there are two high level functions: `rcTrainTranslation` and
`rcTranslate`.  see the [vignette](https://github.com/stnava/rcissus/blob/master/vignettes/rcissus.Rmd)
for more details: [here](https://htmlpreview.github.io/?https://github.com/stnava/rcissus/blob/master/inst/doc/rcissus.html).

## notes on parameters

* more complex / detailed modeling may require more hidden nodes / layers -- *not necessarily more basis vectors* ... a 2-layer model with 1200 nodes produced good results in 3D. not a coincidence - this is roughly the size of the information that we want to encode.  a rule of thumb could be to divide the number of voxels you want to predict into the number of layers.   nodes per layer would equal sqrt n-voxels for 2 layers, cube root for 3 layers, etc.

* it is reasonable to pilot models in 2D and then extend to 3D - e.g. with slices of a translation problem

* h2o is easier to start with, in my experience, but keras can yield better results

## to do list

* read / write basis vectors ( library of basis sizes, given radius and problem ).  this will mean we do not need to recompute bases.

* read / write models to allow reuse - not sure possible with h2o. see `?save_model_hdf5`
