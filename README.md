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
