# Home field advantage (HFA) analysis of local adaptation in crop populations.

Home Field Advantage (HFA) is the fitness gain an entry (genotype) realizes by growing in its home location, which is empirically defined as the location of highest fitness relative to other entries (Ewing et al., 2019).

## Requirements

Input data should be analogous to common gardens. Ex. breeding trials, multi environment trials, yield trials. This package was tested on Linux, Windows, and macOS running R version 4.4.x

## Brief how-to:

## Getting Started

To learn how to use `mypackage`, check out our [introduction vignette](vignettes/rhfa-guide.Rmd), which walks you through the core functionality.

1.  id_home() to identify home site (optional)
2.  permute_hfa() to run hfa analysis. This will work at the population and genotype level and may work across years (for populations) - in progress. You currently need to manually calculate p-values for genotypes from the permutations (first column is observed; try .two-tailed()). If really slow, try setting BLUP=FALSE.
3.  temporal_hfa() for calculationg hfa across years at the population level.

See analysis code for publications below and methods described in those publications - especially Ewing et al. (2019). All functions are relatively documented in source code.

The median (quantile) regression approach of Ewing et al. (2019) has not yet been transferred to this package - only ordinary least squares (`lm()`-based) and BLUP (`lmer()`-based) are. Of the two, BLUP is more robust to low sample sizes (site-years per genotype) and outliers but also is substantially slower.

Of course, message with questions.

Publications:

-   Ewing, P. M., Runck, B. C., Kono, T. Y., & Kantar, M. B. (2019). The home field advantage of modern plant breeding. PloS one, 14(12), e0227079. [Link](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0227079)
-   MacQueen, Alice H., et al. (2022). "Local to continental‐scale variation in fitness and heritability in common bean." *Crop Science 62*(2): 767-779. [Link](https://acsess.onlinelibrary.wiley.com/doi/full/10.1002/csc2.20694) [Code](https://github.com/Alice-MacQueen/cdbn-home-away)
-   Ewing, P. M., et al. (2024). Local adaptation and broad performance are synergistic to productivity in modern barley. *Crop Science, 64*(1), 192-199. [Link](https://acsess.onlinelibrary.wiley.com/doi/full/10.1002/csc2.21168) [Code](https://doi.org/10.5281/zenodo.10267964)

*Disclaimer: Mention of trade names or commercial products in this publication is solely for the purpose of providing specific information and does not imply recommendation or endorsement by the U.S. Department of Agriculture. The USDA is an equal opportunity provider and employer.*
