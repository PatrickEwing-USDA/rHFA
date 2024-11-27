## Home field advantage (HFA) - Analysis of local adaptation in crop populations.

Home Field Advantage (HFA) is the fitness gain an entry (genotype) realizes by growing in its home location, which is empirically defined as the location of highest fitness relative to other entries (Ewing *et al*., 2019).

### Requirements

Input data should be analogous to common gardens. Ex. breeding trials, multi-environment trials, yield trials.

### Brief how-to:

1.  `permute_hfa()` to run hfa analysis. This will work at the population, genotype, site, and year level.
2.  `id_home()` to return the home site of each genotype
3.  `specialist_test()` and `generalist_test()` to identify specialists and generalists.

The rHFA Guide vignette contains a step-by-step analysis of sample data. When installing the package from github, use this command to ensure the vignette is included.

`> remotes::install_github("PatrickEwing-USDA/rHFA", build_vignettes = TRUE)`

Then:

`> vignette('rhfa-guide', package='rHFA')`

See analysis code for publications below and methods described in those publications - especially Ewing *et al*. (2019).

### Publications:

-   Ewing, P. M., Runck, B. C., Kono, T. Y., & Kantar, M. B. (2019). The home field advantage of modern plant breeding. *PloS one 14*(12), e0227079. [Link](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0227079)
-   MacQueen, Alice H., et al. (2022). Local to continental‐scale variation in fitness and heritability in common bean. *Crop Science 62*(2): 767-779. [Link](https://acsess.onlinelibrary.wiley.com/doi/full/10.1002/csc2.20694) [Code](https://github.com/Alice-MacQueen/cdbn-home-away)
-   Ewing, P. M., et al. (2024). Local adaptation and broad performance are synergistic to productivity in modern barley. *Crop Science, 64*(1), 192-199. [Link](https://acsess.onlinelibrary.wiley.com/doi/full/10.1002/csc2.21168) [Code](https://doi.org/10.5281/zenodo.10267964)

*Disclaimer: Mention of trade names or commercial products in this publication is solely for the purpose of providing specific information and does not imply recommendation or endorsement by the U.S. Department of Agriculture. The USDA is an equal opportunity provider and employer.*

## Issues and Development ##
rHFA is still under development! The following are known issues:

1. When `method="blup"`, permutation testing with `permute_hfa()` may stop on error or the significance test may produce an NA. These are related issues. Either increase `min_times` in `*filter_instances()` or try a different method, such as `method='median'`.
2. If using `remotes::install_github()`, you may need to run `install.packages('ps')` and restart R.

We are also actively developing:

1. Friendly `plot()` methods to easily recreate graphs as in the above-cited publications.
2. Methods to calculate additional breeding metrics.

Find another bug? Have ideas for enhancements? General questions? Raise an [issue](https://github.com/PatrickEwing-USDA/rHFA/issues). 
