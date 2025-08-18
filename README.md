
<!-- README.md is generated from README.Rmd. Please edit that file -->

# liaendmembers

<!-- badges: start -->

<!-- badges: end -->

This packege provides a convenient way to identify endmembers for the
purpouse of provenance analysis, using Principle component analys
(Albarede et al. 2024), and to group the components according to their
distribution along the slope of the geochron (Eshel et al. 2019).

## Installation

You can install the development version of liaendmembers like so:

``` r
# devtools::install_github()
```

## Datasets

This package comes with Lead Isotope analyse of Iron Age silver hoards
from Tel Dor, Akko, Arad, Ein Hofez and Eshtemoa in Israel (Eshel et al.
2019).

``` r
head(liaendmembers::tel_dor)
#>     pb64   pb74   pb84
#> 1 18.987 15.702 39.043
#> 2 18.971 15.694 39.023
#> 3 18.981 15.703 39.053
#> 4 18.988 15.695 39.029
#> 5 18.976 15.695 39.028
#> 6 17.910 15.645 38.003
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(liaendmembers)
## Using the endmemebrs function
data("tel_dor")
end_members <- endmembers(
     tel_dor,
     colnames(tel_dor),
     tolerance = c(0.01, 0.01),
     clamp = c(Inf, Inf)
)
#> PC2 or PC3 are not normally distributed. This may indicate that their variation may not be random noise.

## Print summary of the liaendmembers object
summary(end_members)
#> Summary of End memebers:
#> 
#> Tolarance: 0.01 0.01 
#> Clamp: Inf Inf 
#> 
#> PCA Endmembers
#>      pb64   pb74   pb84
#> 9  17.906 15.643 37.991
#> 10 19.034 15.732 39.382
#> 
#>        Group1 Group2 Mixing Total
#> Counts      3     12     17    32
#> - - - - - - - - - - - - - - - - - - 
#> Importance of components:
#>                           PC1     PC2     PC3
#> Standard deviation     0.4957 0.05492 0.01592
#> Proportion of Variance 0.9869 0.01211 0.00102
#> Cumulative Proportion  0.9869 0.99898 1.00000
```

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-A.D.G+24" class="csl-entry">

Albarede, Francis, Gillan Davis, Liesel Gentelli, Janne Blichert-Toft,
Haim Gitler, Marine Pinto, and Philippe Telouk. 2024. “Bullion Mixtures
in Silver Coinage from Ancient Greece and Egypt.” *Journal of
Archaeological Science* 162 (February): 105918.
<https://doi.org/10.1016/j.jas.2023.105918>.

</div>

<div id="ref-E.E.Y+19" class="csl-entry">

Eshel, T\>, Y. Erel, N. Yahalom-Mack, O. Tirosh, and A. Gilboa. 2019.
“Lead Isotopes in Silver Reveal Earliest Phoenician Quest for Metals in
the West Mediterranean.” *Proceedings of the National Academy of
Sciences* 116 (13): 6007–12. <https://doi.org/10.1073/pnas.1817951116>.

</div>

</div>
