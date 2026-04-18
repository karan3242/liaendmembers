
<!-- README.md is generated from README.Rmd. Please edit that file -->

# liaendmembers

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16909256.svg)](https://doi.org/10.5281/zenodo.16909256)
<!-- badges: end -->

This package provides convenient tools for provenance analysis using
lead isotope data. Specifically, the software enables researchers to:

- Identify endmembers for provenance analysis using Principal Component
  Analysis (Albarede, Davis, Gentelli, et al. 2024) and group components
  by their distribution along the geochron slope (Eshel et al. 2019).

- Calculate nearest distances between sample and reference data within
  ore isotope databases using Euclidean metrics with or without
  mass-dependent fractionation corrections (Albarede, Davis,
  Blichert-Toft, et al. 2024).

- Determine probabilities of samples originating from specific
  subregional groups as defined by established isotope databases
  (**S.K.D+26?**).

## Installation

You can install the development version of liaendmembers like so:

    devtools::install_github("karan3242/liaendmembers")

## Datasets

This package comes with Lead Isotope analyse of Iron Age silver hoards
from Tel Dor, Akko, Arad, Ein Hofez and Eshtemoa in Israel (Eshel et al.
2019),

``` r
head(liaendmembers::tel_dor)
#>     pb64   pb74   pb84
#> 1 18.987 15.702 39.043
#> 2 18.971 15.694 39.023
#> 3 18.981 15.703 39.053
#> 4 18.988 15.695 39.029
#> 5 18.976 15.695 39.028
#> 6 17.910 15.645 38.003
head(liaendmembers::ref_data)
#>    groups     pb64     pb74     pb84
#> 1 Bezarak 59.66587 18.58592 41.85561
#> 2 Bezarak 18.82176 15.62394 39.00056
#> 3 Bezarak 44.84305 18.07175 63.92825
#> 4 Bezarak 24.98751 15.70215 46.43428
#> 5 Bezarak 18.85014 15.71725 39.09331
#> 6 Bezarak 18.53225 15.65048 38.76575
```

along with a pretrained preovence `xgb.Booster` model for prevenance
prediction.

``` r
head(liaendmembers::model, 3)
#> $Bezarak
#> ##### xgb.Booster
#> call:
#>   xgboost::xgb.train(params = params, data = dtrain, nrounds = .nrounds)
#> # of features: 3 
#> # of rounds:  100 
#> 
#> $Kabul
#> ##### xgb.Booster
#> call:
#>   xgboost::xgb.train(params = params, data = dtrain, nrounds = .nrounds)
#> # of features: 3 
#> # of rounds:  100 
#> 
#> $`Touissit Bou Beker`
#> ##### xgb.Booster
#> call:
#>   xgboost::xgb.train(params = params, data = dtrain, nrounds = .nrounds)
#> # of features: 3 
#> # of rounds:  100
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(liaendmembers)
data("tel_dor")
end_members <- endmembers(
  tel_dor,
  colnames(tel_dor),
  tolerance = c(0.01, 0.01),
  clamp = c(Inf, Inf)
)
#> PC2 or PC3 are not normally distributed. This may indicate that their variation may not be random noise.
# Print summary of the liaendmembers object
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
# Get Uclidian distance
isoprov_dist(end_members, ref_data, dist_type = "ed")
#> $group1
#>     pb64   pb74   pb84        dist ref_groups ref_pb64 ref_pb74 ref_pb84
#> 6 17.910 15.645 38.003 0.005477226 Iglesiente 17.90500 15.64300 38.00200
#> 9 17.906 15.643 37.991 0.007842831 Iglesiente 17.89869 15.64346 37.98819
#> 8 17.910 15.639 37.988 0.009188576 Domusnovas 17.90510 15.64011 37.98030
#> 
#> $group2
#>        pb64     pb74     pb84        dist       ref_groups ref_pb64 ref_pb74
#> 2  18.97100 15.69400 39.02300 0.006491991        Taurus 1A 18.96813 15.68854
#> 3  18.98100 15.70300 39.05300 0.006633250        Taurus 1A 18.98700 15.70100
#> 1  18.98700 15.70200 39.04300 0.008062258        Taurus 1A 18.98700 15.70100
#> 27 18.99659 15.71336 39.07771 0.010183619        Taurus 1A 18.98800 15.70900
#> 5  18.97600 15.69500 39.02800 0.010599280        Taurus 1A 18.96813 15.68854
#> 24 18.94409 15.67360 38.93039 0.011136904 North South Axis 18.95320 15.67998
#> 30 18.97662 15.69148 39.01145 0.012303060            Thera 18.96800 15.68700
#> 4  18.98800 15.69500 39.02900 0.014060897        Taurus 1A 18.97893 15.69368
#> 26 19.00919 15.71697 39.09145 0.014279849        Taurus 1A 19.01105 15.70598
#> 29 18.98511 15.69619 39.01657 0.019574223            Thera 18.96800 15.68700
#> 31 18.99229 15.69722 39.02254 0.022002198        Taurus 1A 18.97893 15.69368
#> 10 19.03400 15.73200 39.38200 0.050163228  Mitterberg area 19.05700 15.68925
#>    ref_pb84
#> 2  39.02504
#> 3  39.05100
#> 1  39.05100
#> 27 39.08100
#> 5  39.02504
#> 24 38.92987
#> 30 39.01900
#> 4  39.03967
#> 26 39.08253
#> 29 39.01900
#> 31 39.03967
#> 10 39.39463
# Get Mass Fractionnation Corrected distance
isoprov_dist(end_members, ref_data, dist_type = "mfd")
#> $group1
#>     pb64   pb74   pb84    dist_sq ref_groups ref_pb64 ref_pb74 ref_pb84
#> 9 17.906 15.643 37.991 0.02100244 Domusnovas  17.9051 15.64011  37.9803
#> 6 17.910 15.645 38.003 0.03137617 Domusnovas  17.9051 15.64011  37.9803
#> 8 17.910 15.639 37.988 0.27828533     Oridda  17.9000 15.62100  37.9180
#> 
#> $group2
#>        pb64     pb74     pb84    dist_sq                      ref_groups
#> 29 18.98511 15.69619 39.01657 0.01133854                              BS
#> 5  18.97600 15.69500 39.02800 0.02314524                       Taurus 1A
#> 3  18.98100 15.70300 39.05300 0.05952820                       Taurus 1A
#> 1  18.98700 15.70200 39.04300 0.06574969                       Taurus 1A
#> 27 18.99659 15.71336 39.07771 0.06693365                              BS
#> 4  18.98800 15.69500 39.02900 0.07421114          South Eastern Anatolia
#> 30 18.97662 15.69148 39.01145 0.07846808                       Taurus 1A
#> 24 18.94409 15.67360 38.93039 0.09509082 Baia Borşa orefield, NW Romania
#> 31 18.99229 15.69722 39.02254 0.16471001          South Eastern Anatolia
#> 2  18.97100 15.69400 39.02300 0.20874441                       Taurus 1A
#> 26 19.00919 15.71697 39.09145 0.22758868          South Eastern Anatolia
#> 10 19.03400 15.73200 39.38200 0.95180218                  Alaverdi Kapan
#>    ref_pb64 ref_pb74 ref_pb84
#> 29 18.96500 15.67078 38.93704
#> 5  18.98800 15.70900 39.08100
#> 3  18.98800 15.70900 39.08100
#> 1  18.98700 15.70100 39.05100
#> 27 18.96500 15.67078 38.93704
#> 4  18.98000 15.68700 38.99100
#> 30 18.98700 15.70100 39.05100
#> 24 18.93300 15.65800 38.89200
#> 31 18.98000 15.68700 38.99100
#> 2  18.96813 15.68854 39.02504
#> 26 18.98000 15.68700 38.99100
#> 10 18.92620 15.58650 38.88880
# Get XGBOOST predicted resutls
isoprov_predict(end_members, model)
#> $group1
#>          group      prob
#> 425     Oridda 0.9991462
#> 406 Iglesiente 0.9971285
#> 408 Iglesiente 0.9970817
#> 399 Domusnovas 0.8086246
#> 407 Iglesiente 0.4738236
#> 397 Domusnovas 0.1697372
#> 
#> $group2
#>                                  group      prob
#> 3877                         Taurus 1A 0.9997184
#> 1155                             Thera 0.9996513
#> 3885                         Taurus 1A 0.9990544
#> 3786                Eastern Taurides_1 0.9872018
#> 2395 Baia Borşa orefield  NW Romania_1 0.9809610
#> 2396 Baia Borşa orefield  NW Romania_1 0.9737222
```

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-A.D.B+24" class="csl-entry">

Albarede, Francis, Gillan Davis, Janne Blichert-Toft, Liesel Gentelli,
Haim Gitler, Marine Pinto, and Philippe Telouk. 2024. “A New Algorithm
for Using Pb Isotopes to Determine the Provenance of Bullion in Ancient
Greek Coinage.” *Journal of Archaeological Science* 163 (March): 105919.
<https://doi.org/10.1016/j.jas.2023.105919>.

</div>

<div id="ref-A.D.G+24" class="csl-entry">

Albarede, Francis, Gillan Davis, Liesel Gentelli, Janne Blichert-Toft,
Haim Gitler, Marine Pinto, and Philippe Telouk. 2024. “Bullion Mixtures
in Silver Coinage from Ancient Greece and Egypt.” *Journal of
Archaeological Science* 162 (February): 105918.
<https://doi.org/10.1016/j.jas.2023.105918>.

</div>

<div id="ref-E.E.Y+19" class="csl-entry">

Eshel, T, Y. Erel, N. Yahalom-Mack, O. Tirosh, and A. Gilboa. 2019.
“Lead Isotopes in Silver Reveal Earliest Phoenician Quest for Metals in
the West Mediterranean.” *Proceedings of the National Academy of
Sciences* 116 (13): 6007–12. <https://doi.org/10.1073/pnas.1817951116>.

</div>

</div>
