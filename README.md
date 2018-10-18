Gamma GP: a polynomial approximation of Neutral Density
=======================================================

[![DOI](https://zenodo.org/badge/153560501.svg)](https://zenodo.org/badge/latestdoi/153560501)

The $\gamma_GP$ variable is a collection of polynomials that builds an approximate form of Neutral Density $\gamma_n$ over the world's oceans. [Sérazin (2011)](https://www.researchgate.net/publication/304023500_An_approximate_Neutral_Density_variable_for_the_world_oceans)Constructed this variable during his master thesis, under the supervision of [Trevor McDougall](https://research.unsw.edu.au/people/scientia-professor-trevor-mcdougall) and Paul Barker at CSIRO. The density variable $\gamma_GP$ was designed to replace the variable $\sigma_2$, currently used in ocean models. The world's oceans were decomposed into segments that covered each of the major oceanic basins. This allowed for accurate polynomial functions, $\gamma_{poly}$, to be fitted to a Neutral Density labelled version of the WOCE climatology. The global polynomial $\gamma_GP$ is the result of combining each of the $\gamma_poly$ functions. The polynomials have been built in terms of Practical Salinity and potential temperature to be able be applied in current ocean models. The polynomials $gamma_{poly}$ give 30 % of improvement compared to sigma_2 in the amount of fictitious diapycnal diffusivity on the Southern Ocean and at least 40 % on the other ocean basins.

Functions to compute $\gamma_GP$ are made available on this github project for the Matlab and Python language. This still experimental variable has not been fully approved and tested. As a result, care should be taken when interpreted the results coming from those functions.

Please cite the DOI associated to this project when you publish work including $\gamma_GP$.


