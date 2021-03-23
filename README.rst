Asymptotic approximation of the MLE of DPPs
===========================================

R code accompanying the paper `Asymptotic approximation of the likelihood of 
stationary determinantal point processes <https://arxiv.org/abs/2103.02310>`_ with 
`Frederic Lavancier <https://github.com/lavancier-f>`_.

The file "MLE_DPP.R" contains the main function computing the asymptotic approximation of the MLE
of common parametric families of stationnary DPPs. The file "Estimator Comparison.R" gives examples
of utilisation of this function as well as a way to reproduce the results of out paper.

Instructions
------------

Syntax
~~~~~~

The full syntax of the MLEDPP function is the following

MLEDPP = function(ppp, DPPfamily, startpar=NULL, sigma=NULL, edgecorr=FALSE, Trunc=50)

The main arguments are:

- ppp -> The observed point pattern, an object of class "ppp".
- DPPfamily -> The DPP family that is fitted to the data, it has to be either "Gauss", "Bessel", "Cauchy" or "Whittle-Matern".

The additional arguments are:

- startpar -> Optional. Initial value of alpha used in the optimization procedure. By default it is alpha_max/2 where alpha_max is the highest value of alpha for which the DPP is well-defined.
- sigma -> Must be specified for "Whittle-Matern" families only. The shape parametric of the family (called nu in `dppMatern <https://rdrr.io/cran/spatstat.core/man/dppMatern.html>`__)
- edgecorr -> Logical. If 'TRUE' and the observation window is rectangular, it computes the periodic edge correction. If 'TRUE' and the observation window isn't rectangular, it computes the experimental edge correction described in Section 5.3 of our paper.
- Trunc -> Optional. Truncation used for the computation of L_0 for all DPP families except Bessel.

The function returns an object of class "detpointprocfamily" corresponding to the DPP family fitted
on the observed point pattern by the asymptotic approximation of the maximum likelihood estimator.

Limitations
~~~~~~~~~~~

Currently, the approximated MLE is coded for stationnary DPPs on any window of R^2 with either a Gaussian-type kernel, a Bessel-type kernel with sigma=0, a Cauchy-type kernel with nu=1/2 or a Whittle-Matern-type kernel with a fixed shape parameter.

Dependencies
~~~~~~~~~~~~

This project depends on the following packages:

-  `spatstat <https://github.com/spatstat/spatstat>`__
-  stats

The following dependencies are optional, and only needed for parallelizing the task of testing the
performance of the various estimators on DPPs that we consider.

-  foreach
-  doSNOW
-  parallel

