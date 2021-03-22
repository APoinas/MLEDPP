Asymptotic approximation of the MLE of DPPs
===========================================

R code accompanying the paper "Asymptotic approximation of the 
likelihood of stationary determinantal point processes" <https://arxiv.org/abs/2103.02310> with 
Frederic Lavancier <https://github.com/lavancier-f>.

The file "MLE_DPP.R" contains the main function computing the asymptotic approximation of the MLE
of common parametric families of stationnary DPPs. The file "Estimator Comparison.R" gives examples
of utilisation of this function as well as a way to reproduce the results of out paper.

Instructions
------------

Syntax
~~~~~~

TO BE COMPLETED

Dependencies
~~~~~~~~~~~~

This project depends on the following packages:

-  `Spatstats <https://github.com/spatstat/spatstat>`__
-  `stats`__
-  `Matplotlib httpmatplotlib.org`__

The following dependencies are optional, and only needed for parallelizing the task of testing the
performance of the various estimators on DPPs that we consider.

-  `foreach`__
-  `doSNOW`__
-  `parallel`__

