# lineqGPR
## Gaussian processes regression models with linear inequality constraints

**Updates:**
  The beta version 0.2.0 is now available at Github. It contains new implementations
  based on the MaxMod algorithm proposed in (Bachoc et al., 2020)


**Description:**
  *lineqGPR* is a package for Gaussian process interpolation, regression and
  simulation under linear inequality constraints based on (López-Lopera et
  al., 2017). The constrained models are given as objects with "lineqGP" S3
  class. Implementations according to (Maatouk and Bay, 2017) are also
  provided as objects with "lineqDGP" S3 class.

**Note:**
  *lineqGPR* was developed within the frame of the Chair in Applied
  Mathematics OQUAIDO, gathering partners in technological research (BRGM,
  CEA, IFPEN, IRSN, Safran, Storengy) and academia (CNRS, Ecole Centrale
  de Lyon, Mines Saint-Etienne, University of Grenoble, University of Nice,
  University of Toulouse) around advanced methods for Computer Experiments.
  
**Authors:** Andrés Felipe López-Lopera (Mines Saint-Étienne)
  with contributions from
  Olivier Roustant (Mines Saint-Étienne) and
  Yves Deville (Alpestat).
  
**Maintainer:** Andrés Felipe López-Lopera, <anfelopera@utp.edu.co>

**References**

  A. F. López-Lopera, F. Bachoc, N. Durrande and O. Roustant (2018),
  "Finite-dimensional Gaussian approximation with linear inequality constraints".
  \emph{SIAM/ASA Journal on Uncertainty Quantification}, 6(3): 1224–1255.
  [[link]](https://doi.org/10.1137/17M1153157)

  F. Bachoc, A. Lagnoux and A. F. Lopez-Lopera (2019),
  "Maximum likelihood estimation for Gaussian processes under inequality constraints".
  \emph{Electronic Journal of Statistics}, 13 (2): 2921-2969.
  [[link]](https://doi.org/10.1214/19-EJS1587)
  
  F. Bachoc, A. F. Lopez-Lopera and O. Roustant (2020),
  "Sequential construction and dimension reduction of Gaussian processes under inequality constraints".
  \emph{ArXiv e-prints}
  [[link]](https://arxiv.org/abs/2009.04188)

  Maatouk, H. and Bay, X. (2017),
  "Gaussian process emulators for computer experiments with inequality constraints".
  *Mathematical Geosciences*, 49(5): 557-582.
  [[link]](https://link.springer.com/article/10.1007/s11004-017-9673-2)

  Roustant, O., Ginsbourger, D., and Deville, Y. (2012),
  "DiceKriging, DiceOptim: Two R Packages for the Analysis of
  Computer Experiments by Kriging-Based Metamodeling and Optimization".
  *Journal of Statistical Software*, 51(1): 1-55.
  [[link]](http://www.jstatsoft.org/v51/i01/)
