[![Travis-CI Build Status](https://travis-ci.org/zeemkr/ncpen.svg?branch=master)](https://travis-ci.org/zeemkr/ncpen)

# ncpen
`ncpen` package fits the generalized linear models with various nonconvex penalties.
Supported regression models are Gaussian (linear), binomial Logit (logistic),
multinomial Logit, Poisson and Cox proportional hazard.
A unified algorithm is implemented based on the convex concave procedure
and the algorithm can be applied to most of the existing nonconvex penalties.
The algorithm also supports convex penalty:
least absolute shrinkage and selection operator (LASSO).
Supported nonconvex penalties include
smoothly clipped absolute deviation (SCAD),
minimax concave penalty (MCP), truncated LASSO penalty (TLP),
clipped LASSO (CLASSO), sparse ridge (SRIDGE),
modified bridge (MBRIDGE) and modified log (MLOG).
This package accepts a design matrix **X** and vector of responses **y**,
and produces the regularization path over a grid of values for the tuning parameter lambda.
Also provides user-friendly processes for plotting, selecting tuning parameters using cross-validation or generalized information criterion (GIC),
*l*2-regularization, penalty weights, standardization and intercept.
For a data set with many variables (high-dimensional data),
the algorithm selects relevant variables producing a parsimonious regression model.

Related research paper can be found at [ncpen paper](http://arxiv.org/abs/1811.05061).
A recent manual is avaialbe at [ncpen manual](https://github.com/zeemkr/ncpen_resources/blob/master/ncpen.pdf) and for
an example use, see [ncepn example](https://github.com/zeemkr/ncpen_resources/tree/master/example_mortgage).

(This research is funded by Julian Virtue Professorship from Center for Applied Research at Pepperdine
Graziadio Business School and the National Research Foundation of Korea (NRF) funded
by Korean government (No. 2017R1C1B2010113 and 2017R1D1A1B03031239).)

**Authors**

Dongshin Kim, Sunghoon Kwon, Sangin Lee

**References**
* Kim, D., Lee, S. and Kwon, S. (2018). A unified algorithm for the non-convex penalized estimation: The ncpen package <http://arxiv.org/abs/1811.05061>.
* Kwon, S., Lee, S. and Kim, Y. (2015) <https://doi.org/10.1016/j.csda.2015.07.001>,
* Lee, S., Kwon, S. and Kim, Y. (2016) <https://doi.org/10.1016/j.csda.2015.08.019>.
