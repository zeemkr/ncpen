[![Travis-CI Build Status](https://travis-ci.org/zeemkr/ncpen.svg?branch=master)](https://travis-ci.org/zeemkr/ncpen)

# ncpen

This package fits the generalized linear models with various non-convex penalties.
A unified algorithm is implemented in **ncpen** based on the convex concave procedure or difference convex algorithm that can be applied to most of existing non-convex penalties.
The available penalties in the pacakge are
the least absolute shrinkage and selection operator(LASSO),
smoothly clipped absolute deviation (SCAD),
minimax concave penalty (MCP),
truncated *l*1-penalty (TLP),
clipped LASSO (CLASSO),
sparse bridge (SRIDGE),
modified bridge (MBRIDGE),
and modified log (MLOG) penalites.
This package accepts a design matrix **X** and vector of responses **y**,
and produces the regularization path ovaer a grid of values for the tuning parameter *lambda*.
Also provides user-friendly processes for plotting, selecting tuning parameters using cross-validation or generalized information criterion (GIC),
*l*2-regularization, penalty weights, standardization and intercept.

For an example, see [ncepn example](https://github.com/zeemkr/ncpen/blob/master/ncepn_example.pdf).
