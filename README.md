[![Travis-CI Build Status](https://travis-ci.org/zeemkr/ncpen.svg?branch=master)](https://travis-ci.org/zeemkr/ncpen)

# ncpen

`ncpen` package fits the generalized linear models with various nonconvex penalties.
A unified algorithm is implemented in `ncpen` based on the convex concave procedure or difference convex algorithm and can be applied to most of the existing nonconvex penalties.
The available penalties in the pacakge are
the least absolute shrinkage and selection operator(LASSO),
smoothly clipped absolute deviation (SCAD),
minimax concave penalty (MCP),
truncated *l*1-penalty (TLP),
clipped LASSO (CLASSO),
sparse bridge (SRIDGE),
modified bridge (MBRIDGE),
and modified log (MLOG) penalties.
This package accepts a design matrix **X** and vector of responses **y**,
and produces the regularization path over a grid of values for the tuning parameter lambda.
Also provides user-friendly processes for plotting, selecting tuning parameters using cross-validation or generalized information criterion (GIC),
*l*2-regularization, penalty weights, standardization and intercept.
For a data set with many variables (high-dimensional data),
the algorithm selects relevant variables producing a parsimonious regression model.

For an example use, see [ncepn example](https://github.com/zeemkr/ncpen/blob/master/ncepn_example.pdf).

(This project is funded by Julian Virtue Professorship from Center for Applied Research at
Graziadio School of Business and Management at Pepperdine University.)

**References**
Kwon, S., Lee, S. and Kim, Y. (2015) <https://doi.org/10.1016/j.csda.2015.07.001>,
Lee, S., Kwon, S. and Kim, Y. (2016) <https://doi.org/10.1016/j.csda.2015.08.019>.
