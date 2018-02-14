# ncpen

[![Travis-CI Build Status](https://travis-ci.org/zeemkr/ncpen.svg?branch=master)](https://travis-ci.org/zeemkr/ncpen)

We introduce an R package called ncpen for estimating generalized linear models (GLM)
with various nonconvex penalties. We consider a class of penalties which includes the least
absolute shrinkage and selection operator (LASSO), smoothly clipped absolute deviation
(SCAD), minimax concave penalty (MCP), moderately clipped LASSO (MCL), truncated
LASSO penalty (TLP), and modied bridge penalty (MBP). In this paper, we develop a
unied algorithm for penalized estimation using these penalties and the R package ncpen.
The proposed algorithm includes recent technical advances for the optimization of penalized
estimation; convex-concave procedure for nonconvex optimization, local quadratic approximation
(LQA), a descent directional modication of LQA for GLM, coordinate-wise descent
algorithm for quadratic `1-penalized problems, active set optimization for speed up of the
algorithm, warm-start strategy for pathwise optimization, etc. The R package ncpen provides
various functions and options for user-specic choices of an initial value, penalty terms,
`2-regularization, observation and penalty weights, standardization and intercept.
