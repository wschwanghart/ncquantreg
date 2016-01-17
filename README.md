# Non-crossing polynomial quantile regression

ncquantreg finds the coefficients of a polynomial p(x) of degree n that fits the data in vector x to the quantiles tau of y. 

ncquantreg(x,y) performs median regression (tau = 0.5) using a polynomial of degree n=1. 

ncquantreg(x,y,n,tau) fits a numel(tau) polynomials with degree n. The algorithm uses a stepwise multiple quantile regression estimation using non-crossing constraints (Wu and Liu, 2009). The approach is stepwise in a sense that a quantile function is estimated so that it does not cross with a function fitted in a previous step. The algorithm starts from the middle quantile (i.e. the one closest to 0.5) and than progressivly works through the quantiles with increasing distance from the middle.

ncquantreg(x,y,n,tau,pn,pv,...) takes several parameter name value pairs that control the algorithm and plotting. 

## Reference

Wu, Y., Liu, Y., 2009. Stepwise multiple quantile regression estimation using non-crossing constraints. Statistics and its Interface 2, 299–310.
