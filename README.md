# Mortality-Rate-Hedging

This project is based on the model by (Luciano & Vigna, 2005). Data of mortality rate is from Human Mortality Database(https://www.mortality.org/) and National Bureau of Statistics of China(http://data.stats.gov.cn/english/). Zero-coupon rate is used for calculating the present value.

‘Main.m’ runs all the functions used in this paper and get a table including all the results.

‘rParameterCalibration.m’ calibrates the parameters of interest rate model, which are CIR(Cox, Ingersoll Jr, & Ross, 2005) and Vasicek(Vasicek, 1977) models.

‘qParameterCalibration.m’ calibrates the parameters of mortality rate model, which are OU(Uhlenbeck & Ornstein, 1930) and Feller(Feller, 1951) Models.

‘RiskAdjustment.m’ calibrates the risk premium, delta, in order to compensate those companies exposed to the longevity risk.

‘HedgingSim.m’ runs Monte Carlo simulation for 100,000 times to get the incomes of pension fund, of which the mean, std, skewness, Kurtosis, VaREs and ES are all calculated in ‘AllTests.m’, and the function ‘VaREs.m’ is called to calulate VaREs. 

‘CapPricing.m’ is called in ‘HedgingSim.m’ to get the price of caplets.

‘KQD.m’ is based on a new way to arrange hedges(Li & Luo, 2012). But unfortunately, it doesn’t work well ☹ so its code are commented in throughout all functions.


---------------------------------------------------------------------------------------------------------------------------------

Reference:

Cox, J. C., Ingersoll Jr, J. E., & Ross, S. A. (2005). A theory of the term structure of interest rates. In Theory of valuation (pp. 129-164): World Scientific.

Feller, W. (1951). Two singular diffusion problems. Annals of mathematics, 173-182. 

Li, J. S.-H., & Luo, A. (2012). Key q-duration: A framework for hedging longevity risk. ASTIN Bulletin: The Journal of the IAA, 42(2), 413-452. 

Luciano, E., & Vigna, E. (2005). Non mean reverting affine processes for stochastic mortality. 

Uhlenbeck, G. E., & Ornstein, L. S. (1930). On the theory of the Brownian motion. Physical review, 36(5), 823. 

Vasicek, O. (1977). An equilibrium characterization of the term structure. Journal of financial economics, 5(2), 177-188. 

