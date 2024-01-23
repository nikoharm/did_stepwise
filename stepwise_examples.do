/*

The code below simulates data from a DID/event study design and provides examples of
how to use did_stepwise to estimate treatment effects from this data, including various
additional functionalities. did_stepwise implements the Stepwise DID estimator of 
Harmon (2023)

The simulation setup and associated code is due originally to Kirill Borusyak as supporting
material for Borusyak et al. (2023).

*/

*** Initialize

clear all
set seed 020617
set sortseed 050619

*** Simulate an appropriate data set

* Simulate data on 250 units across 6 timer periods

global T = 6
global I = 250

* Let the errors be AR(1) with an autocorrelation coefficient of rho (=1 implies Random Walk, =0 impliess spherical errors)

global rho 1

* Set total observations and create a unit and time identifer

loc obs=$I*$T
set obs `obs'
gen i = int((_n-1)/$T )+1 					
gen t = mod((_n-1),$T )+1					

* Simulate the time of first treatment as uniform on 2,3,...,7

gen Ei = 1+ ceil(runiform()*$T) if t==1	
bys i(t): replace Ei = Ei[1]

* Create a "time since treated" variable and a treatment dummy

gen K = t-Ei 								
gen D = K>=0 & Ei!=. 

* Define the treatment effect for each unit (heterogeneous across time here)

gen tau = cond(D==1, (t+2.5-$T), 0)

* Simulate errors depending on the value of rho

if $rho <1 {
gen eta=rnormal()*(1-$rho^2)^0.5
gen eps=rnormal() if t==1
}
else {
gen eta=rnormal()*(2/5)^0.5
gen eps=0 if t==1
}
by i: replace eps=$rho * eps[_n-1]+eta if _n!=1

* Generate the outcome variable according to the following model:

gen Y = -E + 3*t + tau*D + eps 

* Finally to illustrate clustering and covariates, we create some additional variables

gen X1=log(i)
gen X2=ceil(i/50)*50
gen clust=ceil(i/10)*10

*** Analyze the data using SWDD

* Estimating treatment effects at all horizons

did_stepwise Y i t Ei

* Also compute an aggregate effect across horizons

did_stepwise Y i t Ei, agg

* Produce some pretrend estimates 

did_stepwise Y i t Ei, pretrends(4)

* Alternatively add also the pretrend estimate corresponding to the baseline period (mechanically this equals zero)

did_stepwise Y i t Ei, pretrends(4) includepre1

* Condition on the continous covariate X1 (e.g. conditional parallel trends based on a linear outcome model)

did_stepwise Y i t Ei, contcov(X1)

* Condition also on the categorical covariate X2 (e.g. conditional parallel trends based on a linear outcome model)

did_stepwise Y i t Ei, contcov(X1) catcov(X2)

* Use clustering on the variable clust

did_stepwise Y i t Ei, cluster(clust)

*** The code below shows how to create a standard event study plot of SWDD with estimates using the event_plot package
*** (this requires the event_plot package to be installed)

* Estimate treatment effects and store estimates and variances

did_stepwise Y i t Ei

matrix blag=e(b)
matrix Vlag=vecdiag(e(V))

* Estimate pretrend coefficients and store estimates and variances

did_stepwise Y i t Ei, pretrends(4) includepre1

matrix blead=e(b)
matrix Vlead=vecdiag(e(V))

* Stack the treatment effects and pretreatment results

matrix b=blag,blead
matrix V=Vlag,Vlead

* Produce an event plot

event_plot b#V, stub_lag(tau_horizon#) stub_lead(tau_pretrend#)  default_look graph_opt(xtitle("Time since event") ytitle("Coefficients"))


