*! did_stepwise: 
*! Treatment effect estimation in differnce-in-difference designs using the Stepwise DID estimator of Harmon (2023)
*! The program is mostly a wrapper script for did_imputation which appropriately differences the data and creates
*! weights for the relevant ATEs.
*! Version: November 7, 2023
*! Author: Nikolaj Harmon
*! Citation:
*! Please cite both of the below.
*! Harmon, "Difference-in-difference and the efficient estimation of treatment effects" (2023)
*! Borusyak, Jaravel, and Spiess, "Revisiting Event Study Designs: Robust and Efficient Estimation" (2023)

cap program drop did_stepwise
program define did_stepwise, eclass sortpreserve

syntax varlist(min=4 max=4) [if] [aw /] [, iwtr(varname) minn(integer 30)  ///
	AVGEFFectsby(varlist)  ///
	CLUSter(varname) tol(real 0.000001) maxit(integer 100) verbose nose ///
	fast AGGregate pretrends(integer 0) ///
	CONTCOVariates(varlist) CATCOVariates(varlist) ///
	includepre1 ///
	]

	* Instruct Stata that the data should be preserved (if the program terminates prematurely)

	qui preserve
	
	* Show launch message program
	
	noi disp ""
	noi disp "Stepwise DID estimator for average treatment effects"
	noi disp "(refer to Harmon (2023) for estimator details)"
	noi disp ""
	
	* Check that did_imputation is installed

	capture {
	which did_imputation
	}
	
	if _rc!=0 {
		di as error "Package did_imputation not found; please install"
		di as error ""
		error 199
	}
		
	* Grab variable names
	
		tokenize `varlist'
	local Y `1'
	local i `2'
	local t `3'
	local ei `4'
	
	* Check that specified variable names do not clash with variable names used later

	foreach name in `varlist' {
		if strpos("`name'","horizon")!=0 {
				di as error "Ilegal variable name specified. With this command, specified variables cannot contain the word 'horizon'."
				di as error ""
				error 459
		}
		if strpos("`name'","pretrend")!=0 {
				di as error "Ilegal variable name specified. With this command, specified variables cannot contain the word 'pretrend'."
				di as error ""
				error 459
		}
		if "`name'"=="average" {
				di as error "Ilegal variable name specified. With this command, specified variables cannot be named 'average'."
				di as error ""
				error 459
		}
	}
	
	* Prepare appropriate locals for implementing some options later
	
	if "`avgeffectsby'"!="" {
	loc long_avgeffectsby avgeffectsby(`avgeffectsby')
	}
	
	if "`cluster'"!="" {
	loc long_cluster cluster(`cluster')
	}
	
	if "`contcovariates'"!="" {
	loc long_contcov timecontrols(`contcovariates')
	}
	
	if "`catcovariates'"!="" {
	loc icatcov
	foreach x in `catcovariates' {
	loc icatcov `icatcov' `t'#`x'
	}
	}
	
	if "`weight'"!="" {
	loc long_weight [`weight'=`exp']
	}
	
	* Mark observations to be included in the analysis
	
	tempvar insample
	
	* We will not include observations with missings or strings in the following variables
	
	loc varlist `Y' `t' `contcovariates' `iwtr'
	marksample touse1, 
	qui gen `insample'=`touse1'
	
	* We will not include observations with missings in the following variables (but strings are ok)
	
	loc varlist `i' `catcovariates' `cluster'
	marksample touse2,  strok
	qui replace `insample'=`insample'*`touse2'
	
	
	* Check for holes (if not doing fast version)
	
	if "`fast'"=="" {
	
	sort `insample' `Ymiss' `i' `t'
	qui count if (`i'==`i'[_n-1] & `t'!=`t'[_n-1]+1 & `insample'==1) 
	if `r(N)'!=0 {
	di as error "Holes detected in the data. Stepwise DID relies on the panel having no holes."
	di as error ""
	error 459
	}
	
	}

	* Ensure that data is sorted
	
	qui sort `insample' `i' `t'
	
	* Check that unit specific weight and covariates do not vary within unit if set (if not doing fast version)

	if "`fast'"=="" {

	if "`iwtr'"!="" {
	qui count if (`i'==`i'[_n-1] & `iwtr'!=`iwtr'[_n-1] & `insample'==1) 
	
	if `r(N)'!=0 {
	di as error "Unit-specific weight varies across time within unit. Unit-specific weights must be time-invariant"
	di as error ""
	error 459
	}
	
	}
	
	if "`contcovariates'"!="" {
	
	foreach x of varlist `contcovariates' {
	
	qui count if (`i'==`i'[_n-1] & `x'!=`x'[_n-1] & `insample'==1) 
	
	if `r(N)'!=0 {
	di as error "Continuous covariate `x' varies across time within unit. Covariates must be time-invariant"
	di as error ""
	error 459
	
	}
	
	}
	}
	
		if "`catcovariates'"!="" {
	
	foreach x of varlist `catcovariates' {
	
	qui count if (`i'==`i'[_n-1] & `x'!=`x'[_n-1] & `insample'==1) 
	if `r(N)'!=0 {
	di as error "Categorical covariate `x' varies across time within unit. Covariates must be time-invariant"
	di as error ""
	error 459
	
	}
	
	}
	}
	
	
	
	}
	
    * If no unit weight set then normalize everyone to a weight of 1
	
	if "`iwtr'"=="" {
	tempvar iwtr
	qui gen `iwtr'=1 if `insample'==1
	}
	
	
	
	
	
	*** Now we use an if-statement to branch on whether we are estimating treatment effects or pretrends
	
	
	
	
	if `pretrends'<=0 {
	
	
	
	
	
	** If estimating TREATMENT EFFECTS we...
	
	* Generate a differenced outcome variable
	
	tempvar Y_orig
	qui gen `Y_orig'=`Y'
	qui replace `Y'=.
	qui sort `i' `t'
	qui by `i': replace `Y'=`Y_orig'-`Y_orig'[_n-1] if `t'==`t'[_n-1]+1 & `insample'==1 & `insample'[_n-1]==1

	* Generate variables measuring number of periods since treatment and highest 
	* observed number of posttreatment periods
	
	tempvar K maxK
	
	qui gen `K' = `t'-`ei' if `insample'==1
	
	qui egen `maxK'=max(`K') if `insample'==1, by(`i')
	
	
	* If Stata weights are specified then check that treated observations all receive a weight of 1 (if not doing fast version)
	
	if "`fast'"=="" & "`exp'"!="" {
	
	qui count if (`exp'!=1 & `K'>=0 & `K'!=. & `insample'==1) 
	if `r(N)'!=0 {
	di as error "Some specified regression weights (`exp') are different from 1 for treatment unit in the post period."
	di as error ""
	error 459
	
	}
	
	}

	* Determine the highest number of observed posttreatment periods for any unit
	* with positive weight

	qui su `K' if `iwtr'>0 & `insample'==1, meanonly
	loc maxmaxK `r(max)'
	
	* Generate weights for SWDD estimate of ATT_h (for use in did_imputation)
 
	* It will be convenient to permanently name the weight variables horizon0, horizon1, etc.
	* but an unlucky user could have variables with these names already. Next bit
	* makes sure we store those variables and can restore them later
	
	loc exist_h_vars
	tempvar h_backup
	
	forvalues x = 0/`maxmaxK' {
 
	capture confirm variable horizon`x'
		if _rc==0 {
			loc exist_h_vars `exist_h_vars' `x'
			qui gen `h_backup'_`x'=horizon`x'
			qui drop horizon`x'
		}
 
	}
	
	* Now generate appropriate weights for effects at the different horizons
	
	loc horizon_vars

	forvalues x = 0/`maxmaxK' {
	
	qui su `iwtr' if `K'==`x' & `Y'!=. & `insample'==1, meanonly
	qui sort `i' `t'
	qui by `i': gen horizon`x' = ( (`K'>=0 & `K'<=`x')  & (`maxK'>=`x' ) & `Y'!=. ) * (`iwtr'/r(sum)) if `insample'==1 & `Y'!=.
	
	loc horizon_vars `horizon_vars' horizon`x'
	}
	
	* If specified, create weights appropriate to compute an aggregate treatment 
	* effect; name this variable average and deal with naming issues as before
	
	if "`aggregate'"!="" {
	
	tempvar avg_backup
	
	capture confirm variable average
		if _rc==0 {
			qui gen `avg_backup'=average 
			qui drop average
			loc exist_average yes
		}
 
	
	loc agghor average
	qui gen `agghor'=0 if `Y'!=. & `insample'==1
	
	tempname posttreats
	qui su `iwtr' if `K'>=0 & `K'!=. & `Y'!=.  & `insample'==1 , meanonly
	matrix `posttreats'=r(sum)
	
	forvalues x = 0/`maxmaxK' {
	
	qui su `iwtr' if `K'==`x' & `Y'!=. & `insample'==1, meanonly
	qui replace `agghor'=`agghor'+horizon`x'*r(sum)/`posttreats'[1,1] if `insample'==1 & `Y'!=.
	
	}
	
	}
	
	* Run did_imputation
	
	noi disp "Estimating average treatment effects for all possible horizons"
	noi disp "
	noi disp "Running did_imputation on first-differenced data with appropriate weights:" 
	noi disp "(did_imputation due to Borusyak, Jaravel, and Spiess (2023))"
	noi disp ""
	/*
	noi disp "Basic command being run (ignoring if-statements, weights and other detailed options):" 
	noi disp ""
	noi disp as input ". did_imputation `Y' `i' `t' `ei', wtr(`horizon_vars' `agghor') sum fe(t `icatcov') `long_contcov' `long_cluster'" 
	noi disp ""
	*/
	
	
	capture noisily did_imputation `Y' `i' `t' `ei' if `insample'==1 `long_weight' , wtr(`horizon_vars' `agghor') sum minn(`minn') ///
	`long_avgeffectsby' fe(`t' `icatcov')  ///
	`long_cluster' tol(`tol') maxit(`maxit') `verbose' `nose' `long_contcov'
	
	if _rc!=0 {
	disp as error""
	disp as error "did_imputation exited with above warning/error."
	di as error ""
	}
	

	* Restore the outcome variable to non-differenced version
	
	qui replace `Y'=`Y_orig'
	
	* Restore horizon variables and average variable (if applicable)
	
	qui drop `horizon_vars' `agghor'
	
	if "`exist_h_vars'"!="" {
		foreach x in  `exist_h_vars' {
				qui gen horizon`x'= `h_backup'_`x'
				qui drop `h_backup'_`x'
		}
	}
	
	if "`exist_average'"!="" {
				qui gen average= `avg_backup'
				qui drop `avg_backup'
	}
	
	* If we made it this far, there is no need to restore the data (we have cleaned up after ourselves)
	
	qui restore, not

	* If estimating treatment effects, the program terminates here
	
	
	
	
	}
	else {
		
		
		
		
	** If instead we produce PRETREND estimates we...
		
	* Determine whether a unit is ever treated
	
	tempvar D ever_treat
	qui gen `D'=(`t'>=`ei') if `insample'==1
	qui egen `ever_treat'=max(`D') if `insample'==1, by(`i')
	
	* Exclude all treated periods from the analysis sample
	
	qui replace `insample'=0 if `D'==1
	
	* Create duplicate observations for eventually treated (e.g. placebo treated preperiods)
	
	tempvar copy
	qui expand=2 if abs(`ei' - `t')<=`pretrends'+1 & `insample'==1 & `ever_treat'==1 , generate(`copy')
	
	* Create a new i variable so that the duplicate observations are treated as belonging to a different unit
	
	tempvar i_new 
	qui egen `i_new'=group(`i' `copy')
	
	* If no clustering is specified, make sure we eventually cluster on the old i when
	* estimating pretrends
	
		if "`long_cluster'"=="" {
	loc long_cluster cluster(`i')
	}
	
	* Create new t and ei variables so that time is in reverse and the placebo treated periods 
	* are marked as treatment periods (save the old values)
	
	tempvar t_new ei_new
	
	qui gen `t_new'=-`t'
	qui gen `ei_new'=99
	qui replace `ei_new'=-(`ei'-2) if `copy'==1 & `ever_treat'==1 
	
	* With this modified data, we then proceed similar to the treatment effect case:
	
	* Generate differenced outcome variable
	
	tempvar Y_orig
	qui gen `Y_orig'=`Y'
	qui replace `Y'=.
	qui sort `i_new' `t_new'
	qui by `i_new': replace `Y'=`Y_orig'-`Y_orig'[_n-1] if `t_new'==`t_new'[_n-1]+1 & `insample'==1 & `insample'[_n-1]==1
	
	* Generate variables measuring number of periods since treatment and highest 
	* observed number of posttreatment periods
	
	tempvar K maxK
	
	qui gen `K' = `t_new'-`ei_new' if `insample'==1
	
	qui egen `maxK'=max(`K') if `insample'==1, by(`i_new')
	
	* If Stata weights are specified then check that the placebo treated preperiod observations all receive a weight of 1 (if not doing fast version)
	
	if "`fast'"=="" & "`exp'"!="" {
	
	qui count if (`exp'!=1 & `K'>=0 & `K'!=. & `insample'==1) 
	if `r(N)'!=0 {
	di as error "Some specified regression weights (`exp') are different from 1 for an eventually treatment unit one of the pretrend periods"
	di as error ""
	error 459
	
	}
	
	}

	* Determine the highest number of observed posttreatment periods for any unit
	* with positive weight

	qui su `K' if `iwtr'>0 & `insample'==1, meanonly
	loc maxmaxK `r(max)'
	
	* Generate weights for SWDD estimation of a pretrend ATT_h (e.g. for h<0) (for use in did_imputation)
 
	* It will be convenient to permanently name the weight variables pretrend1, pretrend2, etc.
	* but an unlucky user could have variables with these names already. Next bit
	* makes sure we store those variables and can restore them later. Note that
	* the loop below starts at -1 so that we also save pretrend1 if it exists 
	
	loc exist_h_vars
	tempvar h_backup
	
	forvalues x = -1/`maxmaxK' {
		
		loc xp2=`x'+2
 
	capture confirm variable pretrend`xp2'
		if _rc==0 {
			loc exist_h_vars `exist_h_vars' `xp2'
			qui gen `h_backup'_`xp2'=pretrend`xp2'
			qui drop pretrend`xp2'
		}
 
	}
	
	* Now generate appropriate weights for effects at the different horizons. If the user
	* has asked us to exclude it, we report a pretrend1 which is identically zero
	* (and so we set the corresponding weights to be identically zero)
	
	if "`includepre1'"!="" {
	gen pretrend1=0
	loc horizon_vars pretrend1
	}
	else {
	loc horizon_vars
	}

	forvalues x = 0/`maxmaxK' {
		
				loc xp2=`x'+2
	
	qui su `iwtr' if `K'==`x' & `Y'!=. & `insample'==1, meanonly
	qui sort `i_new' `t_new'
    qui by `i_new': gen pretrend`xp2' = ( (`K'>=0 & `K'<=`x')  & (`maxK'>=`x' ) & `Y'!=. ) * (`iwtr'/r(sum)) if `insample'==1 & `Y'!=.
	
	loc horizon_vars `horizon_vars' pretrend`xp2'
	}
	
	* If specified, create weights appropriate to compute an aggregate pretreatment effect;
	* name this variable average and dealing with naming issues as before
	
	if "`aggregate'"!="" {
	
	tempvar avg_backup
	
	capture confirm variable average
		if _rc==0 {
			qui gen `avg_backup'=average 
			qui drop average
			loc exist_average yes
		}
 
	
	loc agghor average
	qui gen `agghor'=0 if `Y'!=. & `insample'==1
	
	tempname posttreats
	qui su `iwtr' if `K'>=0 & `K'!=. & `Y'!=.  & `insample'==1 , meanonly
	matrix `posttreats'=r(sum)
	
	forvalues x = 0/`maxmaxK' {
		
		loc xp2=`x'+2
	
	qui su `iwtr' if `K'==`x' & `Y'!=. & `insample'==1, meanonly
	qui replace `agghor'=`agghor'+pretrend`xp2'*r(sum)/`posttreats'[1,1] if `insample'==1 & `Y'!=.
	
	}
	
	}
	
	* Run did_imputation
	
	noi disp "PRETRENDS: Estimating pretrends back to `pretrends' periods before baseline (if available)"
	noi disp "
	noi disp "Running did_imputation on appropriately modified data and appropriate weights:" 
	noi disp "(did_imputation due to Borusyak, Jaravel, and Spiess (2023))"
	if "`includepre1'"!="" {
	noi disp ""
	noi disp "Note: Inclusion of a mechanically zero first pretrend requested;"
	noi disp "command will produce a corresponding error:" 
	noi disp ""
	}
	
	capture noisily did_imputation `Y' `i_new' `t_new' `ei_new' if `insample'==1 `long_weight'  , wtr(`horizon_vars' `agghor') sum minn(`minn') ///
	`long_avgeffectsby' fe(`t' `icatcov')  ///
	`long_cluster' tol(`tol') maxit(`maxit') `verbose' `nose' `long_contcov'
	
	if _rc!=0 {
	disp as error""
	disp as error "did_imputation exited with above warning/error."
	di as error ""
	}
	

	* Restore the outcome variable to non-differenced version
	
	qui replace `Y'=`Y_orig'
	
	* Restore horizon variables and average variable (if applicable)
	
	qui drop `horizon_vars' `agghor'
	
	if "`exist_h_vars'"!="" {
		foreach x in  `exist_h_vars' {
				qui gen pretrend`x'= `h_backup'_`x'
				qui drop `h_backup'_`x'
		}
	}
	
	if "`exist_average'"!="" {
				qui gen average= `avg_backup'
				qui drop `avg_backup'
	}
	
	* Drop duplicate observations that was created to estimate pretrends
	
	qui drop if `copy'==1
	
	* If we made it this far, there is no need to restore the data (we have cleaned up after ourselves)
	
	qui restore, not
	
	* If producing pretrends, then the program terminates here
	
	}

	
	
	
	
	end
