{smcl}
{* *! version 1 2022-08-16}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "event_plot" "help event_plot"}{...}
{viewerjumpto "Syntax" "did_imputation##syntax"}{...}
{viewerjumpto "Description" "did_imputation##description"}{...}
{viewerjumpto "Options" "did_imputation##options"}{...}
{viewerjumpto "Weights" "did_imputation##weights"}{...}
{viewerjumpto "Stored results" "did_imputation##results"}{...}
{title:Title}

{pstd}
{bf:did_stepwise} - Treatment effect estimation in Difference-in-difference design using the Stepwise DID estimator of Harmon (2023)

{marker syntax}{...}
{title:Syntax}

{phang}
{cmd: did_stepwise} {it:Y i t Ei} [if]  [{help did_stepwise##weights:estimation weights}] [{cmd:,} {help did_stepwise##options:options}]
{p_end}

{synoptset 8 tabbed}{...}
{synopt : {it:Y}}outcome variable {p_end}
{synopt : {it:i}}variable for unique unit id{p_end}
{synopt : {it:t}}variable for calendar period {p_end}
{synopt : {it:Ei}}variable for unit-specific date of treatment (missing = never-treated) {p_end}

{phang}
{it: Note}: {cmd:did_stepwise} requires a recent version of {help did_imputation} as well as {help reghdfe} {p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:did_stepwise} implements the Stepwise DID estimator of Harmon (2023). The command is mostly a wrapper script for the command {help did_imputation} written by Borusyak et al. (2023). Before calling {help did_imputation}, {cmd:did_stepwise} applies appropriate differencing to the data and constructs weights necessary for computing average treatment effects. By default, {cmd:did_stepwise} will produce average treatment effects for all possible horizons after treatment, where averages are over all treated units observed at the corresponding horizon. See to the helpfile for {help did_imputation} for additional details of this command.
{p_end}

{pstd}
Please report bugs and/or unusual behavior to {browse "mailto:nikolaj.harmon@econ.ku.dk":nikolaj.harmon@econ.ku.dk}.
{p_end}

{marker options}{...}
{title:Options}

{dlgtab:Estimands}

{phang}{opt iwtr(varlist)}: variable specifying unit-specific (time-invariant) weights that will be applied when computing ATEs.
{p_end}

{phang}{opt agg:regate}: option to have the command also output a treatment effect aggregated over time (averaged across time horizon). Default is only to produce ATEs for different time horizons. If produced, the aggregate over time is weighted by the number of observations that experience each horizon in the sample (e.g. if there are few units observed 3 periods after treatment then treatment effects at this horizon receive small weight).  {p_end}

{phang}{opt pretrends(integer)}: If specified, the command will produce estimates of pretrends for the specified number of periods prior to baseline, instead of estimating treatment effects. The pretrend estimates can be interpreted as estimates of the average treatment effect 2, 3, 4, etc. periods BEFORE the treatment actually occurs (note that SWDD estimates of the treatment effect 1 period before treatment would always be identically zero because this is the baseline period). If the identifying assumptions hold, the true value of such preperiod treatment effects should be zero. Note that when estimating pretrends, the relevant preperiod observations are being used both as control and placebo treatment observations, which affects the reported total number of observations (see {help did_stepwise##results:stored results}) {p_end}

{phang}{opt includepre1}: If specified when estimating pretrends, the command will produce a pretrend estimate for the treatment effect 1 period before treatment. This estimate will mechanically equal zero and will induce a warning message from the estimation command. The option is useful if the pretend estimates are to be plotted, since such pretrend plots customarily include a mechanical 0 estimate at time -1 (e.g. at baseline). {p_end}


{dlgtab:Covariates}

{phang}{opt contcov:ariates(varlist)}: list of continuous, unit-specific (time-invariant) covariates to condition on (e.g. conditional parallel trends imposed via a linear outcome equation)

{pmore}{it:In the parlance of {help did_imputation}, these are "timecontrols". See corresponding help file for more discussion.}
{p_end}

{phang}{opt catcov:ariates(varlist)}: list of categorial (e.g. fixed effects), unit-specific (time-invariant) covariates to condition on (e.g. conditional parallel trends imposed via a linear outcome equation) {p_end}

{pmore}{it:In the parlance of {help did_imputation}, these are FEs formed by interacting the categorial variables with time. See corresponding help file for more discussion.}


{dlgtab:Standard errors}

{phang}{opt clus:ter(varname)}: cluster SE within groups defined by this variable. Default is {it:i}. {p_end}

{phang}{opt avgeff:ectsby(varlist)}: Use this option (and/or {opt leaveout}) if you have small cohorts of treated observations, and after reviewing
Section 4.3 of Borusyak et al. (2021). In brief, SE computation requires averaging the treatment effects by groups of treated observations.{p_end}
{phang2}- These groups should be large enough, so that there is no downward bias from overfitting. {p_end}
{phang2}- But the larger they are, the more conservative SE will be, unless treatment effects are homogeneous within these groups. {p_end}
{phang2}- The varlist in {opt avgeffectsby} defines these groups.{p_end}
{phang2}- The default is cohort-years {opt avgeffectsby(Ei t)}, which is appropriate for large cohorts.{p_end}
{phang2}- With small cohorts, specify coarser groupings: e.g. {opt avgeffectsby(K)} (to pool across cohorts) or 
{opt avgeffectsby(D)} (to pool across cohorts and periods when computing the overall ATT). {p_end}
{phang2}- The averages are computed using the "smart" formula from Section 4.3, adjusted for any clustering and any choice of {opt avgeffectsby}. {p_end}

{phang}{opt nose}: do not produce standard errors (much faster).

{pmore}{it:These option are inherited directly from {help did_imputation}. See corresponding help file for more discussion.}
{p_end}


{dlgtab:Computation and implementation}

{phang}{opt fast}: option to speed up the computation by skipping various consistency checks on the data. Not recommended.  {p_end}

{phang}{opt minn(#)}: the minimum effective number (i.e. inverse Herfindahl index) of treated observations,
below which a coefficient is suppressed and a warning is issued. 
The inference on coefficients which are based on a small number of observations is unreliable. The default is {opt minn(30)}.
Set to {opt minn(0)} to report all coefficients nevertheless.

{pmore}{it:This option is inherited directly from {help did_imputation}.}: Check corresponding help file for more discussion.
{p_end}

{phang}{opt tol(real)}, {opt maxit(integer)}: tolerance and the maximum number of iterations.
This affects the iterative procedure used to search for the weights underlying the estimator (to produce SE).
Defaults are 10^-6 and 100, resp. If convergence is not achieved otherwise, try increasing them.

{pmore}{it:This option is inherited directly from {help did_imputation}.}: Check corresponding help file for more discussion.
{p_end}


{marker weights}{...}
{title:Estimation weights}

{phang} 
{cmd: did_stepwise} allows aweights. Moreover, the weights for all of the included treatment units are required to be 1 in all post-treatment periods (i.e. the weights can only be used to reweight control observations). See {help did_imputation##weights:did_imputation with weights} for additional discussion.
 {p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:did_stepwise} inherits the following from {help did_imputation} and stores it in {cmd:e()}:

{synoptset 10 tabbed}{...}
{p2col 5 15 15 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}A row-vector of (i) the estimates, (ii) pre-trend coefficients, and (iii) coefficients on controls from Step 1:{p_end}
{pmore3}- If {opt horizons} is specified, the program returns {bf:tau}{it:h} for each {it:h} in the list of horizons.{p_end}
{pmore3}- If multiple {opt wtr} are specified, the program returns {bf:tau_}{it:v} for each {it:v} in the list of {opt wtr} variables. {p_end}
{pmore3}- Otherwise the single returned coefficient is called {bf:tau}. {p_end}
{pmore3}- In addition, if {cmd:pretrends} is specified, the command returns {bf:pre}{it:h} for each pre-trend coefficient {it:h}=1..{opt pretrends}. {p_end}
{pmore3}- And if {cmd:controls} is specified, the command returns the coefficients on those controls. (Estimated fixed effects, {cmd:unitcontrols}, and {cmd:timecontrols} are not reported in the {cmd:e(b)}.){p_end}
{synopt:{cmd:e(V)}}Corresponding variance-covariance matrix {p_end}
{synopt:{cmd:e(Nt)}}A row-vector of the number of treated observations used to compute each estimator {p_end}

{p2col 5 15 15 2: Scalars}{p_end}
{synopt:{cmd:e(N)}} the total number of observations. Note that if estimating pretrends, this observation total will double count the preperiod observations on which pretrends are being estimated since these observation effectively serve both as (placebo) treatment and control observations. (scalar) {p_end}
{synopt:{cmd:e(Nc)}} the number of control observations used in imputation (scalar) {p_end}
{synopt:{cmd:e(pre_F), e(pre_p), e(pre_df)}} if {opt pretrends} is specified, the F-statistic, pvalue, and dof for the joint test for no pre-trends {p_end}
{synopt:{cmd:e(Niter)}} the # of iterations to compute SE {p_end}

{p2col 5 15 15 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}} {cmd:did_imputation} {p_end}
{synopt:{cmd:e(droplist)}} the set of coefficients suppressed to zero because of insufficient effective sample size (see the {opt minn} option){p_end}
{synopt:{cmd:e(autosample_drop)}} the set of coefficients suppressed to zero because treatment effects could not be imputed for any observation (if {opt autosample} is specified) {p_end}
{synopt:{cmd:e(autosample_trim)}} the set of coefficients where the sample was partially reduced because treatment effects could not be imputed for some observations (if {opt autosample} is specified) {p_end}

{p2col 5 15 15 2: Functions}{p_end}
{synopt:{cmd:e(sample)}} Marks the estimation sample: all treated observations for which imputation was successful and the weights are non-zero for at least one coefficient + all non-treated observations used in Step 1 {p_end}




{title:References}

{phang}
Harmon, Nikolaj (2023). "Difference-in-differences and efficient estimation of treatment effects," Working paper.
{p_end}

{phang}
Borusyak, Kirill, Xavier Jaravel, and Jann Spiess (2023). "Revisiting Event Study Designs: Robust and Efficient Estimation," Review of Economic Studies.
{p_end}

{marker author}{...}
{title:Author}

{pstd}
Nikolaj A. Harmon (University of Copenhagen), {browse "mailto:nikolaj.harmon@econ.ku.dk":nikolaj.harmon@econ.ku.dk}

