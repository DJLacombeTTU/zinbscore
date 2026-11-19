{smcl}
{* *! version 1.0.2  18nov2026}{...}
{title:Title}

{phang}
{bf:zinbscore} {hline 2} Score test of a negative binomial model against a zero-inflated
negative binomial model, with covariates in the inflation probability


{title:Syntax}

{p 8 17 2}
{cmd:zinbscore} {depvar} [{indepvars}] {ifin}{cmd:,}
{opt inf:late(varlist | _cons)}
[{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt inf:late(varlist)}}candidate variables for the inflation equation; {cmd:_cons} for none{p_end}
{synopt:{opt d:ispersion(string)}}{cmd:best} (default), {cmd:mean} (NB2) or {cmd:constant} (NB1){p_end}
{synopt:{opt exp:osure(varname)}}exposure variable{p_end}
{synopt:{opt off:set(varname)}}offset variable{p_end}
{synopt:{opt reps(#)}}parametric bootstrap replications for the p-value; default 0 (none){p_end}
{synopt:{opt seed(#)}}random-number seed for the bootstrap{p_end}
{synoptline}
{p 4 6 2}
{it:indepvars} and the {cmd:inflate()} varlist may contain {cmd:i.}{it:factor} terms but not
interactions. Weights are not allowed.


{title:Description}

{pstd}
{cmd:zinbscore} fits {cmd:nbreg} (NB2 and NB1 by default) and computes the score test of
Jansakul and Hinde (2009) of the negative binomial model against a zero-inflated negative
binomial model whose inflation probability depends on the {cmd:inflate()} variables through an
identity link. Only {cmd:nbreg} is fitted; no zero-inflated model is needed. The test is
two-sided: it detects too many or too few zeros relative to the baseline.

{pstd}
With {cmd:dispersion(best)} the test is reported for both baselines and the one with the larger
log likelihood is marked and used for the one-step table. For each inflation variable the table
gives a one-step estimate and a t statistic, conditional on the other variables. The
intercept statistic z (const) gives the direction: positive means the baseline predicts too few
zeros.


{title:Reading the results}

{pstd}
A rejection means the chosen negative binomial model predicts the zeros wrongly along these
variables. It does not show that never-users exist, and a flagged variable may drive the
dispersion rather than the zeros. Before calling a variable an inflation variable, compare
{cmd:nbreg}, {cmd:nbreg, dispersion(constant)}, {cmd:gnbreg} with the variable in
{cmd:lnalpha()} and {cmd:zinb} with the variable in {cmd:inflate()} by AIC
({cmd:estimates stats}); see the course tutorial.


{title:Stored results}

{pstd}
{cmd:zinbscore} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations{p_end}
{synopt:{cmd:r(q)}}inflation parameters, constant included{p_end}
{synopt:{cmd:r(lm)}, {cmd:r(p)}, {cmd:r(p_F)}}LM statistic, chi2 and F p-values, chosen baseline{p_end}
{synopt:{cmd:r(z_int)}}intercept statistic, chosen baseline{p_end}
{synopt:{cmd:r(ll)}}log likelihood of the chosen baseline{p_end}
{synopt:{cmd:r(lm_mean)} ...}the same results for each baseline fitted ({cmd:_mean}, {cmd:_constant}){p_end}
{synopt:{cmd:r(p_boot)}, {cmd:r(reps_ok)}}bootstrap p-value and usable replications, with {cmd:reps()}{p_end}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(dispersion)}}chosen baseline: {cmd:mean} or {cmd:constant}{p_end}
{synopt:{cmd:r(inflate)}}inflation variables used{p_end}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(onestep)}}one-step estimate, standard error and t, chosen baseline{p_end}
{synopt:{cmd:r(onestep_mean)} ...}the same for each baseline fitted{p_end}


{title:Reference}

{phang}
Jansakul, N., and J. P. Hinde. 2009. Score tests for extra-zero models in zero-inflated
negative binomial models. {it:Communications in Statistics: Simulation and Computation}
38(1): 92-108.


{title:Note}

{pstd}
Written for PFP 6320 at Texas Tech University. The expected information is used throughout,
with the dispersion block summed over the count distribution.
