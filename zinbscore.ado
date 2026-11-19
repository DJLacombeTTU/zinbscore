*! zinbscore 1.0.2  18nov2026
* FILENAME: zinbscore.ado
* Score test of a negative binomial model against a zero-inflated negative binomial model
* with covariates in the inflation probability (Jansakul and Hinde 2009, Communications in
* Statistics: Simulation and Computation 38(1): 92-108), with:
*   - NB2 (dispersion(mean)) and NB1 (dispersion(constant)) baselines; by default both are
*     fitted and the test is reported for the one with the larger log likelihood;
*   - the efficient information with every block expected (the dispersion block summed over y);
*   - one-step estimates and t statistics for each inflation variable, conditional on the others;
*   - the intercept-only signed statistic (positive: NB predicts too few zeros);
*   - an F(q, N - p) reference and an optional parametric bootstrap that refits only nbreg.
* Syntax:
*   zinbscore depvar [indepvars] [if] [in], inflate(varlist | _cons)
*       [dispersion(best|mean|constant) exposure(varname) offset(varname) reps(#) seed(#)]
* indepvars and the inflate() varlist may contain i.factor terms (no interactions). No weights.
* 1.0.1: factor levels are built as indicators directly (1.0.0 lost the first non-base level).
* 1.0.2: the one-step table is labeled from the original names (display only).

program define zinbscore, rclass
    version 15
    syntax varlist(min=1 numeric fv) [if] [in], INFlate(string) ///
        [ Dispersion(string) EXPosure(varname numeric) OFFset(varname numeric) ///
          REPS(integer 0) SEED(string) ]

    gettoken depvar xvars : varlist
    _fv_check_depvar `depvar'
    if "`dispersion'" == "" local dispersion "best"
    if !inlist("`dispersion'", "best", "mean", "constant") {
        di as err "dispersion() must be best, mean or constant"
        exit 198
    }
    if "`exposure'" != "" & "`offset'" != "" {
        di as err "specify exposure() or offset(), not both"
        exit 198
    }
    marksample touse
    if "`exposure'" != "" markout `touse' `exposure'
    if "`offset'"   != "" markout `touse' `offset'

    * offset
    tempvar off
    if "`exposure'" != "" quietly gen double `off' = ln(`exposure') if `touse'
    else if "`offset'" != "" quietly gen double `off' = `offset' if `touse'
    else quietly gen double `off' = 0 if `touse'
    markout `touse' `off'

    * count design: expand factor variables, drop base and omitted levels, then collinear
    * columns. The indicators are temporary variables of this program, so they survive.
    * Factor levels become indicator variables built here (fvrevar applied to a list of
    * single levels treated the first level as a base in version 1.0.0 and dropped it).
    local xv ""
    _zs_expand `xvars' if `touse'
    local xkeep `r(keep)'
    if "`xkeep'" != "" {
        local xraw ""
        foreach t of local xkeep {
            if regexm("`t'", "^([0-9]+)\.([A-Za-z_][A-Za-z0-9_]*)$") {
                local lev = regexs(1)
                local fvv = regexs(2)
                tempvar tv
                quietly gen double `tv' = (`fvv' == `lev') if !missing(`fvv')
                local xraw `xraw' `tv'
            }
            else if regexm("`t'", "^c\.([A-Za-z_][A-Za-z0-9_]*)$") {
                local cv = regexs(1)
                local xraw `xraw' `cv'
            }
            else if regexm("`t'", "^[A-Za-z_][A-Za-z0-9_]*$") local xraw `xraw' `t'
            else {
                di as err "zinbscore: `t' is not supported; use variables and i.factor terms (no interactions)"
                exit 198
            }
        }
        _zs_keep `xraw' if `touse', names(`xkeep')
        local xv `r(tvars)'
    }

    * inflation design
    local infl = trim("`inflate'")
    if "`infl'" == "_cons" local infl ""
    local gv ""
    local gnames ""
    if "`infl'" != "" {
        _zs_expand `infl' if `touse'
        local gkeep `r(keep)'
        if "`gkeep'" != "" {
            local graw ""
            foreach t of local gkeep {
                if regexm("`t'", "^([0-9]+)\.([A-Za-z_][A-Za-z0-9_]*)$") {
                    local lev = regexs(1)
                    local fvv = regexs(2)
                    tempvar tv
                    quietly gen double `tv' = (`fvv' == `lev') if !missing(`fvv')
                    local graw `graw' `tv'
                }
                else if regexm("`t'", "^c\.([A-Za-z_][A-Za-z0-9_]*)$") {
                    local cv = regexs(1)
                    local graw `graw' `cv'
                }
                else if regexm("`t'", "^[A-Za-z_][A-Za-z0-9_]*$") local graw `graw' `t'
                else {
                    di as err "zinbscore: `t' is not supported; use variables and i.factor terms (no interactions)"
                    exit 198
                }
            }
            _zs_keep `graw' if `touse', names(`gkeep')
            local gv `r(tvars)'
            local gnames `r(names)'
            markout `touse' `gv'
        }
    }
    if "`xv'" != "" markout `touse' `xv'

    quietly count if `touse'
    local N = r(N)
    capture assert `depvar' >= 0 & `depvar' == floor(`depvar') if `touse'
    if _rc {
        di as err "`depvar' must be a nonnegative integer count"
        exit 459
    }
    quietly count if `depvar' == 0 & `touse'
    local nzero = r(N)

    local k : word count `xv'
    local q : word count `gv'
    local q = `q' + 1
    local dfden = `N' - (`k' + 2)

    if "`dispersion'" == "best" local disps "mean constant"
    else local disps "`dispersion'"

    foreach d of local disps {
        tempname b_`d' LM_`d' Z_`d' O_`d' M_`d'
        capture quietly nbreg `depvar' `xv' if `touse', offset(`off') dispersion(`d') iterate(300)
        if _rc & _rc != 430 {
            di as err "nbreg, dispersion(`d') failed with rc = " _rc
            exit _rc
        }
        if e(converged) != 1 di as txt "(warning: nbreg, dispersion(`d') did not converge)"
        local ll_`d' = e(ll)
        matrix `b_`d'' = e(b)
        if colsof(`b_`d'') != `k' + 2 {
            di as err "unexpected e(b) from nbreg, dispersion(`d')"
            exit 498
        }
        mata: _zs_compute("`depvar'", "`xv'", "`gv'", "`off'", "`touse'", "`b_`d''", "`d'", ///
                          "`LM_`d''", "`Z_`d''", "`O_`d''", "`M_`d''")
        if missing(scalar(`LM_`d'')) {
            di as err "the information matrix is not positive definite (dispersion(`d'))"
            exit 498
        }
        if scalar(`M_`d'') < 1 - 1e-6 {
            di as err "the summation over y for the expected information is incomplete (dispersion(`d'))"
            exit 498
        }
        matrix rownames `O_`d'' = _cons `gnames'
        matrix colnames `O_`d'' = onestep se t
    }

    if "`dispersion'" == "best" {
        local chosen = cond(`ll_constant' > `ll_mean', "constant", "mean")
    }
    else local chosen "`dispersion'"

    * ---- bootstrap ------------------------------------------------------
    if `reps' > 0 {
        if "`seed'" != "" set seed `seed'
        tempvar mu g ystar
        tempname bs LMs Zs Os Ms
        quietly gen double `mu' = .
        mata: _zs_mu("`xv'", "`off'", "`touse'", "`b_`chosen''", "`mu'")
        local phi = exp(el(`b_`chosen'', 1, `k' + 2))
        quietly gen double `g' = .
        quietly gen double `ystar' = .
        local ge 0
        local ok 0
        forvalues r = 1/`reps' {
            if "`chosen'" == "mean" {
                quietly replace `g' = rgamma(1/`phi', `phi') if `touse'
                quietly replace `ystar' = rpoisson(`mu' * `g') if `touse'
            }
            else {
                quietly replace `g' = rgamma(`mu'/`phi', `phi') if `touse'
                quietly replace `ystar' = rpoisson(`g') if `touse'
            }
            quietly replace `ystar' = 0 if `touse' & missing(`ystar')
            capture quietly nbreg `ystar' `xv' if `touse', offset(`off') dispersion(`chosen') iterate(300)
            if _rc continue
            if e(converged) != 1 continue
            matrix `bs' = e(b)
            mata: _zs_compute("`ystar'", "`xv'", "`gv'", "`off'", "`touse'", "`bs'", "`chosen'", ///
                              "`LMs'", "`Zs'", "`Os'", "`Ms'")
            if missing(scalar(`LMs')) continue
            local ++ok
            if scalar(`LMs') >= scalar(`LM_`chosen'') local ++ge
        }
        local pboot = (1 + `ge') / (1 + `ok')
    }

    * ---- display ----------------------------------------------------------
    di as txt _n "Score test: negative binomial against zero-inflated negative binomial"
    di as txt "(Jansakul and Hinde 2009; identity link for the inflation probability)"
    di as txt "Number of obs = " as res `N' as txt "   zeros = " as res `nzero' ///
       as txt "   inflation parameters q = " as res `q' as txt " (constant + `=`q'-1')"
    di as txt "{hline 78}"
    di as txt "Baseline" _col(14) "log lik" _col(29) "LM" _col(38) "df" _col(44) "p chi2" ///
       _col(56) "p F" _col(66) "z (const)"
    di as txt "{hline 78}"
    foreach d of local disps {
        local lab = cond("`d'" == "mean", "NB2 (mean)", "NB1 (const)")
        local star = cond("`d'" == "`chosen'", "*", " ")
        di as txt "`lab'`star'" _col(12) as res %10.3f `ll_`d'' _col(24) %8.3f scalar(`LM_`d'') ///
           _col(37) %3.0f `q' _col(43) %8.4f chi2tail(`q', scalar(`LM_`d'')) ///
           _col(53) %8.4f Ftail(`q', `dfden', scalar(`LM_`d'')/`q') _col(65) %8.3f scalar(`Z_`d'')
    }
    di as txt "{hline 78}"
    if "`dispersion'" == "best" di as txt "* chosen baseline: larger log likelihood"
    if `reps' > 0 di as txt "Parametric bootstrap p (" as res `ok' as txt " replications): " as res %6.4f `pboot'
    if scalar(`Z_`chosen'') < 0 {
        di as txt "Note: the intercept statistic is negative: the baseline predicts more zeros than"
        di as txt "      observed, so a rejection does not indicate zero inflation."
    }
    di as txt _n "One-step estimates for the inflation probability (" ///
       cond("`chosen'" == "mean", "NB2", "NB1") " baseline), each conditional on the others"
    di as txt "{hline 60}"
    di as txt "Variable" _col(20) "one-step" _col(34) "std. err." _col(48) "t" _col(54) "P>|t|"
    di as txt "{hline 60}"
    * row labels from the original names (matrix rownames turn a lone first level into "2bn.f3")
    local rn _cons `gnames'
    local i 0
    foreach nm of local rn {
        local ++i
        di as txt abbrev("`nm'", 16) _col(16) as res %12.5f `O_`chosen''[`i', 1] ///
           _col(30) %12.5f `O_`chosen''[`i', 2] _col(44) %7.2f `O_`chosen''[`i', 3] ///
           _col(53) %6.3f 2*normal(-abs(`O_`chosen''[`i', 3]))
    }
    di as txt "{hline 60}"

    * ---- stored results ---------------------------------------------------
    return scalar N = `N'
    return scalar q = `q'
    return scalar df_r = `dfden'
    foreach d of local disps {
        return scalar ll_`d' = `ll_`d''
        return scalar lm_`d' = scalar(`LM_`d'')
        return scalar p_`d' = chi2tail(`q', scalar(`LM_`d''))
        return scalar pF_`d' = Ftail(`q', `dfden', scalar(`LM_`d'')/`q')
        return scalar z_int_`d' = scalar(`Z_`d'')
        return scalar minmass_`d' = scalar(`M_`d'')
    }
    return scalar ll = `ll_`chosen''
    return scalar lm = scalar(`LM_`chosen'')
    return scalar p = chi2tail(`q', scalar(`LM_`chosen''))
    return scalar p_F = Ftail(`q', `dfden', scalar(`LM_`chosen'')/`q')
    return scalar z_int = scalar(`Z_`chosen'')
    if `reps' > 0 {
        return scalar p_boot = `pboot'
        return scalar reps_ok = `ok'
    }
    return local dispersion "`chosen'"
    return local inflate "`gnames'"
    foreach d of local disps {
        return matrix onestep_`d' = `O_`d'', copy
    }
    return matrix onestep = `O_`chosen'', copy
end

* expand a varlist that may hold factor variables; drop base and omitted levels
program define _zs_expand, rclass
    syntax [varlist(default=none fv)] [if]
    local keep ""
    if "`varlist'" != "" {
        fvexpand `varlist' `if'
        foreach t in `r(varlist)' {
            if strpos("`t'", "b.") | strpos("`t'", "o.") continue
            local keep `keep' `t'
        }
    }
    return local keep `keep'
end

* drop columns collinear with the others and a constant; keep each survivor's original name
program define _zs_keep, rclass
    syntax varlist(numeric) [if], names(string)
    quietly _rmcoll `varlist' `if', forcedrop
    local kept `r(varlist)'
    local tvars ""
    local nms ""
    local j 0
    foreach t of local varlist {
        local ++j
        local pos : list posof "`t'" in kept
        if `pos' == 0 continue
        local nm : word `j' of `names'
        local tvars `tvars' `t'
        local nms `nms' `nm'
    }
    return local tvars `tvars'
    return local names `nms'
end

version 15
mata:

real matrix _zs_info_nb2(real matrix X, real colvector mu, real scalar r, real scalar minmass)
{
    real scalar k, i, m, sd, ymax, Edd, mass
    real colvector yy, lp, p, u
    real matrix I

    k = cols(X)
    I = J(k + 1, k + 1, 0)
    I[|1, 1 \ k, k|] = cross(X, r :* mu :/ (r :+ mu), X)
    Edd = 0
    minmass = 1
    for (i = 1; i <= rows(mu); i++) {
        m = mu[i]
        sd = sqrt(m + m^2 / r)
        ymax = ceil(m + 40 * sd + 100)
        yy = (0::ymax)
        lp = lngamma(yy :+ r) :- lngamma(r) :- lngamma(yy :+ 1) :+ r * ln(r / (r + m)) :+ yy :* ln(m / (r + m))
        p = exp(lp)
        mass = sum(p)
        if (mass < minmass) minmass = mass
        u = -r :* (digamma(yy :+ r) :- digamma(r) :+ ln(r / (r + m)) :+ 1 :- (r :+ yy) :/ (r + m))
        Edd = Edd + sum(p :* u :^ 2)
    }
    I[k + 1, k + 1] = Edd
    return(I)
}

real matrix _zs_info_nb1(real matrix X, real colvector mu, real scalar lnd, real scalar minmass)
{
    real scalar k, i, m, a, d, l1d, sd, ymax, mass
    real colvector yy, lp, p, dg, ue, ud, Eee, Eed, Edd
    real matrix I

    k = cols(X)
    d = exp(lnd)
    l1d = ln(1 + d)
    Eee = J(rows(mu), 1, 0)
    Eed = J(rows(mu), 1, 0)
    Edd = J(rows(mu), 1, 0)
    minmass = 1
    for (i = 1; i <= rows(mu); i++) {
        m = mu[i]
        a = m / d
        sd = sqrt(m * (1 + d))
        ymax = ceil(m + 40 * sd + 100)
        yy = (0::ymax)
        lp = lngamma(yy :+ a) :- lngamma(a) :- lngamma(yy :+ 1) :+ yy :* lnd :- (yy :+ a) :* l1d
        p = exp(lp)
        mass = sum(p)
        if (mass < minmass) minmass = mass
        dg = digamma(yy :+ a) :- digamma(a)
        ue = a :* (dg :- l1d)
        ud = -a :* (dg :- l1d) :+ yy :- d :* (yy :+ a) :/ (1 + d)
        Eee[i] = sum(p :* ue :^ 2)
        Eed[i] = sum(p :* ue :* ud)
        Edd[i] = sum(p :* ud :^ 2)
    }
    I = J(k + 1, k + 1, 0)
    I[|1, 1 \ k, k|] = cross(X, Eee, X)
    I[|1, k + 1 \ k, k + 1|] = cross(X, Eed)
    I[|k + 1, 1 \ k + 1, k|] = cross(X, Eed)'
    I[k + 1, k + 1] = sum(Edd)
    return(I)
}

real matrix _zs_X(string scalar xn, string scalar tousen, real scalar n)
{
    if (xn == "") return(J(n, 1, 1))
    return((st_data(., tokens(xn), tousen), J(n, 1, 1)))
}

void _zs_mu(string scalar xn, string scalar offn, string scalar tousen,
            string scalar bname, string scalar muname)
{
    real colvector off, mu
    real matrix X, b
    real scalar k

    off = st_data(., offn, tousen)
    X = _zs_X(xn, tousen, rows(off))
    b = st_matrix(bname)
    k = cols(X)
    mu = exp(X * b[1, (1..k)]' + off)
    st_store(., muname, tousen, mu)
}

void _zs_compute(string scalar yn, string scalar xn, string scalar gn, string scalar offn,
                 string scalar tousen, string scalar bname, string scalar disp,
                 string scalar lmn, string scalar zn, string scalar onen, string scalar massn)
{
    real colvector y, off, mu, f0, s, U, lrm, a, one, se
    real matrix X, G, b, Gt, Itt, Igg, Igt, V, Vi
    real scalar n, k, lnphi, r, d, l1d, minmass

    y = st_data(., yn, tousen)
    off = st_data(., offn, tousen)
    n = rows(y)
    X = _zs_X(xn, tousen, n)
    if (gn == "") G = J(n, 1, 1)
    else G = (J(n, 1, 1), st_data(., tokens(gn), tousen))
    b = st_matrix(bname)
    k = cols(X)
    lnphi = b[1, k + 1]
    mu = exp(X * b[1, (1..k)]' + off)
    minmass = 1
    if (disp == "mean") {
        r = exp(-lnphi)
        lrm = -ln(1 :+ mu :/ r)
        f0 = exp(r :* lrm)
        Gt = (X :* (-r :* mu :/ (r :+ mu)), -r :* (lrm :+ mu :/ (r :+ mu)))
        Itt = _zs_info_nb2(X, mu, r, minmass)
    }
    else {
        d = exp(lnphi)
        a = mu :/ d
        l1d = ln(1 + d)
        f0 = exp(-a :* l1d)
        Gt = (X :* (-a :* l1d), a :* l1d :- mu :/ (1 + d))
        Itt = _zs_info_nb1(X, mu, lnphi, minmass)
    }
    s = (y :== 0) :/ f0 :- 1
    U = cross(G, s)
    Igg = cross(G, (1 :- f0) :/ f0, G)
    Igt = cross(G, Gt)
    V = Igg - Igt * cholsolve(Itt, Igt')
    V = (V + V') / 2
    Vi = cholinv(V)
    st_numscalar(massn, minmass)
    if (hasmissing(Vi)) {
        st_numscalar(lmn, .)
        st_numscalar(zn, .)
        st_matrix(onen, J(cols(G), 3, .))
        return
    }
    one = Vi * U
    se = sqrt(diagonal(Vi))
    st_numscalar(lmn, U' * Vi * U)
    st_numscalar(zn, U[1] / sqrt(V[1, 1]))
    st_matrix(onen, (one, se, one :/ se))
}

end
