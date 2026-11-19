# zinbscore

A Stata command for the score test of a negative binomial model against a zero-inflated
negative binomial model, with covariates in the inflation probability (Jansakul and Hinde
2009). It tests every candidate for `inflate()` at once using only negative binomial fits,
so it runs even when `zinb` would not converge, and it reports a one-step *t* statistic for
each candidate to guide which variables belong in the inflation equation.

## Installation

In Stata 15 or later:

```stata
net install zinbscore, from("https://raw.githubusercontent.com/DJLacombeTTU/zinbscore/main/") replace
help zinbscore
```

To remove it: `ado uninstall zinbscore`. Rerun the install line with `replace` to update.

Without installing, you can also put `zinbscore.ado` and `zinbscore.sthlp` in your working
folder; Stata finds commands there.

## Syntax

```stata
zinbscore depvar [indepvars] [if] [in], inflate(varlist | _cons)
    [dispersion(best|mean|constant) exposure(varname) offset(varname) reps(#) seed(#)]
```

| Option | Description |
|---|---|
| `inflate()` | Candidate variables for the inflation equation (factor variables allowed, no interactions); `_cons` for a constant only |
| `dispersion()` | `best` (default): fit NB2 and NB1 and test under the one with the larger log likelihood; `mean`: NB2 (`nbreg`'s default); `constant`: NB1 |
| `exposure()`, `offset()` | As in `nbreg` |
| `reps()`, `seed()` | Optional parametric bootstrap p-value (refits only `nbreg`) |

Weights are not supported.

## Example

```stata
zinbscore loans expshock homeowner age, inflate(hascard incvol age expshock commute jobchg)
```

The output has two parts:

1. **The joint test** for each baseline: log likelihood, LM statistic (chi-squared with q
   degrees of freedom, q = candidates plus the constant), chi-squared and F p-values, and the
   intercept statistic z (const). The chosen baseline is marked with an asterisk.
2. **One-step estimates** of the inflation coefficients on the probability scale, with
   standard errors and *t* statistics, each conditional on the other candidates.

## Interpreting the results

- **The joint test** says whether the negative binomial misses the zeros somewhere. It does
  not say where.
- **z (const)** positive means more zeros than the baseline predicts. When it is negative,
  the command prints a note: the baseline already predicts too many zeros, so a rejection
  points to something other than zero inflation.
- **The one-step *t*** (|t| > 1.96) flags candidates whose high values go with more zeros
  than the negative binomial predicts.
- **A flag is not proof of a never-zero group.** A group with a low but positive rate, or a
  dispersion that changes with a covariate, also produces excess zeros tracking a variable,
  and is flagged the same way. Each *t* also has a 5% false-flag rate, so several noise
  candidates raise the chance of at least one false flag.
- **The baseline matters.** Testing under NB2 when the data are NB1 (or the reverse) can
  reject "no inflation" when there is none. The default `dispersion(best)` guards against this.

A suggested workflow: fit `zinb` with the full list of theory candidates only if something
is flagged, then compare the average marginal effects on the expected count
(`margins, dydx(*) predict(n)`) with those from `nbreg` with every variable in the count
equation. When they agree, the choice of inflation variables did not change the answer.

## Stored results

`r(N)`, `r(q)`, `r(df_r)`; for the chosen baseline `r(lm)`, `r(p)`, `r(p_F)`, `r(z_int)`,
`r(ll)` and the matrix `r(onestep)`; the same for each baseline fitted with suffixes `_mean`
and `_constant`; `r(dispersion)` (the chosen baseline) and `r(inflate)`; with `reps()`,
`r(p_boot)` and `r(reps_ok)`. See `help zinbscore`.

## Method

The score for each observation is $s_i = 1\{y_i = 0\}/\hat f_0(x_i) - 1$, where
$\hat f_0$ is the negative binomial probability of a zero. The statistic is
$LM = U'\hat V^{-1}U$ with $U = \sum_i z_i s_i$ and $\hat V$ the efficient information, every
block of it expected information (the dispersion block computed by summation over y). Under
no zero inflation LM is asymptotically chi-squared with q degrees of freedom. The one-step
estimates are $\hat V^{-1}U$.

## Validation

The command was checked against an independent Python implementation on seven simulated
datasets (NB2 and NB1 counts, no, weak and strong inflation, a factor candidate, and
n = 300): 192 of 192 checks agreed. The Python implementation was checked against numerical
derivatives of the zero-inflated likelihood and by Monte Carlo simulation of the test's size
and power.

## Reference

Jansakul, N., and J. P. Hinde. 2009. Score tests for extra-zero models in zero-inflated
negative binomial models. *Communications in Statistics: Simulation and Computation*
38(1): 92-108. https://doi.org/10.1080/03610910802421632

## Author

Donald J. Lacombe, School of Financial Planning, Texas Tech University.

## License

MIT. See [LICENSE](LICENSE).
