# Note 7: Convergence

Convergence is where it is easiest to fool oneself, so the platform reports it
conservatively and refuses to certify a run that does not meet the thresholds.

The scalar diagnostics are the rank-normalized split-Rhat and the Geyer
initial-monotone-sequence ESS, with bulk and tail variants. They behave correctly
on known cases: iid draws give ESS near the sample size and Rhat near one; an
AR(1) series with autocorrelation 0.9 gives an ESS near the analytic value of
`n(1-phi)/(1+phi)`; chains at different means give a large Rhat. These checks are
in the test suite, so a regression in the estimator would be caught.

The hard part is topology. A scalar summary can look settled while the discrete
tree has not mixed, so the platform also reports the number of distinct
topologies, the switch rate between consecutive draws, and the clade agreement
between independent chains. On the diagnostic-profile run of the 112-tip primary
corpus, the scalar endpoints had Rhat between 1.1 and 2.8 and ESS in the tens,
and the run was reported as not converged. The topology was exploring actively
(60 to 72 distinct trees per chain, switch rate near 0.55, clade agreement 0.76),
which tells us the chains are moving but have not yet agreed.

The honest conclusion recorded here is a negative one: a short run of a
hundred-tip tree space does not converge, and the acceptance gate correctly
refuses to present its draws as a posterior. More iterations are required, and no
amount of scalar-diagnostic massaging substitutes for them.
