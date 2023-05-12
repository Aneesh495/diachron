# Chapter 12: Sensitivity and robustness

Twelve predeclared configurations each change one lever and refit, so stable
findings can be separated from prior- and calibration-dependent ones. The levers
are: coalescent population size at three values (tree-prior strength),
ascertainment on versus off, singleton inclusion, loan removal, a half-calibration
subset, no calibration, and the four observation models.

## Headline result: chronology is calibration-dependent

The root-age estimate shifts materially under the calibration-related levers
(halving the calibrations, removing calibration, and the tree-only model) and
under the geography model, while it is stable to the tree-prior strength and to
the ascertainment correction. This confirms hypothesis H4: without calibration the
time scale is not identifiable and the chronology is driven by the prior. The
practical consequence is that any dated conclusion from this platform must be
reported together with its calibration, and a chronology presented without stating
its calibration dependence would be misleading.

## Ascertainment and tree-prior stability

That the root age is stable to turning ascertainment off and to a factor-of-nine
change in the coalescent population size is itself a finding: the character
likelihood, not those priors, is driving the age estimate given a fixed
calibration. Stability here does not mean the ascertainment correction is
unnecessary, only that it does not dominate the dated age; it does change
likelihood values and the relative fit of constant-heavy versus variable-heavy
regimes, which matters for model comparison.

## Deliberate misspecification

Beyond the twelve levers, recovery includes datasets generated with borrowing and
missingness regimes the fitted model does not match, so the study tests whether a
confidently narrow posterior survives model violation. A posterior that narrows
below its honest width under misspecification is reported as over-confident.

## Reporting rule

A conclusion appears in the summary as stable only if it holds across the relevant
levers. Calibration-dependent and geography-dependent conclusions are labeled as
such. The full configuration table is written to
`runs/sensitivity/sensitivity.json`.
