# ADR 0012: Separate CPU diagnostic and full acceptance profiles

## Status
Accepted.

## Context
Quick feedback during development needs short chains, but a short chain on a large
tree does not converge. If a quick run could set the acceptance flag, an
unconverged result could masquerade as certified.

## Decision
Declare two named profiles. The CPU diagnostic profile uses short chains for
feedback and cannot set the acceptance flag. The full acceptance profile uses long
chains, writes a sealed evidence directory, and is verified independently. The
convergence criteria (split-Rhat below 1.01, ESS above 400 per endpoint) are the
same for both; only the budget differs.

## Consequences
A diagnostic run is honestly labeled as diagnostic and its non-convergence is
reported rather than hidden. The acceptance campaign is resumable and checkpointed,
so the long chains survive interruption. The cost is that the full campaign is
expensive on a hundred-tip tree; this is inherent to the problem, and the platform
reports non-convergence rather than shortening the campaign to force a pass.
