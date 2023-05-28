# ADR 0005: Derive one aggregate root-age calibration from attestations

## Status
Accepted (revised from an earlier per-tip design).

## Context
Dated attestations must enter the model as node-age calibrations to make the time
scale identifiable. An initial design placed one calibration on each attested
tip's parent node.

## Decision
Derive a single aggregate root-age calibration: the oldest reliable attestation
bounds the root age below, and a margin above it caps it. The MRCA of all tips
(the root) is special-cased so no ancestor walk is needed.

## Consequences
The time scale is identifiable without over-constraining the tree. The earlier
per-tip design placed forty tight constraints on internal nodes, which the sampler
could almost never satisfy simultaneously, driving acceptance to zero; it was also
slow because each constraint rebuilt the children map. The aggregate calibration
fixed both. The limitation is that a single root calibration is coarse; richer
designs with defended internal-node calibrations are possible but require per-node
evidence the synthetic corpus does not provide.
