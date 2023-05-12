# Chapter 3: Evidence, annotation, and the character transform

## Data contract

Lexical data enters as typed records: concepts, languages (with Glottolog codes
and modern reference coordinates), lexemes (with original script and
transliteration provenance), and cognate judgments (with confidence, loan and
partial-cognate flags, annotator, and source). Dated attestations carry a dating
interval in years before present, a source reference and page, an observation
type (attestation, inscription, reconstruction, trade), a location with
uncertainty, a language-assignment confidence, and explicit conflicting-record
links. Referential integrity and value ranges are checked on ingestion.

## CLDF ingestion

The CLDF reader locates the Wordlist tables by their component URL and resolves
each role-bearing column by its CLDF property, so it works across datasets that
name columns differently. It records the dataset license and DOI. Restricted
lexical collections are ingested from a directory acquired under their own
license and are never redistributed; the public release ships an acquisition
recipe and derived artifacts the license permits.

## The character transform and its audit trail

The transform from judgments to the character matrix is explicit and reversible
on paper. Each cognate class within a concept becomes one binary character whose
state for a language is present, absent, or missing by the rules in Chapter 4.
Every present and absent cell carries an evidence link back to the lexemes and
judgment sources that produced it, so an audit sample can be traced from a matrix
cell to the underlying annotation. Singleton classes are excluded and counted;
loans are recorded and, in the borrowing sensitivity analysis, removed from the
inherited characters to measure their effect.

## Evidence profile

The primary study corpus has 112 eligible varieties and 110 concepts with 11333
observed cells and 140 dated attestations across four observation types; the
second study has 36 varieties and 90 concepts with 60 attestations and its own
calibration regime. Both corpora are labeled synthetic known-truth: their history
is known, which is what the recovery and coverage studies need, and the synthetic
flag propagates into every artifact and view so they are never presented as
observed. The real-data path ingests CLDF collections of comparable scale when a
license permits.

## Correlated missingness and ascertainment

The transform preserves the information needed to test whether missingness,
replacement rate, and concept selection correlate with age or geographic
sampling. The corpora deliberately include uneven coverage and region-biased
attestation sampling so inference must reflect that bias rather than assume even
coverage.
