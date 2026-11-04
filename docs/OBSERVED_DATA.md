# Observed data inputs

The primary lexical study uses IE-CoR v1.2.1 at commit
`dd1ebc997bc38c824fb18a1ef234bbf2c229420b`. The independent lexical study
uses a bounded Austronesian subset of ABVD v0.3 at commit
`81748db0c5fcf3bf1d8c8372b86df75c375e1f22`. Both wordlists are CC BY 4.0.
Diachron records Glottolog 5.3 separately, including its CC BY 4.0 license and
release DOI. The lexical repositories and Glottolog classify at different
granularities. The Arrow crosswalk keeps the source IDs, labels and catalog
parent separately. Future catalog changes do not rewrite source IDs.

- [IE-CoR v1.2.1](https://github.com/lexibank/iecor/releases/tag/v1.2.1) and [dataset paper](https://doi.org/10.1038/s41597-025-05445-3)
- [ABVD v0.3](https://doi.org/10.5281/zenodo.23558745) and [source database](https://abvd.eva.mpg.de/austronesian/)
- [Glottolog 5.3](https://doi.org/10.5281/zenodo.18840935) and [release downloads](https://glottolog.org/meta/downloads)

## Acquire and freeze

Run from the Diachron checkout after `make bootstrap` and `make build`.
The freeze checks byte counts and SHA-256 against the committed source manifests
before parsing. The output directory must be new.

```bash
mkdir -p data/cache
git clone --depth 1 --branch v1.2.1 https://github.com/lexibank/iecor.git data/cache/iecor
git clone --depth 1 --branch v0.3 https://github.com/lexibank/abvd.git data/cache/abvd
curl -fL --max-time 60 \
  https://cdstar.eva.mpg.de/bitstreams/EAEA0-608B-9919-A962-0/glottolog_languoid.csv.zip \
  -o data/cache/glottolog53.zip
make cldf CLDF_ARGS='freeze --study iecor --source data/cache/iecor --catalog data/cache/glottolog53.zip --out runs/observed/iecor'
make cldf CLDF_ARGS='freeze --study abvd --source data/cache/abvd --catalog data/cache/glottolog53.zip --out runs/observed/abvd'
```

IE-CoR retains varieties with judgments for at least 80 meanings. ABVD retains
Austronesian varieties with judgments for at least 100 meanings and a
Glottocode, then takes the highest-coverage variety per Glottocode. It selects
40 of those using a fixed SHA-256 ordering. The sampled set is for bounded
lexical comparisons and carries no Indo-European calibrations. The selection
was set before measuring any predictive or historical endpoints.

| Observed study | Varieties | Concepts | Lexical cells | Lexically missing | Uncoded lexical cells | Loan judgments | Doubt flags |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| IE-CoR v1.2.1 | 152 | 170 | 24,989 | 851 | 0 | 2,255 | 0 |
| ABVD v0.3 subset | 40 | 210 | 7,355 | 1,045 | 1,169 | 70 | 129 |

IE-CoR excludes eight low-coverage varieties and includes 44 historical
varieties. ABVD has 386 eligible duplicate Glottocode records and omits 716
eligible representatives under the bounded-sample rule. Ten IE-CoR forms and
1,031 ABVD forms have multiple judgment rows. That count covers multiple
assignments and does not imply every such assignment is a contradiction.

Original script and lexical values are retained as supplied. Display forms use
NFC normalization. Cognate classes come from source judgments. For IE-CoR, the
loan flag is derived from membership in a class listed in `loans.csv`. For ABVD,
the source flag is attached to the form and carried through to its judgments.
A source doubt flag is categorical and carries no calibrated probability.
The local Arrow tables retain source licensing and identity in schema metadata.
Raw wordlists and Arrow tables stay outside the public source release.

## Dated and located inscriptions

The Epigraphic Database Heidelberg cohort uses the first 250 XML filenames in
`inscriptions/1/1` at the pinned EDH snapshot. The reader retains records with a
Latin edition, a two-sided date, a single named origin locality, and bibliography.
The resulting 175 entries each link to a distinct EDH record and carry its
bibliography, date bounds, and a source section locator. Of those, 61 also have a
unique EDH gazetteer reference point matched by ancient name and Trismegistos ID.
An ambiguous point is left unset. The other 75 source entries were excluded for
missing or ambiguous locality (59), incomplete date (14), or edition language (2).

The date display uses year units before 1950. The original date strings are
retained, including month and day where available. The reader retains multiple
structured dates as linked alternatives. Competing dates present only in prose
require individual review. Location precision is qualitative. An inscription
origin records where an object is reported to originate and cannot locate a
speech community. The cohort is a Latin-heavy epigraphic slice. Its records are
indexed to the IE-CoR Latin variety without asserting that all these inscriptions
sample that variety's lexical or temporal distribution. It provides inspection
evidence and is excluded from automatic root-age calibration.

The derived attestation table carries EDH's CC BY-SA 4.0 terms and attribution.
It contains bibliographic citations and provenance with no inscription text.

- [Epigraphic Database Heidelberg](https://edh.ub.uni-heidelberg.de/)
- [Pinned EDH snapshot](https://github.com/epigraphic-database-heidelberg/data/tree/45f166654ab4551a1954617a0df8e56fa7724ccb)
- [EDH data reuse guide](https://doi.org/10.11588/heidok.00026599)

```bash
mkdir -p data/cache/edh
edh_base=https://raw.githubusercontent.com/epigraphic-database-heidelberg/data/45f166654ab4551a1954617a0df8e56fa7724ccb
jq -r '.files | keys[]' data/manifests/edh.json |
  xargs -P 6 -I '{}' curl -fsSL --retry 2 --max-time 30 \
    "$edh_base/inscriptions/1/1/{}" -o 'data/cache/edh/{}'
curl -fsSL --max-time 60 "$edh_base/geography/edhGeographicData.json" -o data/cache/edh-geography.json
make attestations ATTESTATION_ARGS='edh --source data/cache/edh --geography data/cache/edh-geography.json --out runs/observed/attestations.json'
```

IE-CoR's lexical tip distributions use years before 2000 CE. The EDH
derived records retain their independent 1950 reference year. The
workbench adds 50 years to EDH ages on the linked IE-CoR timeline and
keeps the original 1950-based interval in the inspector. The conversion
is for display. EDH ages do not enter the IE-CoR posterior.
[IE-CoR dataset paper](https://doi.org/10.1038/s41597-025-05445-3).
