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
