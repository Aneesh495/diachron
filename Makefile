# Diachron build and research orchestration.
#
# The targets fall into three groups: environment and build (bootstrap, doctor,
# build), correctness and demo (test, demo), and the research campaign (cldf,
# attestations, likelihood-check, sample, diagnose, posterior-predictive,
# sensitivity, thesis) plus the acceptance/verify/release pipeline.

PYTHON ?= python3
PY := .venv/bin/python
PIP := .venv/bin/pip
MATURIN := .venv/bin/maturin
VENV := .venv

.DEFAULT_GOAL := help

.PHONY: help
help:
	@echo "Diachron make targets:"
	@echo "  bootstrap            create venv and install Python tooling"
	@echo "  doctor               report toolchain, CPU, storage, data availability"
	@echo "  build                compile the Rust core and install the extension"
	@echo "  test                 run Rust and Python correctness suites"
	@echo "  data                 build the labeled-synthetic corpus and manifests"
	@echo "  cldf                 ingest CLDF lexical data into typed tables"
	@echo "  attestations         build the dated/located attestation table"
	@echo "  likelihood-check     pruning vs enumeration oracle (gate G02)"
	@echo "  sample               run the MCMC campaign (four chains per model)"
	@echo "  diagnose             convergence and mixing diagnostics (gate G06)"
	@echo "  posterior-predictive posterior predictive checks (gate G07 support)"
	@echo "  sensitivity          predeclared sensitivity configurations (gate G08)"
	@echo "  recovery             known-truth recovery study (gate G05)"
	@echo "  predict              held-out concept/language prediction (gate G07)"
	@echo "  thesis               assemble thesis and research notes (gate G14)"
	@echo "  demo                 bounded end-to-end demonstration"
	@echo "  bench / stress       performance and load campaigns"
	@echo "  docs                 render diagrams and check docs"
	@echo "  acceptance           full campaign into a new evidence directory"
	@echo "  verify               independent verification of an evidence bundle"
	@echo "  release              clean-checkout release packaging"

.PHONY: bootstrap
bootstrap:
	$(PYTHON) -m venv $(VENV)
	$(PIP) install -r requirements-dev.txt

.PHONY: doctor
doctor:
	$(PY) -m diachron.cli.main doctor

.PHONY: build
build:
	$(MATURIN) develop --locked --release

.PHONY: test
test:
	cd rust/diachron-core && cargo test --locked --lib
	$(PY) -m pytest tests -q

.PHONY: data
data:
	$(PY) -m diachron.cli.main data build

.PHONY: cldf
cldf:
	$(PY) -m diachron.cli.main cldf $(CLDF_ARGS)

.PHONY: attestations
attestations:
	$(PY) -m diachron.cli.main attestations $(ATTESTATION_ARGS)

.PHONY: likelihood-check
likelihood-check:
	$(PY) -m diachron.cli.main likelihood-check --cases 1000

.PHONY: sample
sample:
	$(PY) -m diachron.cli.main sample --profile cpu $(SAMPLE_ARGS)

.PHONY: diagnose
diagnose:
	$(PY) -m diachron.cli.main diagnose

.PHONY: posterior-predictive
posterior-predictive:
	$(PY) -m diachron.cli.main ppc

.PHONY: sensitivity
sensitivity:
	$(PY) -m diachron.cli.main sensitivity

.PHONY: recovery
recovery:
	$(PY) -m diachron.cli.main recovery --datasets 200

.PHONY: predict
predict:
	$(PY) -m diachron.cli.main predict

.PHONY: thesis
thesis:
	$(PY) -m diachron.cli.main thesis build

.PHONY: demo
demo:
	$(PY) -m diachron.cli.main demo

.PHONY: bench
bench:
	$(PY) -m diachron.cli.main bench

.PHONY: stress
stress:
	$(PY) -m diachron.cli.main stress

.PHONY: docs
docs:
	$(PY) -m diachron.cli.main docs render

.PHONY: acceptance
acceptance:
	$(PY) -m diachron.cli.main acceptance --run-id $(RUN_ID) --profile $(PROFILE)

.PHONY: verify
verify:
	$(PY) -m diachron.cli.main verify --evidence $(EVIDENCE)

.PHONY: release
release:
	$(PY) -m diachron.cli.main release
