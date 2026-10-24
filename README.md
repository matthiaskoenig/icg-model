[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5552397.svg)](https://doi.org/10.5281/zenodo.5552397)
[![Release](https://img.shields.io/github/v/release/matthiaskoenig/icg-model)](https://github.com/matthiaskoenig/icg-model/releases)
[![Documentation](https://github.com/matthiaskoenig/icg-model/actions/workflows/docs.yml/badge.svg)](https://matthiaskoenig.github.io/icg-model/)
[![CI](https://github.com/matthiaskoenig/icg-model/actions/workflows/ci.yml/badge.svg)](https://github.com/matthiaskoenig/icg-model/actions/workflows/ci.yml)
[![Docker](https://img.shields.io/docker/v/matthiaskoenig/icg?sort=semver&label=docker)](https://hub.docker.com/r/matthiaskoenig/icg)
[![Python](https://img.shields.io/badge/python-3.14-blue)](https://www.python.org)
[![Code license: MIT](https://img.shields.io/badge/code%20license-MIT-blue)](https://opensource.org/license/MIT)
[![Model license: CC BY-SA 4.0](https://img.shields.io/badge/model%20license-CC%20BY--SA%204.0-lightgrey)](https://creativecommons.org/licenses/by-sa/4.0/)

# Physiologically based pharmacokinetic (PBPK) model of indocyanine green

This repository provides the indocyanine green (ICG) physiologically based pharmacokinetic (PBPK) model of the ICG liver function tests: the distribution of ICG in the plasma, its hepatic uptake, its reflux into the plasma and its biliary excretion, with the effects of cirrhosis, hepatectomy, hepatic blood flow and bilirubin.

Documentation: <https://matthiaskoenig.github.io/icg-model/>

## Quick start

```bash
git clone https://github.com/matthiaskoenig/icg-model.git && cd icg-model
uv sync
uv run run_simulation pkdb_models.models.icg --action all -r results
```

## How to cite

> Köller, A. & König, M. (2026). *Physiologically based pharmacokinetic (PBPK) model of indocyanine green*. Zenodo. https://doi.org/10.5281/zenodo.5552397

## License

* Source Code: [MIT](https://opensource.org/license/MIT)
* Documentation: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

&copy; 2021-2026 Adrian Köller and Matthias König
