[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5552397.svg)](https://doi.org/10.5281/zenodo.5552397)
[![GitHub Action](https://github.com/matthiaskoenig/icg-model/actions/workflows/python.yml/badge.svg)](https://github.com/matthiaskoenig/icg-model/actions/workflows/python.yml)
[![GitHub Action](https://github.com/matthiaskoenig/icg-model/actions/workflows/docker.yml/badge.svg)](https://github.com/matthiaskoenig/icg-model/actions/workflows/docker.yml)

# icg model
This repository provides the indocyanine green (ICG) physiologically based pharmacokinetic (PBPK) model of the ICG liver function tests: the distribution of ICG in the plasma, its hepatic uptake by OATP1B3 (inhibited by bilirubin), its reflux from the liver into the plasma and its biliary excretion by MRP2, with the effects of cirrhosis, hepatectomy, hepatic blood flow and bilirubin. The fitted parameters come from the adopted fit of the model and carry their uncertainty (`distrib:uncertainty` in the SBML).

The model is distributed as [SBML](http://sbml.org) format available from [`icg_body_flat.xml`](./models/icg_body_flat.xml) with 
corresponding [SBML4humans model report](https://sbml4humans.de/report?url=https%3A%2F%2Fraw.githubusercontent.com%2Fmatthiaskoenig%2Ficg-model%2Fmain%2Fmodels%2Ficg_body_flat.xml) and [model equations](./models/icg_body_flat.md).
The COMBINE archive is available from [`icg_model.omex`](./icg_model.omex).

![model overview](./figures/icg_model.png)

### Comp submodels
* **liver** submodel [`icg_liver.xml`](./models/icg_liver.xml) with [SBML4humans report](https://sbml4humans.de/report?url=https%3A%2F%2Fraw.githubusercontent.com%2Fmatthiaskoenig%2Ficg-model%2Fmain%2Fmodels%2Ficg_liver.xml) and [equations](./models/icg_liver.md).
* **whole-body** submodel [`icg_body.xml`](./models/icg_body.xml) with [SBML4humans report](https://sbml4humans.de/report?url=https%3A%2F%2Fraw.githubusercontent.com%2Fmatthiaskoenig%2Ficg-model%2Fmain%2Fmodels%2Ficg_body.xml) and [equations](./models/icg_body.md).


## How to cite
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5552397.svg)](https://doi.org/10.5281/zenodo.5552397)

> Köller, A., & König, M. (2021).
> *Physiologically based pharmacokinetic (PBPK) model of indocyanine green.*   
> Zenodo. [https://doi.org/10.5281/zenodo.5552397](https://doi.org/10.5281/zenodo.5552397)

## License

* Source Code: [MIT](https://opensource.org/license/MIT)
* Documentation: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

## Run simulations
### python
Clone the repository 
```bash
git clone https://github.com/matthiaskoenig/icg-model.git
cd icg-model
```

#### uv
Run the complete analysis with uv (https://docs.astral.sh/uv/getting-started/installation/):
```bash
uv run run_simulation pkdb_models.models.icg -a all -r results
```

#### pip
If you use pip install the package via
```bash
pip install -e .
```
Run the complete analysis in the environment via:
```bash
run_simulation pkdb_models.models.icg -a all -r results
```

### docker
Simulations can also be run within a docker container:

```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/icg:latest /bin/bash
```

Run the complete analysis:
```bash
uv run run_simulation pkdb_models.models.icg -a all -r /results
```
The results are written into the mounted `/results` folder on the host.

In case of permission issues with the mounted folder, adjust ownership and access rights with:
```bash
sudo chown $(id -u):$(id -g) -R "${PWD}/results"
sudo chmod 775 "${PWD}/results"
```

## Funding
Adrian Köller and Matthias König were supported by the Federal Ministry of Education and Research (BMBF, Germany) within the research network Systems Medicine of the Liver (LiSyM) by grant number 031L0054. Matthias König was supported by the Federal Ministry of Research, Technology and Space (BMFTR, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B).

© 2021-2026 Matthias König, [Systems Medicine of the Liver](https://livermetabolism.com)
