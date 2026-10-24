# Indocyanine green

Physiologically based pharmacokinetic (PBPK) model of the liver function tests with indocyanine green (ICG): an intravenous bolus or infusion of ICG, its distribution in the plasma (ICG is bound to the plasma proteins and does not enter the erythrocytes or the tissues), its hepatic uptake by OATP1B3 (saturable, competitively inhibited by bilirubin), its reflux from the liver into the sinusoidal plasma (the slow secondary phase of the plasma curve), its export into the bile by MRP2 (saturable) and the transit of the bile into the feces. ICG is not metabolized and not excreted by the kidneys. The model describes the liver function test of the plasma disappearance rate (ICG-PDR), the retention after 15 min (ICG-R15), the clearance and the hepatic extraction ratio, and the effects of cirrhosis (the loss of functional tissue with the shunts around the liver), hepatectomy (the resection rate), the hepatic blood flow, the cardiac output, the body weight, the hematocrit, bilirubin and the activity of OATP1B3.

The published version of the model is [icg-model](https://github.com/matthiaskoenig/icg-model) ([doi:10.5281/zenodo.5552397](https://doi.org/10.5281/zenodo.5552397)).

## Structure

- `__init__.py`: the paths of the model and its model object `MODEL` with its release metadata and its adopted fit, which the commands below run (pkpdsim).
- `models/`: the liver model (`model_liver.py`) and the whole-body model (`model_body.py`). `models/factory.py` creates the SBML models, the flattened whole-body model which is simulated, the ODEs and the COMBINE archive in `models/output/`; the fitted values come from the adopted fit (`models/fits.py`, `models/output/fits/`).
- `experiments/`: the simulation experiments of the clinical studies (`studies/`, the hepatectomies in `studies/hepatectomy/`) and the scans of dose, body weight, cirrhosis, resection, hepatic blood flow and cardiac output, liver volume, bilirubin and OATP1B3, the mapping of the Child-Turcotte-Pugh classes on the degree of cirrhosis and the reference timecourse (`misc/`). The fit mappings are described with the vocabulary of `pkpdsim.metadata`.
- `fitting/`: the fit problem of the pharmacokinetics (`PK`), its parameters and the selection of the training, validation, outlier and excluded data.
- `sensitivity/`: the sensitivity analysis of the pharmacokinetics and the liver function test.
- `icg_pk.py`: the pharmacokinetic parameters of the scans, the liver function test (the log-linear regression of the plasma curve in the window of a study, 2 - 15 min by default) and the non-compartmental analysis of `pkpdutils`.
- `classification/`: the prediction of the survival after a hepatectomy from the post-operative ICG-R15 of the model (Seyama2009, Wakabayashi2004), see `classification/README.md`.
- `nhanes/`: the sampling of virtual individuals from the NHANES data.
- `tolerances.py`: the absolute tolerances of the integrator.
- `simulations.py`: the simulation experiments and their groups, `MODEL` runs them.
- `data/`: the data of the studies, a subset of the PK-DB studies (`data/collect_data.py`).
- `TODO.md`: the open issues of the model.

## Commands

```bash
# models, simulations
run_simulation pkdb_models.models.icg --action factory
run_simulation pkdb_models.models.icg --action list_experiments
run_simulation pkdb_models.models.icg --action simulate --experiments all
run_simulation pkdb_models.models.icg --action simulate --experiments all --uncertainty

# fits, adoption of a fit, reports of stored parameters, identifiability
run_fit pkdb_models.models.icg --subset=PK --runs=3 --cores=3 --method=DE
run_fit pkdb_models.models.icg --adopt PK_<fit id>
run_report pkdb_models.models.icg --subset=PK results/fit/PK_<fit id>/parameters.json
run_identifiability pkdb_models.models.icg --subset=PK results/fit/PK_<fit id>/parameters.json

# sensitivity analysis (local by default)
python -m pkdb_models.models.icg.sensitivity.sensitivity_analysis

# classification of the survival after a hepatectomy (after two runs of ClassificationData)
python -m pkdb_models.models.icg.classification.classification
```

The results are written into `results/` of the model.

## Adopted fits and uncertainty

The values of the fitted parameters come from the adopted fit of the model,
not from numbers in the code. `run_fit pkdb_models.models.icg --adopt <fit id>`
copies the record of a finished fit (its parameters, problem, profile
likelihood and Fisher information) into `models/output/fits/<fit id>/`, which is
tracked; the ids of the adopted fits are listed in `fits` of the model object
(`__init__.py`). The factory takes the fitted values and their uncertainty from
the records (`models/fits.py`) and writes the uncertainty into the SBML
(`distrib:uncertainty`). `run_simulation ... --uncertainty` runs the
experiments over draws of these distributions and draws the simulations as the
median with a band of quantiles into `results/simulation_uncertainty`; the
classification and the sampled populations (Ohwada2006, Thomas2015,
Sakka2006) run with the fitted values only. A refit is: fit, adopt,
replace the id in `fits`, rebuild the model. See the documentation of the
repository on fitting (`docs/fitting.md`).

The adopted fit `PK_20261020_010006__506c2` fits the hepatic uptake (Vmax), the
reflux, the Vmax and Km of the export by MRP2 and the transit of the bile to the
plasma curves and extraction ratios of healthy subjects and to the plasma and
the biliary excretion of Meijer1988 (three doses). At 0.5 mg/kg the model gives
an ICG-PDR of 18 %/min, a half-life of 3.8 min, an ICG-R15 of 7 % and a
clearance of 560 - 640 ml/min, below the hepatic plasma flow of 790 ml/min.
