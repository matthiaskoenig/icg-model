"""Sampling of in silico individuals.

The distributions of the parameters of the model which vary between individuals
(the amount of OATP1B3, body weight, hepatic blood flow, resection rate) and
the samples drawn from them. The simulation experiments sample their changes
with `sampled_changes`; the individuals of the NHANES cohort (see `nhanes.py`)
are sampled with `sample_nhanes_individuals`. Run as a script, the parameters
of the distributions are stored as JSON in the NHANES data folder.
"""

import json
from typing import NamedTuple, Protocol, TypedDict

import numpy as np
import pandas as pd
from sbmlsim import Q
from sbmlsim.data import load_pkdb_dataframe
from sbmlsim.units import Quantity
from scipy.stats import lognorm, norm, uniform

from pkdb_models.models.icg import (
    DATA_PATH_GALACTOSE,
    DATA_PATH_ICG,
    NHANES_DIR,
    NHANES_PROCESSED_PATH,
    SAMPLES_INDIVIDUALS_NOPK_PATH,
)
from pkdb_models.models.icg.nhanes.fbloodflow import f_bloodflow


class FrozenDistribution(Protocol):
    """The methods of a frozen distribution of `scipy.stats` which are used."""

    def rvs(self, size: int = 1, random_state: int | None = None) -> np.ndarray:
        """Draw random samples of the distribution."""

    def pdf(self, x: np.ndarray) -> np.ndarray:
        """Get the probability density function at `x`."""

    def mean(self) -> float:
        """Get the mean of the distribution."""


class Distribution(NamedTuple):
    """Distribution of a parameter of the model."""

    distribution: FrozenDistribution
    unit: str
    parameters: dict[str, float]


class LiverVolumeBloodflow(TypedDict):
    """Joint normal distribution of liver volume and blood flow of an age group."""

    group: str
    cov: np.ndarray  # covariance matrix of volume and blood flow
    volbw_mean: float  # liver volume per body weight [ml/kg]
    volbw_sd: float
    bfbw_mean: float  # liver blood flow per body weight [ml/min/kg]
    bfbw_sd: float


def oatp1b3_distribution() -> Distribution:
    """Get the lognormal distribution of OATP1B3, scaled to a mean of 1.0.

    The distribution is fitted to the data of Peng2015.
    """
    df = load_pkdb_dataframe("Peng2015_Fig5", data_path=DATA_PATH_ICG)
    data = df.value

    # fit lognorm distribution to data with fixed loc=0.0
    s, loc, scale = lognorm.fit(data, floc=0)
    dist = lognorm(s, loc, scale)

    # new distribution, scaled to mean=1.0
    mean = dist.mean()
    return Distribution(
        lognorm(s=s, scale=scale / mean),
        "dimensionless",
        {"s": s, "scale": scale / mean, "mean": mean},
    )


def bodyweight_distribution(
    reference_bw: float = 75, variation: float = 0.5
) -> Distribution:
    """Get the uniform distribution of the body weight.

    Args:
        reference_bw: body weight in the center of the distribution [kg].
        variation: half width of the distribution relative to `reference_bw`.
    """
    loc = reference_bw * (1 - variation)
    scale = (reference_bw - loc) * 2

    return Distribution(
        uniform(loc=loc, scale=scale),
        "kg",
        {"loc": loc, "scale": scale},
    )


def bloodflow_distribution(
    reference_f: float = 1, variation: float = 0.5
) -> Distribution:
    """Get the uniform distribution of the hepatic blood flow (`f_bloodflow`).

    Args:
        reference_f: `f_bloodflow` in the center of the distribution.
        variation: half width of the distribution relative to `reference_f`.
    """
    loc = reference_f * (1 - variation)
    scale = (reference_f - loc) * 2

    return Distribution(
        uniform(loc=loc, scale=scale),
        "dimensionless",
        {"loc": loc, "scale": scale},
    )


def postop_bloodflow_distribution(
    reference_f: float = 1, variation: float = 0.6
) -> Distribution:
    """Get the uniform distribution of the change of blood flow after hepatectomy.

    The reference data is Kawasaki1991.

    Args:
        reference_f: `f_bloodflow` in the center of the distribution.
        variation: half width of the distribution relative to `reference_f`.
    """
    loc = reference_f * (1 - variation)
    scale = (reference_f - loc) * 2

    return Distribution(
        uniform(loc=loc, scale=scale),
        "dimensionless",
        {"loc": loc, "scale": scale},
    )


def resection_distribution() -> Distribution:
    """Get the uniform distribution of the resection rate between 0 and 0.9."""
    loc = 0
    scale = 0.9
    return Distribution(
        uniform(loc=loc, scale=scale),
        "dimensionless",
        {"loc": loc, "scale": scale},
    )


def cov_liver_volume_bloodflow() -> dict[str, LiverVolumeBloodflow]:
    """Get the distributions of liver volume and blood flow by age group.

    The distributions are fitted to the data of Wynne1989.

    Returns:
        The distribution of every age group by its name.
    """
    data = load_pkdb_dataframe("Wynne1989_TabA", data_path=DATA_PATH_GALACTOSE)

    agegroups = {
        "18-40": data[data.age <= 40],
        "41-65": data[(data.age > 40) & (data.age <= 65)],
        "66-84": data[data.age > 65],
    }

    results: dict[str, LiverVolumeBloodflow] = {}
    for group, df in agegroups.items():
        vol = df["livvolumekg"]
        bf = df["livbloodflowkg"]

        # covariance matrix
        cov = np.cov(m=vol, y=bf)

        # fit marginal distributions
        vol_mean, vol_sd = norm.fit(vol)
        bf_mean, bf_sd = norm.fit(bf)

        results[group] = {
            "group": group,
            "cov": cov,
            "volbw_mean": vol_mean,
            "volbw_sd": vol_sd,
            "bfbw_mean": bf_mean,
            "bfbw_sd": bf_sd,
        }

    return results


def sample_liver_volume_bloodflow(samples: pd.DataFrame) -> pd.DataFrame:
    """Add the liver volume and blood flow of the age group to the samples.

    Args:
        samples: individuals with their age group (`AGECLASS`), changed in place.

    Returns:
        The samples with the liver volume (`LIVVOLKG`) and blood flow
        (`LIVBFKG`) per body weight.
    """
    # read covariance information
    cov_info = cov_liver_volume_bloodflow()

    livvolkg = []
    livbfkg = []

    for _, row in samples.iterrows():
        group = row["AGECLASS"]
        dist_info = cov_info[group]

        vol_bf = np.random.multivariate_normal(
            (dist_info["volbw_mean"], dist_info["bfbw_mean"]),
            cov=dist_info["cov"],
            size=1,
        )
        livvolkg.append(vol_bf[0][0])
        livbfkg.append(vol_bf[0][1])

    samples["LIVVOLKG"] = livvolkg
    samples["LIVBFKG"] = livbfkg

    return samples


def sample_nhanes_individuals(n: int | None = 100000, random_seed: int = 42) -> None:
    """Sample individuals of the NHANES cohort and write them into a table.

    The NHANES data is processed in `nhanes.py`. The samples are written to
    `SAMPLES_INDIVIDUALS_NOPK_PATH`.

    Args:
        n: number of samples, the NHANES individuals are repeated as often as
            needed; without a number the individuals of the processed NHANES
            cohort are used once.
        random_seed: random seed of the sampling.
    """
    np.random.seed(random_seed)

    nhanes = pd.read_csv(NHANES_PROCESSED_PATH, sep="\t")
    n_nhanes = len(nhanes)

    # copy the nhanes data and use the first n entries
    if n is not None and n > 0:
        n_copy = int(np.ceil(n / n_nhanes))
        samples = pd.concat([nhanes for _ in range(n_copy)])
        samples = samples.head(n=n)
    else:
        samples = nhanes

    # oatp
    samples["FOATP1B3"] = oatp1b3_distribution().distribution.rvs(size=len(samples))

    # bsa via formula
    samples["BSA"] = (
        0.024265
        * np.power(samples["BMXWT"], 0.5378)
        * np.power(samples["BMXHT"], 0.3964)
    )

    # dose
    samples["IVDOSE_icg"] = 0.5 * samples["BMXWT"]

    # liver volume and hbf
    samples = sample_liver_volume_bloodflow(samples=samples)
    samples["f_bloodflow"] = f_bloodflow(samples["LIVBFKG"])

    SAMPLES_INDIVIDUALS_NOPK_PATH.parent.mkdir(parents=True, exist_ok=True)
    samples.to_csv(SAMPLES_INDIVIDUALS_NOPK_PATH, sep="\t", index=False)


def sampled_changes(
    pars: list[str], random_state: int, hepatectomy: bool = False, size: int = 200
) -> dict[str, Quantity]:
    """Sample the changes of parameters of the model.

    Args:
        pars: the parameters to sample, `LI__f_oatp1b3`, `BW`, `f_bloodflow`
            and `resection_rate`.
        random_state: random seed of the sampling.
        hepatectomy: sample `f_bloodflow` from the change of blood flow after a
            hepatectomy.
        size: number of samples.

    Returns:
        The samples of every parameter.
    """
    np.random.seed(random_state)
    distributions = {
        "LI__f_oatp1b3": oatp1b3_distribution(),
        "BW": bodyweight_distribution(),
        "f_bloodflow": bloodflow_distribution(),
        "resection_rate": resection_distribution(),
    }
    if hepatectomy:
        distributions["f_bloodflow"] = postop_bloodflow_distribution()
    changes: dict[str, Quantity] = {}
    for par in pars:
        distribution = distributions[par]
        changes[par] = Q(distribution.distribution.rvs(size=size), distribution.unit)
    return changes


def store_distributions() -> None:
    """Store the parameters of the distributions as JSON in the NHANES data folder.

    These are `oatp1b3_distribution` and `cov_liver_volume_bloodflow`.
    """
    with (NHANES_DIR / "data" / "oatp1b3_distribution.json").open("w") as f_json:
        json.dump(oatp1b3_distribution().parameters, fp=f_json, indent=2)

    with (NHANES_DIR / "data" / "cov_liver_volume_bloodflow.json").open("w") as f_json:
        pars_vol_bf = cov_liver_volume_bloodflow()
        print(pars_vol_bf)

        # the covariance matrices as lists, which are written as JSON
        json.dump(
            {
                group: {**results, "cov": results["cov"].tolist()}
                for group, results in pars_vol_bf.items()
            },
            fp=f_json,
            indent=2,
        )


if __name__ == "__main__":
    store_distributions()
