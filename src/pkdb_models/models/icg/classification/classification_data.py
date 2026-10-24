"""Data of the classification of the survival after hepatectomy (Koeller2021).

The survival after hepatectomy of Seyama2009 and Wakabayashi2004 with the
resection rate and the pre-operative ICG-R15, extended by the degree of
cirrhosis (`f_cirrhosis`) and the post-operative ICG-R15 which the model
predicts. The dataset is the input of `classification.py`.
"""

import logging
from collections.abc import Callable
from typing import ClassVar, override

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from sbmlsim import Q
from sbmlsim.data import Data, DataSet
from sbmlsim.simulation import Dimension, Scan, Simulation
from scipy.optimize import curve_fit

from pkdb_models.models.icg import RESULTS_PATH_CLASSIFICATION
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.icg_pk import PK_SELECTIONS
from pkdb_models.models.icg.models.fits import FITS

logger = logging.getLogger(__name__)

#: the survival data with the degree of cirrhosis, written by every run
SURVIVAL_DATASET_PATH = RESULTS_PATH_CLASSIFICATION / "survival_dataset.tsv"

#: the survival data with the post-operative ICG-R15 of the model, written by a
#: run which finds the survival data of a previous run
CLASSIFICATION_DATASET_PATH = RESULTS_PATH_CLASSIFICATION / "classification_dataset.tsv"

# FIXME: works for now, fitted function divides values by 100 (i think because of
#  unit=percent), therefore ICGR15 values from seyama is multiplied by a hundred
#  before converted into f_cirrhosis values


def fit_func(x: np.ndarray, a: float, c: float) -> np.ndarray:
    """Logarithmic dependency of `f_cirrhosis` on ICG-R15 with offset.

    Negative values are set to 0.0, there is no negative degree of cirrhosis.
    A factor of ICG-R15 inside the logarithm, `a * log(b * x) + c`, is the
    offset `a * log(b)`, which `c` already is: with it the fit had no unique
    optimum and tried negative factors, whose logarithm is `nan`.

    Args:
        x: ICG-R15 [-].
        a: factor of the logarithm.
        c: offset.

    Returns:
        The degree of cirrhosis `f_cirrhosis` [-].

    Raises:
        ValueError: if `x` is not a numpy array.
    """
    if not isinstance(x, np.ndarray):
        raise ValueError(f"Interpolation data must be an np.array, but '{type(x)}'")

    y = a * np.log(x) + c
    if np.any(y < 0.0):
        logger.warning("Negative values set to 0.0")
        y[y < 0.0] = 0.0

    return y


class ClassificationData(ICGSimulationExperiment):
    """Preparation of the classification data.

    **studies**
    - Seyama2009
    - Wakabayashi2004

    Survivors/Nonsurvivors with resection rate and pre-operative ICG kel.
    Calculation of f_cirrhosis (estimation liver disease) using the model.
    Heatmap with survival status.

    The post-operative ICG-R15 of the model is simulated for the survival data
    of a previous run (`survival_dataset.tsv`), so the classification dataset
    is written by the second run.

    The classification is an analysis of the data with the fitted model: in a
    run over the draws of the fitted parameters (`--uncertainty`) its
    simulations set the fitted parameters to their fitted values, so they run
    once and not over the draws (about 17000 simulations, times 100 draws).
    """

    n_samples: ClassVar[int] = 130
    resection_rates: ClassVar[np.ndarray] = np.linspace(0, 0.9, num=n_samples)

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for sid in ["Seyama2009_Fig3", "Wakabayashi2004_Fig5"]:
            study, fig_id = sid.split("_", 1)
            df = self.load_dataframe(fig_id, study=study)
            dsets[sid] = DataSet.from_df(df, self.ureg)
        return dsets

    @override
    def data(self) -> dict[str, Data]:
        # the 130 scans of 130 simulations record only what the pharmacokinetic
        # parameters read: every selection of the base experiment takes 5 GB
        self.add_selections_data(selections=PK_SELECTIONS)
        return {}

    @override
    def simulations(self) -> dict[str, Scan]:
        tcscans = {}
        for k, res_rate in enumerate(self.resection_rates):
            # FIXME: better to have a 2D-scan and calculated parameters from them
            # The current iterated 1D scan is just a work-around
            tcscans[f"scan_Fvli_{k}"] = Scan(
                simulation=Simulation(
                    time_unit="min",
                    end=30,
                    preinit_changes={
                        **self.default_changes(),
                        **self.fitted_changes(),
                        "IVDOSE_icg": Q(37.5, "mg"),
                        "FVli": Q((1 - res_rate) * 0.021, "l/kg"),
                    },
                ),
                dimensions=[
                    Dimension(
                        "dim_cirrhosis",
                        values=self.cirrhosis_scan_changes(num=self.n_samples),
                    )
                ],
            )

        # simulation to add model-calculated ICG-R15 to classification dataset
        if not SURVIVAL_DATASET_PATH.exists():
            logger.warning("Survival data does not exist on first run, execute again!")
            return tcscans

        df = pd.read_csv(SURVIVAL_DATASET_PATH, sep="\t")
        for k, (f_cirrhosis, resection_rate) in enumerate(
            zip(df["f_cirrhosis"], df["resection_rate"], strict=True)
        ):
            # this is a scan with length 1, because pk-parameters can only be
            # calculated from Scan
            tcscans[f"postop_r15_model_{k}"] = Scan(
                simulation=Simulation(
                    time_unit="min",
                    end=30,
                    preinit_changes={
                        **self.default_changes(),
                        **self.fitted_changes(),
                        "IVDOSE_icg": Q(37.5, "mg"),
                        "FVli": (1 - resection_rate) * Q(0.021, "l/kg"),
                    },
                ),
                dimensions=[
                    Dimension(
                        "dim_cirrhosis",
                        values={
                            "f_tissue_loss": Q([f_cirrhosis], "dimensionless"),
                            "f_shunts": Q([f_cirrhosis], "dimensionless"),
                        },
                    )
                ],
            )

        return tcscans

    def fitted_changes(self) -> dict[str, float]:
        """Get the changes of the fitted parameters to their fitted values.

        Returns:
            In a run over the draws the fitted values of the drawn parameters
            (in the units of the model), which the draws then leave out;
            nothing otherwise.
        """
        run = self.uncertainty
        if run is None:
            return {}
        return {pid: FITS[pid].value for pid in run.dimension.values}

    def icgr15_as_f_cirrhosis(self) -> Callable[[np.ndarray], np.ndarray]:
        """Fit the dependency of f_cirrhosis on ICG-R15.

        The function `fit_func` is fitted to the ICG-R15 of the scan of the
        degree of cirrhosis without resection; the parameters are logged.

        Returns:
            The fitted function f_cirrhosis(ICG-R15 [-]).
        """
        # icgr15 vector [-]
        xdf = self.pk_dfs["scan_Fvli_0"]  # total liver (no resection)
        xdata = (
            Q(xdf["icgr15"].to_numpy(), xdf["icgr15_unit"][0])
            .to("dimensionless")
            .magnitude
        )

        # f_cirrhosis vector [-]
        ydata = np.linspace(start=0.0, stop=0.9, num=self.n_samples)

        # Use non-linear least squares to fit a function, f, to data.
        pars, _ = curve_fit(f=fit_func, xdata=xdata, ydata=ydata)
        logger.info("f_cirrhosis = a * log(ICG-R15) + c: (a, c) = %s", pars)

        def f_fit(x: np.ndarray) -> np.ndarray:
            """Fitted function f_cirrhosis(ICG-R15 [-])."""
            return fit_func(x, *pars)

        return f_fit

    def survival_data(self) -> pd.DataFrame:
        """Prepare the survival data.

        The identical dataset is used for plotting and classification and
        written to `survival_dataset.tsv`:

        - `study`: study identifier ([Seyama2009|Wakabayashi2004])
        - `nonsurvival`: outcome ([0|1])
        - `resection_rate`: resection rate ([0-1])
        - `preop_r15`: preoperative ICG-R15 ([0-1])
        - `f_cirrhosis`: f_cirrhosis ([0-1])
        - `postop_r15_data`: post-operative ICG-R15 ([0-1]) (calculated from
          resection rate & preop ICGR15 from data)

        Returns:
            The survival data.
        """
        f_prediction = self.icgr15_as_f_cirrhosis()
        df_dict: dict[str, list] = {
            "study": [],
            "nonsurvival": [],
            "resection_rate": [],
            "preop_r15": [],
            "f_cirrhosis": [],
        }

        for key, dset in self._datasets.items():
            df_dict["nonsurvival"].extend(dset["nonsurvival"])
            df_dict["study"].extend(dset["study"])
            if key == "Wakabayashi2004_Fig5":
                df_dict["preop_r15"].extend(dset["y"] / 100)
                df_dict["resection_rate"].extend(1 - dset["x"])
            else:
                df_dict["preop_r15"].extend(dset["x"] / 100)
                df_dict["resection_rate"].extend(dset["y"] / 100)

        df_dict["f_cirrhosis"].extend(f_prediction(np.array(df_dict["preop_r15"])))

        df = pd.DataFrame(df_dict)
        df["postop_r15_data"] = np.multiply(
            1 - np.array(df["resection_rate"]), df["preop_r15"]
        )

        RESULTS_PATH_CLASSIFICATION.mkdir(parents=True, exist_ok=True)
        df.to_csv(SURVIVAL_DATASET_PATH, sep="\t", index=False)
        return df

    def classification_data(self) -> pd.DataFrame | None:
        """Add the post-operative ICG-R15 of the model to the survival data.

        The survival data with the column `postop_r15_model` (post-operative
        ICG-R15 [0-1] of the model) is written to `classification_dataset.tsv`.

        Returns:
            The classification data, None on the first run without survival data.
        """
        if not SURVIVAL_DATASET_PATH.exists():
            return None

        df = self.survival_data()
        # the predicted postop ICG-R15 of every patient
        df["postop_r15_model"] = [
            pk_df["icgr15"].to_numpy()[0] / 100.0  # Hack for dimensionless
            for key, pk_df in self.pk_dfs.items()
            if key.startswith("postop_r15_model_")
        ]
        df.to_csv(CLASSIFICATION_DATASET_PATH, sep="\t", index=False)
        return df

    @override
    def figures_mpl(self) -> dict[str, Figure]:
        self.pk_dfs = self.calculate_icg_pk()
        self.classification_data()
        return self.fig_resection_vs_f_cirrhosis()

    def fig_resection_vs_f_cirrhosis(self) -> dict[str, Figure]:
        """Plot the heatmap of the post-operative ICG-R15 with the survival.

        The post-operative ICG-R15 of the model over the degree of cirrhosis
        and the resection rate, with the survivors and non-survivors of the
        data.

        Returns:
            The figures by id, `fig_icgr15`.
        """
        pk_key = "icgr15"
        cmap = ListedColormap(
            [
                "#fff7fb",
                "#ece7f2",
                "#d0d1e6",
                "#a6bddb",
                "#74a9cf",
                "#3690c0",
                "#0570b0",
                "#045a8d",
                "#023858",
                "#011d2e",
            ]
        )

        f, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.2, 7.2))
        ax.tick_params(axis="x", labelsize=self.tick_font_size)
        ax.tick_params(axis="y", labelsize=self.tick_font_size)

        # simulation heatmap, the resection rate in the rows and the degree of
        # cirrhosis in the columns
        matrix = np.array(
            [
                self.pk_dfs[f"scan_Fvli_{k}"][pk_key].to_numpy()
                for k in range(len(self.resection_rates))
            ]
        )
        heatmap = ax.imshow(
            matrix,
            cmap=cmap,
            extent=(0.0, 0.9, 0, 0.9),
            origin="lower",
            vmin=0,
            vmax=100,
        )

        # contours
        f_cirrhosis = np.linspace(0, 0.9, self.n_samples)
        contours = ax.contour(
            f_cirrhosis,
            self.resection_rates,
            matrix,
            colors="black",
            levels=8,
            linewidths=0.5,
        )

        # colorbar
        cbar = f.colorbar(heatmap, fraction=0.038)
        # the gradient as vector cells: matplotlib rasterizes it into an
        # embedded image, which some SVG viewers do not show; the edges of
        # the cells in their color close the seams between them
        if cbar.solids is not None:
            cbar.solids.set_rasterized(False)
            cbar.solids.set_edgecolor("face")
        cbar.set_label(
            f"Postoperative {self.pk_labels[pk_key]} (predicted) [%]",
            fontdict={"weight": "bold", "size": 10},
        )

        ax.set_xlabel("f_cirrhosis [-]", fontdict=self.font)
        ax.set_ylabel("Resection rate [-]", fontdict=self.font)

        # survival data from data frame
        grouped = self.survival_data().groupby("nonsurvival")
        survivor_df = grouped.get_group(0)
        non_survivor_df = grouped.get_group(1)

        ax.plot(
            survivor_df["f_cirrhosis"],
            survivor_df["resection_rate"],
            marker="s",
            color="white",
            linestyle="",
            label="Survivors",
            markeredgecolor="black",
        )
        ax.plot(
            non_survivor_df["f_cirrhosis"],
            non_survivor_df["resection_rate"],
            marker="s",
            color="black",
            linestyle="",
            label="Non-survivors",
            markeredgecolor="black",
        )

        ax.set_ylim(bottom=0, top=0.9)
        ax.set_xlim(left=0, right=0.9)
        ax.legend(fontsize=self.legend_font_size)
        ax.clabel(
            contours,
            fmt={level: f"{int(level)}%" for level in contours.levels},
            fontsize=14,
        )

        return {f"fig_{pk_key}": f}
