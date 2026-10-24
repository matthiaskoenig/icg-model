"""Reference timecourse of indocyanine green with its parameter sensitivity."""

from typing import NamedTuple, override

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pkpdsim.scan import point_band, point_timecourse
from sbmlsim import Q
from sbmlsim.data import Data
from sbmlsim.result import ScanResult
from sbmlsim.simulation import Scan, Simulation

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure
from pkdb_models.models.icg.icg_pk import PK_SELECTIONS

#: the label of the reference point of a local design (`sampling.local`), the
#: dimension `dim_sens` of a sensitivity scan (`with_sensitivity`)
REFERENCE = "reference"


def perturbed_points(xres: ScanResult) -> ScanResult:
    """Get the perturbed points of a sensitivity scan, without its reference.

    Raises:
        KeyError: if the scan has no reference point.
    """
    return ScanResult(xres.ds.drop_sel(dim_sens=REFERENCE))


def _convert(values: np.ndarray, unit: str, target: str) -> np.ndarray:
    """Convert values of a unit into a target unit."""
    return np.asarray(Q(np.asarray(values, dtype=float), unit).to(target).magnitude)


class Panel(NamedTuple):
    """A panel of the reference timecourse."""

    label: str
    unit: str
    keys: list[str]
    colors: list[str]


class Sensitivity(NamedTuple):
    """The range of the parameter sensitivity of a timecourse."""

    mean: np.ndarray
    sd: np.ndarray
    min: np.ndarray
    max: np.ndarray

    def fill(self, ax: Axes, t_vec: np.ndarray, color: str) -> None:
        """Fill the mean +- SD and the minimum to maximum of the sensitivity.

        Args:
            ax: the axes to plot in.
            t_vec: the time of the timecourse.
            color: the color of the timecourse.
        """
        ax.fill_between(
            x=t_vec,
            y1=self.mean - self.sd,
            y2=self.mean + self.sd,
            color=color,
            alpha=0.1,
        )
        ax.fill_between(x=t_vec, y1=self.min, y2=self.max, color=color, alpha=0.1)


class ReferenceTimecourse(ICGSimulationExperiment):
    """Reference timecourse of an individual.

    The timecourses of an intravenous bolus of 37.5 mg and of an infusion of
    1 mg/min, with the range of the parameter sensitivity, and the
    pharmacokinetic parameters of the reference, which are written to
    `results/reference/ReferenceTimecourse/referencevalues.tsv`. In a run over
    the draws of the fitted parameters (`--uncertainty`) the band of the draws
    takes the place of the parameter sensitivity (which is not simulated), and
    the parameters are the median and the quantiles of the draws
    (`referencevalues_uncertainty.tsv`).
    """

    @override
    def simulations(self) -> dict[str, Simulation | Scan]:
        simulations = self.simulation_icg()
        if self.uncertainty is not None:
            return simulations
        # the reference of the pharmacokinetic parameters has no sensitivity
        scans = self.with_sensitivity(
            {k: v for k, v in simulations.items() if k != "icg_iv_reference"}
        )
        return {**simulations, **scans}

    @override
    def data(self) -> dict[str, Data]:
        # the sensitivity scans record only what the figures and the
        # pharmacokinetic parameters read
        keys = [key for panel in self.panels().values() for key in panel.keys]
        self.add_selections_data(selections=[*PK_SELECTIONS, *keys])
        return {}

    def simulation_icg(self) -> dict[str, Simulation | Scan]:
        """Get the simulations of the intravenous bolus and the infusion.

        Returns:
            The simulations by key.
        """
        return {
            "icg_iv": self.normal_timecourse(end=600),
            "icg_inf": Simulation(
                time_unit="min",
                end=1500,
                preinit_changes={
                    **self.default_changes(),
                    "Ri_icg": Q(1.0, "mg/min"),
                },
            ),
            "icg_iv_reference": self.normal_timecourse(end=40),
            "icg_iv_short": self.normal_timecourse(end=40),
        }

    @override
    def figures_mpl(self) -> dict[str, Figure]:
        figures = self.figure_inset()
        self.write_reference_pk_parameters()
        return figures

    def panels(self) -> dict[str, Panel]:
        """Get the panels of the timecourses by plot.

        The keys of the panels are the selections the figures read.
        """
        return {
            "plot1": Panel(
                label=f"ICG (plasma) [{self.unit_icg}]",
                unit=self.unit_icg,
                keys=["[Cve_icg]"],
                colors=["black"],
            ),
            "plot2": Panel(
                label=f"ICG [{self.unit_icg}]",
                unit=self.unit_icg,
                keys=["[Car_icg]", "[Chv_icg]"],
                colors=["red", "blue"],
            ),
            "plot3": Panel(
                label="Extraction ratio [-]",
                unit=self.unit_er,
                keys=["ER_icg"],
                colors=["black"],
            ),
            "plot4": Panel(
                label=f"ICG (liver) [{self.unit_icg}]",
                unit=self.unit_icg,
                keys=["[LI__icg]"],
                colors=["black"],
            ),
            "plot5": Panel(
                label=f"Bile excretion rate [{self.unit_bile_rate}]",
                unit=self.unit_bile_rate,
                keys=["LI__ICGLI2BI"],
                colors=["black"],
            ),
        }

    def figure_inset(self) -> dict[str, Figure]:
        """Plot the timecourses with the range of the parameter sensitivity.

        The panels of the liver and the bile excretion have an inset.

        Returns:
            The figures by id, `<simulation>_<panel>`.
        """
        figures = {}

        panels = self.panels()
        labels = {"[Car_icg]": "Arterial ICG", "[Chv_icg]": "Hepatic venous ICG"}

        for sim_key in ["icg_iv", "icg_inf"]:
            xres = self._results[f"task_{sim_key}"]
            # the perturbed points of the sensitivity scan, as sbmlsim 0.8 had
            # them; without a sensitivity scan in a run over the draws
            xres_sens = (
                perturbed_points(self._results[f"task_{sim_key}_sensitivity"])
                if self.uncertainty is None
                else None
            )
            uinfo = xres.units

            for plot, panel in panels.items():
                f, ax = base_figure()
                # set xlims for all panels
                ax.set_xlim(left=-1, right=30)

                if plot == "plot3":
                    ax.set_ylim(bottom=0, top=1)

                # define insets for bile excretion and liver icg
                inset_ax: Axes | None = None
                if plot == "plot4":
                    inset_ax = inset_axes(
                        ax,
                        width=3,
                        height=1.7,
                        bbox_to_anchor=(0.1, 0.6),
                        bbox_transform=ax.transAxes,
                        loc=3,
                    )
                    ax.set_ylim(bottom=0, top=50)
                if plot == "plot5":
                    inset_ax = inset_axes(
                        ax,
                        width=3,
                        height=1.7,
                        bbox_to_anchor=(0.11, 0.6),
                        bbox_transform=ax.transAxes,
                        loc=3,
                    )
                    ax.set_xlabel("Time [min]", font=self.font)
                    ax.set_ylim(bottom=0, top=0.1)

                # squared subplots
                ax.tick_params(axis="x", labelsize=self.tick_font_size)
                ax.tick_params(axis="y", labelsize=self.tick_font_size)

                # ylabels only on iv plot because they are next to each other in thesis
                if sim_key == "icg_iv":
                    ax.set_ylabel(panel.label, font=self.font)

                for key, color in zip(panel.keys, panel.colors, strict=True):
                    run = self.uncertainty
                    if run is None:
                        time, values = point_timecourse(xres, key)
                    else:
                        time, lower, values, upper = point_band(
                            xres, key, across=run.id, quantiles=run.quantiles
                        )
                    # the time in the unit of the model [min]
                    t_vec = np.asarray(time, dtype=float)

                    key_vec = _convert(values, uinfo[key], panel.unit)
                    ax.plot(t_vec, key_vec, color=color, label=labels.get(key))

                    if xres_sens is None:
                        # the band of the draws of the fitted parameters
                        band = Sensitivity(
                            key_vec,
                            key_vec * 0,
                            _convert(lower, uinfo[key], panel.unit),
                            _convert(upper, uinfo[key], panel.unit),
                        )
                    else:
                        # mean, sd, min and max of the sensitivity simulations
                        # over the sensitivity dimension at every time of the
                        # reference, interpolated between their own steps
                        sens_values = (
                            xres_sens.interpolate(time)[key]
                            .transpose("time", "dim_sens")
                            .values
                        )
                        band = Sensitivity(
                            *(
                                _convert(v, uinfo[key], panel.unit)
                                for v in (
                                    sens_values.mean(axis=1),
                                    sens_values.std(axis=1),
                                    sens_values.min(axis=1),
                                    sens_values.max(axis=1),
                                )
                            )
                        )
                    band.fill(ax, t_vec=t_vec, color=color)

                    # plot data in insets
                    if inset_ax is not None:
                        inset_ax.plot(t_vec, key_vec, color=color)
                        band.fill(inset_ax, t_vec=t_vec, color=color)

                # legend for arterial and hv icg
                if plot == "plot2":
                    ax.legend(fontsize=self.legend_font_size)
                figures[f"{sim_key}_{plot}"] = f
        return figures

    def write_reference_pk_parameters(self) -> None:
        """Write the pharmacokinetic parameters of the reference.

        The parameters of the liver function test of the reference with the
        mean and SD of the parameter sensitivity are written to
        `results/reference/ReferenceTimecourse/referencevalues.tsv`; in a run
        over the draws the median and the quantiles of the draws to
        `referencevalues_uncertainty.tsv`.
        """
        parameters = ["thalf_kel", "icgpdr", "icgr15", "cl_kel"]
        directory = MODEL.results_dir / "reference" / "ReferenceTimecourse"
        directory.mkdir(parents=True, exist_ok=True)

        if self.uncertainty is not None:
            df_ref = self.calculate_icg_pk(scans=["icg_iv_reference"])[
                "icg_iv_reference"
            ]
            rows = []
            for pk_key in parameters:
                median, lower, upper = self.pk_band(
                    df_ref, pk_key, self.pk_units[pk_key]
                )
                rows.append(
                    {
                        "pk_key": pk_key,
                        "unit": self.pk_units[pk_key],
                        "median": median[0],
                        "lower": None if lower is None else lower[0],
                        "upper": None if upper is None else upper[0],
                    }
                )
            pd.DataFrame(rows).to_csv(
                directory / "referencevalues_uncertainty.tsv", sep="\t", index=False
            )
            return

        pk_dfs = self.calculate_icg_pk(
            scans=["icg_iv_reference", "icg_iv_short_sensitivity"]
        )
        df_ref = pk_dfs["icg_iv_reference"]
        # the perturbed points, without the reference of the local design, as
        # sbmlsim 0.8 had them: the rows are the points of the scan
        labels = self.results["task_icg_iv_short_sensitivity"].ds["dim_sens"].values
        df = pk_dfs["icg_iv_short_sensitivity"].set_axis(labels).drop(index=REFERENCE)
        d: dict[str, list[str | float]] = {
            "pk_key": [],
            "unit": [],
            "reference_value": [],
            "sensitivity_mean": [],
            "sensitivity_sd": [],
        }

        for pk_key in parameters:
            unit = self.pk_units[pk_key]
            source = str(df_ref[f"{pk_key}_unit"].iloc[0])
            pk_ref = Q(df_ref[pk_key].to_numpy(dtype=float), source).to(unit)
            pk_sens = Q(df[pk_key].to_numpy(dtype=float), source).to(unit)

            d["pk_key"].append(pk_key)
            d["unit"].append(unit)
            d["reference_value"].append(float(pk_ref.magnitude[0]))
            d["sensitivity_mean"].append(float(np.mean(pk_sens.magnitude)))
            d["sensitivity_sd"].append(float(np.std(pk_sens.magnitude)))

        pd.DataFrame(data=d).to_csv(directory / "referencevalues.tsv", sep="\t")
