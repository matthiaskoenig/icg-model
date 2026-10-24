"""Simulation experiment of the study Kanstrup1987."""

from typing import ClassVar, override

import numpy as np
from pkpdsim.scan import end_values
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Dimension, Scan, Simulation

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure


class Kanstrup1987(ICGSimulationExperiment):
    """Simulation experiment of Kanstrup1987.

    Change of the clearance of indocyanine green with the change of the
    hepatic blood flow (Fig 1, relative to the control period, Bradley method
    with a constant infusion of 0.2 - 0.3 µmol/min after a priming dose),
    simulated as a scan of the hepatic blood flow at the steady state of an
    infusion of 0.25 µmol/min.

    Group A are 7 healthy men (24 y, 74 kg) 8 - 32 min (period II) and
    140 - 160 min (period III) after ergotamine 0.5 mg i.v., group B 10
    patients (65 y, 54 kg) with a negative investigation for intestinal angina
    10 - 60 min after a standard meal; the 12 points of Rowell 1964 are the
    exercise of 7 healthy subjects (single injections). Only the exercise and
    the ergotamine period II change the blood flow alone; in period III (flow
    restored, clearance -27 %) and after food the uptake changes as well,
    which a scan of the blood flow does not describe.
    """

    #: the series of Fig 1: (group, intervention): label, marker, color
    series: ClassVar[dict[tuple[str, str], tuple[str, str, str]]] = {
        ("Rowell 1964", "exercise, icg_injection"): (
            "Rowell 1964 exercise",
            "o",
            "tab:blue",
        ),
        ("Group A", "ergotamin PII, icg_infusion"): (
            "Ergotamine 8 - 32 min (A)",
            "s",
            "black",
        ),
        ("Group A", "ergotamin PIII, icg_infusion"): (
            "Ergotamine 140 - 160 min (A)",
            "s",
            "white",
        ),
        ("Group B", "food, icg_infusion"): ("Food (B)", "^", "tab:orange"),
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        df = self.load_dataframe("Fig1")
        return {
            self.to_sid(f"Fig1_{group}_{intervention.split(',')[0]}"): DataSet.from_df(
                df[(df.group == group) & (df.intervention == intervention)],
                self.ureg,
            )
            for group, intervention in self.series
        }

    @override
    def simulations(self) -> dict[str, Simulation | Scan]:
        tc = Simulation(
            time_unit="min",
            end=100,
            preinit_changes={
                **self.default_changes(),
                # 0.25 µmol/min, the middle of 0.2 - 0.3 µmol/min
                "Ri_icg": (Q(0.25, "µmol/min") * self.Mr.icg).to("mg/min"),
            },
        )
        return {
            "scan_f_bloodflow": Scan(
                simulation=tc,
                dimensions=[
                    Dimension("dim_bloodflow", values=self.bloodflow_scan_changes()),
                ],
            ),
            "reference_f_bloodflow": tc,
        }

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        return {**self._figure_bloodflow_vs_clearance()}

    def _figure_bloodflow_vs_clearance(self) -> dict[str, FigureMPL]:
        """Change of the clearance with the change of the hepatic blood flow."""
        f, ax = base_figure()

        # simulation results
        xres = self.results["task_scan_f_bloodflow"]
        xres_ref = self.results["task_reference_f_bloodflow"]

        # the change of the steady state clearance of the infusion against the
        # reference, per draw in a run over the draws (the draws of the scan
        # and of the reference are aligned by the dimension of the draws)
        ratio = (
            100
            * end_values(xres, "CLinfusion_icg")
            / end_values(xres_ref, "CLinfusion_icg")
        )
        x = 100 * self.scan_values(xres, "f_bloodflow")
        run = self.uncertainty
        if run is None or run.id not in ratio.dims:
            values = np.asarray(ratio.values, dtype=float)
            lower = upper = None
        else:
            axis = ratio.dims.index(run.id)
            array = np.asarray(ratio.values, dtype=float)
            values = np.nanmedian(array, axis=axis)
            lower = np.nanquantile(array, run.quantiles[0], axis=axis)
            upper = np.nanquantile(array, run.quantiles[1], axis=axis)
        self.plot_band(
            ax, x, values, lower, upper, color="black", label="Blood flow scan"
        )
        ax.set_xlabel("Hepatic blood flow change\n[%]", fontdict=self.font)
        ax.set_ylabel(f"{self.pk_labels['cl']} change [%]", fontdict=self.font)

        # data
        for (group, intervention), (label, marker, fill) in self.series.items():
            dset = self._datasets[
                self.to_sid(f"Fig1_{group}_{intervention.split(',')[0]}")
            ]
            ax.plot(
                dset["hbf_percent"].to_numpy(),
                dset["icg_cl_percent"].to_numpy(),
                marker=marker,
                linestyle="",
                color="black" if fill == "white" else fill,
                markerfacecolor=fill,
                label=label,
            )
        ax.legend(fontsize=self.legend_font_size - 2, loc="upper left")

        ax.set_xticks(ticks=[50, 100, 150, 200])
        ax.set_xticklabels(["-50", "0", "+50", "+100"])
        ax.set_yticks(ticks=[20, 40, 60, 80, 100, 120, 140])
        ax.set_yticklabels(["-80", "-60", "-40", "-20", "0", "+20", "+40"])

        return {"fig_1": f}
