"""Simulation experiment of the study Thomas2015."""

from typing import ClassVar, override

import numpy as np
from matplotlib.axes import Axes
from pkpdsim.scan import start_values
from sbmlsim import Q
from sbmlsim.data import Data, DataSet
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Dimension, Scan

from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment
from pkdb_models.models.icg.experiments.fig_templates import base_figure
from pkdb_models.models.icg.icg_pk import PK_SELECTIONS, summarize
from pkdb_models.models.icg.nhanes.sampling import sampled_changes


class Thomas2015(ICGSimulationExperiment):
    """Simulation experiment of Thomas2015.

    Plasma disappearance rate of indocyanine green (0.25 mg/kg, LiMON) in 20
    patients with liver tumours, all measured during the operation: after the
    induction of the anaesthesia (t0), under the arterial and portovenous trial
    clamping of the segments to be resected (t1) and after the resection (t2)
    (Fig 3a: t0 vs t2, Fig 3b: t1 vs t2, Fig 3c: remnant volume vs PDR t2/t0);
    with the residual PDR on the first postoperative day vs the residual liver
    volume of Stockmann2009 (0.5 mg/kg, Fig 3B).

    The trial clamping takes the clamped segments out of the perfusion while
    the portal and the arterial inflow go on to the remnant, the situation after
    the resection (the premise of the paper): the clamping is simulated as the
    resection of the clamped part (`resection_rate`), so the model predicts the
    same PDR under the clamping and after the resection. The variation of the
    PDR after the resection is the sampled change of the blood flow after a
    hepatectomy.
    """

    #: the dose of 0.25 mg/kg at the 75 kg of the model [mg]
    dose: ClassVar[float] = 0.25 * 75
    #: the resection rates of the scans of the sampled blood flow
    resection_rates_long: ClassVar[np.ndarray] = np.linspace(0, 0.9, num=400)

    #: transparency of the sampled simulations by stage of cirrhosis
    alpha_resection: ClassVar[dict[str, float]] = {
        "Control": 0.7,
        "Mild cirrhosis": 0.7,
        "Moderate cirrhosis": 0.3,
        "Severe cirrhosis": 0.3,
    }
    alpha_clamping: ClassVar[dict[str, float]] = {
        "Control": 0.7,
        "Mild cirrhosis": 0.7,
        "Moderate cirrhosis": 0.2,
        "Severe cirrhosis": 0.1,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df_thomas = self.load_dataframe("Fig3")
        for label, key in [
            ("Fig3c_liver_remnant", "Fig3C_ICG_Thomas2015"),
            ("Fig3b_icgpdrt1", "Fig3B_ICG_Thomas2015"),
            ("Fig3a_preopicgpdr", "Fig3A_ICG_Thomas2015"),
        ]:
            df_subset = df_thomas.copy()[df_thomas.x_label == label]
            dsets[self.to_sid(key)] = DataSet.from_df(df_subset, self.ureg)

        dframes_stockmann = self.load_dataframes_by_substance(
            "Fig3", study="Stockmann2009"
        )
        for df in dframes_stockmann.values():
            dsets["Fig3_ICG_Stockmann2009"] = DataSet.from_df(df, self.ureg)

        return dsets

    @override
    def data(self) -> dict[str, Data]:
        # the scans record only what the pharmacokinetic parameters and the
        # figures (the sampled resection rates) read
        self.add_selections_data(selections=[*PK_SELECTIONS, "resection_rate"])
        return {}

    @override
    def simulations(self) -> dict[str, Scan]:
        tcscans = {}
        tcs = self.cirrhosis_timecourses(dose=self.dose)
        for key in self.cirrhosis_map:
            tc = tcs[key]

            # hepatectomy scan, also the trial clamping
            tcscans[self.to_sid(f"scan_rr_{key}")] = Scan(
                simulation=tc,
                dimensions=[
                    Dimension(
                        "dim_rr",
                        values={
                            "resection_rate": Q(self.resection_rates, "dimensionless")
                        },
                    ),
                ],
            )

            # bloodflow variation on top of resection
            tcscans[self.to_sid(f"sampled_rr_bf_{key}")] = Scan(
                simulation=tc,
                dimensions=[
                    Dimension(
                        "dim_rr_bf",
                        values={
                            **sampled_changes(
                                ["f_bloodflow", "resection_rate"],
                                size=400,
                                random_state=42,
                                hepatectomy=True,
                            ),
                            # the samples run with the fitted values, not over
                            # the draws
                            **self.fitted_values(400),
                        },
                    ),
                ],
            )

            # hepatectomy scan with the sampled postoperative blood flow, for
            # the clamping vs the resection
            tcscans[self.to_sid(f"scan_rr_long_{key}")] = Scan(
                simulation=tc,
                dimensions=[
                    Dimension(
                        "dim_rr",
                        values={
                            "resection_rate": Q(
                                self.resection_rates_long, "dimensionless"
                            ),
                            **sampled_changes(
                                ["f_bloodflow"],
                                size=len(self.resection_rates_long),
                                random_state=42,
                                hepatectomy=True,
                            ),
                            # the samples run with the fitted values, not over
                            # the draws
                            **self.fitted_values(len(self.resection_rates_long)),
                        },
                    ),
                ],
            )

        return tcscans

    @override
    def figures_mpl(self) -> dict[str, FigureMPL]:
        self.pk_dfs = self.calculate_icg_pk()

        return {
            **self.fig_icgpk_vs_remnant_vli(),
            **self.fig_pdrvli_vs_pdrqh(),
            **self.fig_icgpk_vs_remnant_vli_bfvar(),
            **self.fig_pdrvli_vs_pdrqh_bfvar(),
            **self.fig_pre_vs_post(),
        }

    def _pdr(self, sim_key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get the plasma disappearance rate of every point of a scan [%/min].

        Returns:
            The rate (the median of the draws in a run over them), and the
            lower and the upper quantile (the rate without draws).
        """
        df = self.pk_dfs[self.to_sid(sim_key)]
        unit = df["icgpdr_unit"].iloc[0]
        run = self.uncertainty
        if run is None:
            pdr = Q(df["icgpdr"].to_numpy(dtype=float), unit).to("percent/min")
            return pdr.magnitude, pdr.magnitude, pdr.magnitude
        summary = summarize(df, "icgpdr", run.quantiles)
        median, lower, upper = (
            Q(summary[column].to_numpy(), unit).to("percent/min").magnitude
            for column in ["median", "lower", "upper"]
        )
        return median, lower, upper

    def _sampled_resection_rates(self, key: str) -> np.ndarray:
        """Get the resection rates of the scan with the sampled blood flow."""
        xres = self.results[self.to_sid(f"task_sampled_rr_bf_{key}")]
        return (
            start_values(xres, "resection_rate")
            .transpose("dim_rr_bf", ...)
            .to_numpy()
            .reshape(xres.ds.sizes["dim_rr_bf"], -1)[:, 0]
        )

    def _legend(self, ax: Axes) -> None:
        """Draw the legend of a figure right of the axes."""
        ax.legend(
            fontsize=self.legend_font_size, loc="upper left", bbox_to_anchor=(1.02, 1.0)
        )

    def fig_pre_vs_post(self) -> dict[str, FigureMPL]:
        """Postoperative vs preoperative ICG-PDR."""
        f, ax = base_figure()
        ax.set_xlabel("Preoperative ICG-PDR [%/min]", fontdict=self.font)
        ax.set_ylabel("Postoperative ICG-PDR [%/min]", fontdict=self.font)

        for key in self.cirrhosis_map:
            pdr, _, _ = self._pdr(f"scan_rr_{key}")
            ax.plot(
                np.full_like(pdr, pdr[0]),
                pdr,
                marker="o",
                color=self.cirrhosis_colors[key],
                label=f"Sim: {key}",
            )

        # add data
        x = self._datasets["Fig3A_ICG_Thomas2015"].get_quantity("x")
        y = self._datasets["Fig3A_ICG_Thomas2015"].get_quantity("y")
        ax.plot(
            x.magnitude,
            y.magnitude,
            marker="s",
            linestyle="",
            color="dimgrey",
            label="Thomas 2015 (intraoperative)",
        )
        self._legend(ax)

        return {"fig_Thomas_pre_vs_post": f}

    def _plot_remnant_data(self, ax: Axes, white: bool) -> None:
        """Plot the post/preoperative ICG-PDR of Thomas2015 and Stockmann2009.

        Args:
            ax: the axes of the figure.
            white: white markers with a black edge on top of the simulations.
        """
        one = Q(1, "dimensionless")
        dset_thomas = self._datasets["Fig3C_ICG_Thomas2015"]
        x_thomas = one - dset_thomas.get_quantity("x").to("dimensionless")
        y_thomas = dset_thomas.get_quantity("y").to("dimensionless")
        # Stockmann2009: x the residual PDR, y the residual liver volume [%]
        dset_stockmann = self._datasets["Fig3_ICG_Stockmann2009"]
        x_stockmann = one - dset_stockmann.get_quantity("y").to("dimensionless")
        y_stockmann = dset_stockmann.get_quantity("x").to("dimensionless")
        for x, y, marker, label in [
            (x_thomas, y_thomas, "D", "Thomas 2015 (after resection)"),
            (x_stockmann, y_stockmann, "s", "Stockmann 2009 (POD1)"),
        ]:
            if white:
                ax.plot(
                    x.magnitude,
                    y.magnitude,
                    marker=marker,
                    linestyle="",
                    color="white",
                    ms=8,
                    markeredgecolor="black",
                    label=label,
                    zorder=2,
                )
            else:
                ax.plot(
                    x.magnitude,
                    y.magnitude,
                    marker=marker,
                    linestyle="",
                    color="dimgrey",
                    label=label,
                )

    def fig_icgpk_vs_remnant_vli(self) -> dict[str, FigureMPL]:
        """Change in ICG-PDR after hepatectomy."""
        f, ax = base_figure()
        ax.set_xlabel("Resection rate [-]", fontdict=self.font)
        ax.set_ylabel("Post/preoperative ICG-PDR [-]", fontdict=self.font)

        for key in self.cirrhosis_map:
            pdr, lower, upper = self._pdr(f"scan_rr_{key}")
            if self.uncertainty is not None:
                ax.fill_between(
                    self.resection_rates,
                    lower / pdr[0],
                    upper / pdr[0],
                    color=self.cirrhosis_colors[key],
                    alpha=0.2,
                )
            ax.plot(
                self.resection_rates,
                pdr / pdr[0],
                marker="o",
                color=self.cirrhosis_colors[key],
                label=f"Sim: {key}",
            )

        self._plot_remnant_data(ax, white=False)
        self._legend(ax)

        return {"fig_Thomas_Stockmann": f}

    def fig_pdrvli_vs_pdrqh(self) -> dict[str, FigureMPL]:
        """The ICG-PDR under trial clamping vs after the resection, Fig 3b.

        The clamping is the resection of the clamped part, the model predicts
        the same PDR (the line of identity of every stage).
        """
        f, ax = base_figure()
        ax.set_xlabel("ICG-PDR trial clamping [%/min]", fontdict=self.font)
        ax.set_ylabel("Postoperative ICG-PDR [%/min]", fontdict=self.font)

        for key in self.cirrhosis_map:
            pdr, _, _ = self._pdr(f"scan_rr_{key}")
            ax.plot(
                pdr,
                pdr,
                marker="o",
                color=self.cirrhosis_colors[key],
                label=f"Sim: {key}",
            )

        # add data
        x_data = self._datasets["Fig3B_ICG_Thomas2015"].get_quantity("x")
        y_data = self._datasets["Fig3B_ICG_Thomas2015"].get_quantity("y")
        ax.plot(
            x_data.magnitude,
            y_data.magnitude,
            marker="s",
            linestyle="",
            color="dimgrey",
            label="Thomas 2015 (intraoperative)",
        )
        self._legend(ax)

        return {"fig3b": f}

    def fig_icgpk_vs_remnant_vli_bfvar(self) -> dict[str, FigureMPL]:
        """Change in ICG-PDR after hepatectomy with sampled postoperative blood flow."""
        f, ax = base_figure()
        ax.set_xlabel("Resection rate [-]", fontdict=self.font)
        ax.set_ylabel("Post/preoperative ICG-PDR [-]", fontdict=self.font)

        for key in self.cirrhosis_map:
            pdr_ref, _, _ = self._pdr(f"scan_rr_{key}")
            pdr, _, _ = self._pdr(f"sampled_rr_bf_{key}")
            ax.plot(
                self._sampled_resection_rates(key),
                pdr / pdr_ref[0],
                marker="o",
                linestyle="",
                color=self.cirrhosis_colors[key],
                label=f"Sim: {key}",
                alpha=self.alpha_resection[key],
            )

        self._plot_remnant_data(ax, white=True)
        self._legend(ax)

        return {"fig_Thomas_Stockmann_bfvar": f}

    def fig_pdrvli_vs_pdrqh_bfvar(self) -> dict[str, FigureMPL]:
        """The ICG-PDR under trial clamping vs after the resection, sampled.

        As `fig_pdrvli_vs_pdrqh` with a sampled blood flow after the resection:
        the PDR under the clamping is the PDR of the resection with the blood
        flow unchanged (interpolated in `scan_rr`).
        """
        f, ax = base_figure()
        ax.set_xlabel("ICG-PDR trial clamping [%/min]", fontdict=self.font)
        ax.set_ylabel("Postoperative ICG-PDR [%/min]", fontdict=self.font)

        for key in self.cirrhosis_map:
            pdr_clamping, _, _ = self._pdr(f"scan_rr_{key}")
            pdr_post, _, _ = self._pdr(f"scan_rr_long_{key}")
            x = np.interp(self.resection_rates_long, self.resection_rates, pdr_clamping)
            ax.plot(
                x,
                pdr_post,
                marker="o",
                linestyle="",
                color=self.cirrhosis_colors[key],
                label=f"Sim: {key}",
                alpha=self.alpha_clamping[key],
            )

        # add data
        x_data = self._datasets["Fig3B_ICG_Thomas2015"].get_quantity("x")
        y_data = self._datasets["Fig3B_ICG_Thomas2015"].get_quantity("y")
        ax.plot(
            x_data.magnitude,
            y_data.magnitude,
            marker="s",
            linestyle="",
            color="white",
            ms=8,
            markeredgecolor="black",
            label="Thomas 2015 (intraoperative)",
            zorder=2,
        )
        self._legend(ax)

        return {"fig3b_bfvar": f}
