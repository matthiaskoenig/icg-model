"""Simulation experiment of the intravenous application of indocyanine green."""

from typing import ClassVar, override

from sbmlsim import Q
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.icg import MODEL
from pkdb_models.models.icg.experiments.base_experiment import (
    ICGSimulationExperiment,
)


class DoseDependencyExperiment(ICGSimulationExperiment):
    """Intravenous injection of indocyanine green over a day."""

    doses: ClassVar[list[float]] = [37.5]  # [mg]
    routes: ClassVar[list[str]] = ["IV"]

    @override
    def simulations(self) -> dict[str, Simulation]:
        return {
            self.to_sid(f"icg_{route}_{dose}"): Simulation(
                time_unit="min",
                end=24 * 60,
                preinit_changes={
                    **self.default_changes(),
                    f"{route}DOSE_icg": Q(dose, "mg"),
                },
            )
            for route in self.routes
            for dose in self.doses
        }

    @override
    def figures(self) -> dict[str, Figure]:
        figures = {}
        dose_colors = ["tab:blue", "tab:orange", "tab:green"]

        for route in self.routes:
            fig = Figure(
                experiment=self,
                sid=self.to_sid(f"Fig_application_{route}"),
                num_rows=3,
                num_cols=1,
                name=f"ICG's PK dose dependency ({route})",
                settings={"legend_scope": "row"},
            )
            plots = fig.create_plots(xaxis=Axis("time", unit="hr"), legend=True)
            plots[0].set_yaxis(label="ICG plasma", unit="mmol/l")
            plots[1].set_yaxis(label="ICG feces", unit="mmol")
            plots[2].set_yaxis(label=f"ICG dose {route}", unit="mg")
            # every row shows the same doses: one legend, of the first row
            plots[1].legend = False
            plots[2].legend = False

            for kdose, dose in enumerate(self.doses):
                for plot, yid in zip(
                    plots,
                    ["[Cve_icg]", "Afeces_icg", f"{route.upper()}DOSE_icg"],
                    strict=True,
                ):
                    plot.add_data(
                        task=self.to_sid(f"task_icg_{route}_{dose}"),
                        xid="time",
                        yid=yid,
                        label=f"{dose} mg {route}",
                        color=dose_colors[kdose],
                    )

            figures[self.to_sid(f"Fig_application_{route}")] = fig

        return figures


if __name__ == "__main__":
    MODEL.run_experiments(
        DoseDependencyExperiment, output_dir=DoseDependencyExperiment.__name__
    )
