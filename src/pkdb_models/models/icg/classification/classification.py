"""Classification of the survival after hepatectomy (Koeller2021).

Classifiers for the prediction of the survival after hepatectomy based on
indocyanine green data. The following models are fitted

- data based models: survival ~ pre-operative ICG-R15 (Data1A), survival ~
  post-operative ICG-R15 of the data (Data1B) and survival ~ (pre-operative
  ICG-R15, resection rate) (Data2)
- PBPK based models: survival ~ post-operative ICG-R15 of the model (PBPK1)
  and survival ~ (f_cirrhosis, resection rate) (PBPK2)

The classification dataset is written by the experiment `ClassificationData`;
the figures and the results of the cross validation are written to
`RESULTS_PATH_CLASSIFICATION`.
"""

import logging
from collections import Counter
from pprint import pprint
from typing import TypedDict

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from sklearn.metrics import (
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    make_scorer,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import ShuffleSplit, cross_validate, train_test_split
from sklearn.svm import SVC

from pkdb_models.models.icg import RESULTS_PATH_CLASSIFICATION
from pkdb_models.models.icg.classification.classification_data import (
    CLASSIFICATION_DATASET_PATH,
)
from pkdb_models.models.icg.experiments.base_experiment import ICGSimulationExperiment

logger = logging.getLogger(__name__)

RANDOM_SEED = 42

#: the classification models with a single feature, the post-operative or
#: pre-operative ICG-R15
MODELS_1D = ("Data1A", "Data1B", "PBPK1")


class ClassificationModel(TypedDict):
    """Classification model: the features of the classifier and the plot color."""

    features: list[str]
    color: str


def data_prep(df: pd.DataFrame) -> pd.DataFrame:
    """Prepare the classification dataset, the outcome `nonsurvival` as bool.

    Args:
        df: classification dataset, changed in place.

    Returns:
        The classification dataset.
    """
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    print("-" * 80, "\nData preparation\n", "-" * 80)

    df["nonsurvival"] = df["nonsurvival"].astype(bool)
    print(df.head())
    return df


def y_data(df: pd.DataFrame) -> np.ndarray:
    """Outcome (non-survival) of the patients."""
    return np.array(df["nonsurvival"])


def X_data(df: pd.DataFrame, features: list[str]) -> np.ndarray:
    """Features of the patients, one row per patient."""
    data = [df[feature] for feature in features]
    return np.array(data).T


def outcome_ids(df: pd.DataFrame) -> tuple[list[int], list[int]]:
    """Indices of the survivors and of the non-survivors.

    Args:
        df: classification dataset.

    Returns:
        The indices of the survivors and the indices of the non-survivors.
    """
    y = y_data(df)
    survivor_ids = [k for k, outcome in enumerate(y) if not outcome]
    nonsurvivor_ids = [k for k, outcome in enumerate(y) if outcome]
    return survivor_ids, nonsurvivor_ids


def histogram(df: pd.DataFrame, features: list[str], features_inset: list[str]) -> None:
    """Histogram of the post-operative ICG-R15 of survivors and non-survivors.

    The inset shows the histogram of the pre-operative ICG-R15.

    Args:
        df: classification dataset.
        features: feature of the histogram.
        features_inset: feature of the histogram in the inset.
    """
    survivor_ids, nonsurvivor_ids = outcome_ids(df)
    X = X_data(df=df, features=features)
    X_inset = X_data(df=df, features=features_inset)
    bins = (np.linspace(0, 1, 21) * 100).tolist()

    f_hist, ax_hist = create_plot(
        xlabel="Postoperative ICG-R15 (predicted) [%]", ylabel="Count", size=(7.2, 6)
    )

    # create inset
    inset_fontdict = {"weight": "bold", "size": 11}
    inset_hist = inset_axes(
        ax_hist,
        width=2.5,
        height=1.8,
        bbox_to_anchor=(0.53, 0.45, 1, 1),
        bbox_transform=ax_hist.transAxes,
        loc="lower left",
    )
    inset_hist.set_xlabel("Preoperative ICG-R15 [%]", fontdict=inset_fontdict)
    inset_hist.set_ylabel("Count", fontdict=inset_fontdict)

    for ax, x in zip([ax_hist, inset_hist], [X, X_inset], strict=True):
        ax.hist(
            x[survivor_ids].T.flatten() * 100,
            bins=bins,
            label="Survivors",
            color="white",
            alpha=0.7,
            edgecolor="black",
        )
        ax.hist(
            x[nonsurvivor_ids].T.flatten() * 100,
            bins=bins,
            label="Non-survivors",
            color="#045a8d",
            alpha=0.7,
            edgecolor="black",
        )
    ax_hist.legend(fontsize=ICGSimulationExperiment.legend_font_size)
    label_plot(ax_hist, "PBPK1")
    label_plot(inset_hist, "Data1A", fontsize=12, pos="upper right")

    plt.show()
    f_hist.savefig(RESULTS_PATH_CLASSIFICATION / "histogram_PBPK1_Data1A.svg")


def metrics_vs_cutoff(
    cv: ShuffleSplit, X: np.ndarray, y: np.ndarray
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Evaluation metrics of the classification by a cutoff of a single feature.

    The patients with the feature at or above the cutoff are classified as
    non-survivors. The metrics are evaluated on the test data of the cross
    validation splits.

    Args:
        cv: cross validation splits.
        X: feature of the patients, one column.
        y: outcome of the patients.

    Returns:
        The cutoffs, the means of the metrics over the splits and the means
        plus and minus the standard deviations.
    """
    cutoffs = np.linspace(0.01, 1, 100)
    means = []
    sds = []
    for cutoff in cutoffs:
        cutoff_df = pd.DataFrame(
            [
                evaluation_metrics(y_true=y[test], y_pred=X[test][:, 0] >= cutoff)
                for _, test in cv.split(X, y)
            ]
        )
        means.append(cutoff_df.mean())
        sds.append(cutoff_df.std())
    metric_means = pd.DataFrame(means)
    metric_sds = pd.DataFrame(sds)

    return cutoffs, metric_means, metric_means + metric_sds, metric_means - metric_sds


def decision_boundary(classifier: SVC, ax: Axes, stepsize: float = 0.002) -> None:
    """Plot the decision boundary of a 2D classification.

    Args:
        classifier: fitted classifier with two features in [0, 1].
        ax: axes to plot on.
        stepsize: step size of the grid of the features.
    """
    x_min, x_max, y_min, y_max = 0, 1, 0, 1
    xx, yy = np.meshgrid(
        np.arange(x_min, x_max, stepsize), np.arange(y_min, y_max, stepsize)
    )

    cm_bright = ListedColormap(["white", "#045a8d"])
    Z = classifier.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    ax.contourf(xx, yy, Z, cmap=cm_bright, alpha=1.0)
    ax.set_xlim(xx.min(), 1.0)
    ax.set_ylim(yy.min(), 1.0)


def get_cutoff(classifier: SVC) -> float:
    """Cutoff of a 1D classification.

    The classifier is evaluated on 1000 values in [0, 1].

    Args:
        classifier: fitted classifier with a single feature.

    Returns:
        The value before the first value (after the first one) which is
        classified as non-survival.

    Raises:
        ValueError: if no value is classified as non-survival.
    """
    cutoff_samples = np.linspace(0, 1, 1000)
    y_pred = classifier.predict(cutoff_samples.reshape(-1, 1))
    for k, cutoff in enumerate(cutoff_samples[:-1]):
        if y_pred[k + 1]:
            return float(cutoff)
    raise ValueError("No value in [0, 1] is classified as non-survival.")


def npv(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Negative predictive value (NPV), tn / (tn + fn).

    Args:
        y_true: true outcome.
        y_pred: predicted outcome.

    Returns:
        The negative predictive value, 1 without true negatives.
    """
    tn, _fp, fn, _tp = confusion_matrix(y_true=y_true, y_pred=y_pred).ravel()
    # FIXME: npv is not defined for tn, fn = 0
    if tn == 0:
        return 1  # fix for zero division
    return tn / (tn + fn)


def evaluation_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_score: np.ndarray | None = None,
    ndigits: int = 3,
) -> dict[str, float]:
    """Evaluation metrics.

    -> precision, recall, sensitivity, specificity, F1-score (balanced methods,
       e.g. balanced accuracy, f1-score instead of f-score)
    -> matthews_corrcoef
    see:
    - https://scikit-learn.org/stable/modules/model_evaluation.html#
    - https://en.wikipedia.org/wiki/Confusion_matrixClassification
    - https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6413-7

    Args:
        y_true: true outcome.
        y_pred: predicted outcome.
        y_score: score of the outcome, the ROC AUC is calculated if given.
        ndigits: digits the metrics are rounded to.

    Returns:
        The metrics by name.
    """
    metrics = {
        # unbalanced
        "precision": precision_score(y_true=y_true, y_pred=y_pred, zero_division=0),
        "recall": recall_score(y_true=y_true, y_pred=y_pred),
        # balanced
        "balanced_accuracy": balanced_accuracy_score(y_true=y_true, y_pred=y_pred),
        "f1": f1_score(y_true=y_true, y_pred=y_pred),
        "npv": npv(y_true=y_true, y_pred=y_pred),
    }
    if y_score is not None:
        metrics["roc_auc"] = roc_auc_score(y_true=y_true, y_score=y_score)
    return {k: round(v, ndigits=ndigits) for k, v in metrics.items()}


def mean_roc_curve(
    classifier: SVC, cv: ShuffleSplit, X: np.ndarray, y: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Mean ROC curve over the cross validation splits.

    The classifier is fitted to the train data of every split and the ROC
    curve of the test data is interpolated.

    Args:
        classifier: classifier, fitted on every split.
        cv: cross validation splits.
        X: features of the patients.
        y: outcome of the patients.

    Returns:
        The false positive rates, the mean true positive rates and the means
        plus and minus the standard deviations.
    """
    ip_fpr = np.linspace(start=0, stop=1, num=100)
    tprs = []

    # reusing the existing cross-validation splits
    for k, (train, test) in enumerate(cv.split(X, y)):
        classifier.fit(X[train], y[train])
        roc_y_score = classifier.predict_proba(X[test])[:, 1]
        fpr, tpr, _thresholds = roc_curve(y_true=y[test], y_score=roc_y_score)
        counts = Counter(fpr)
        for x, count in counts.items():
            if count > 1 and x < 1:
                for addition in np.linspace(0.001, count * 0.001, count):
                    fpr[list(fpr).index(x)] = x + addition
        ip_tpr = np.interp(x=ip_fpr, xp=fpr, fp=tpr)
        if ip_tpr[-1] != 1.0:
            logger.warning(
                "Interpolated ROC Curve #%s does not end at tpr=1, but instead at y=%s",
                k,
                ip_tpr[-1],
            )
        tprs.append(ip_tpr)

    mean_tpr = np.mean(tprs, axis=0)
    std_tpr = np.std(tprs, axis=0)
    upper_tpr = np.add(mean_tpr, std_tpr)
    lower_tpr = np.subtract(mean_tpr, std_tpr)
    return ip_fpr, mean_tpr, upper_tpr, lower_tpr


def create_plot(
    xlabel: str = "xlabel",
    ylabel: str = "ylabel",
    size: tuple[float, float] = (6, 6),
    roc: bool = False,
) -> tuple[Figure, Axes]:
    """Create a figure with a single axes.

    Args:
        xlabel: label of the x-axis.
        ylabel: label of the y-axis.
        size: size of the figure.
        roc: ROC plot, with the labels of a ROC curve and the diagonal.

    Returns:
        The figure and its axes.
    """
    if roc:
        xlabel = "1 - Specificity [-]"
        ylabel = "Sensitivity [-]"
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=size)
    ax.tick_params(axis="x", labelsize=ICGSimulationExperiment.tick_font_size)
    ax.tick_params(axis="y", labelsize=ICGSimulationExperiment.tick_font_size)
    ax.set_xlabel(xlabel, fontdict=ICGSimulationExperiment.font)
    ax.set_ylabel(ylabel, fontdict=ICGSimulationExperiment.font)
    if roc:
        ax.plot([0, 1], [0, 1], linestyle="-", color="dimgrey", alpha=0.4)
    fig.tight_layout()
    return fig, ax


def label_plot(
    ax: Axes, model_description: str, fontsize: int = 14, pos: str = "upper left"
) -> None:
    """Label the plot with the name of the classification model.

    Args:
        ax: axes to label.
        model_description: name of the classification model.
        fontsize: font size of the label.
        pos: location of the label.
    """
    anchored_text = AnchoredText(
        s=model_description, loc=pos, prop={"size": fontsize, "weight": "bold"}
    )
    ax.add_artist(anchored_text)


def plot_decision_boundary(
    df: pd.DataFrame,
    X: np.ndarray,
    classifier: SVC,
    model_description: str,
    features: list[str],
) -> None:
    """Plot the decision boundary of a 2D classification with the patients.

    Args:
        df: classification dataset.
        X: features of the patients, the second one is the resection rate.
        classifier: fitted classifier.
        model_description: name of the classification model.
        features: features of the classification model.
    """
    labels = {"preop_r15": "Preoperative ICG-R15 [%]", "f_cirrhosis": "f_cirrhosis [-]"}
    survivor_ids, nonsurvivor_ids = outcome_ids(df=df)

    fig_db, ax_db = create_plot(xlabel=labels[features[0]], ylabel="Resection Rate [-]")
    decision_boundary(classifier=classifier, ax=ax_db)
    # Plot the survivors
    ax_db.scatter(
        X[survivor_ids][:, 0],
        X[survivor_ids][:, 1],
        color="white",
        edgecolors="black",
        label="Survivors",
        marker="s",
    )
    # Plot the nonsurvivors
    ax_db.scatter(
        X[nonsurvivor_ids][:, 0],
        X[nonsurvivor_ids][:, 1],
        color="black",
        edgecolors="black",
        label="Non-survivors",
        marker="s",
    )
    label_plot(ax_db, model_description)
    ax_db.legend(fontsize=ICGSimulationExperiment.legend_font_size)
    if model_description == "Data2":
        ax_db.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax_db.set_xticklabels(["0", "20", "40", "60", "80", "100"])
    fig_db.savefig(RESULTS_PATH_CLASSIFICATION / f"{model_description}_db.svg")


def plot_metrics_vs_cutoff(
    cv: ShuffleSplit, X: np.ndarray, y: np.ndarray, model_description: str
) -> None:
    """Plot the evaluation metrics of a 1D classification over the cutoff.

    Args:
        cv: cross validation splits.
        X: feature of the patients, one column.
        y: outcome of the patients.
        model_description: name of the classification model.
    """
    cutoff_labels = {
        "Data1A": "Preoperative ICG-R15 cutoff [%]",
        "Data1B": "Postoperative ICG-R15 cutoff [%]",
        "PBPK1": "Postoperative ICG-R15 cutoff [%]",
    }
    metric_labels = {
        "precision": "Precision (PPV)",
        "recall": "Recall",
        # balanced
        "balanced_accuracy": "Balanced accuracy",
        "f1": "f1-score",
        "npv": "NPV",
    }
    cutoffs, metric_means, metric_upperbounds, metric_lowerbounds = metrics_vs_cutoff(
        cv, X, y
    )

    fig_metrics, ax_metrics = create_plot(
        xlabel=cutoff_labels[model_description], ylabel="Evaluation metrics [-]"
    )
    for metric in metric_means.columns:
        ax_metrics.plot(
            cutoffs * 100, metric_means[metric], label=metric_labels[metric]
        )
        ax_metrics.fill_between(
            cutoffs * 100,
            metric_upperbounds[metric],
            metric_lowerbounds[metric],
            alpha=0.35,
        )
    ax_metrics.legend(fontsize=ICGSimulationExperiment.legend_font_size)
    label_plot(ax_metrics, model_description)
    fig_metrics.savefig(
        RESULTS_PATH_CLASSIFICATION / f"{model_description}_metrics.svg"
    )


def classification(
    df: pd.DataFrame, classification_models: dict[str, ClassificationModel]
) -> None:
    """Classification of the survival status.

    Every classification model is fitted to the train data and to all data,
    evaluated by cross validation and plotted. The results of the cross
    validation are written to `crossvalidation_results.tsv`.

    Args:
        df: classification dataset.
        classification_models: classification models by name.
    """
    classifier = SVC(
        kernel="poly",
        degree=2,
        probability=True,
        class_weight="balanced",
        random_state=RANDOM_SEED,
    )
    optimal_classifier = SVC(
        kernel="poly",
        degree=2,
        probability=True,
        class_weight="balanced",
        random_state=RANDOM_SEED,
    )
    test_size = 0.25
    comparison_table = {}
    fig_roc1, ax_roc1 = create_plot(roc=True)
    fig_roc2, ax_roc2 = create_plot(roc=True)

    for model_description, model in classification_models.items():
        color = model["color"]
        features = model["features"]
        print(model_description)
        print(f"non-survival ~ {features}")

        X = X_data(df, features=features)
        y = y_data(df)

        # test and train data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=RANDOM_SEED
        )

        # --- classification ---
        print(f"\nClassifier: {classifier}")
        # fit classification model
        classifier.fit(X=X_train, y=y_train)
        optimal_classifier.fit(X=X, y=y)
        # predict test data
        y_pred = classifier.predict(X_test)
        y_score = classifier.predict_proba(X_test)[:, 1]
        # predict all data
        optimal_y_pred = optimal_classifier.predict(X)
        optimal_y_score = optimal_classifier.predict_proba(X)[:, 1]
        # confusion matrix
        conf_matrix = confusion_matrix(y_true=y_test, y_pred=y_pred)
        print("\nConfusion matrix")
        print(conf_matrix)
        # evaluation metrics
        metrics = evaluation_metrics(y_true=y_test, y_pred=y_pred, y_score=y_score)
        optimal_metrics = evaluation_metrics(
            y_true=y, y_pred=optimal_y_pred, y_score=optimal_y_score
        )
        print("\nEvaluation metrics")
        pprint(metrics)
        print("\nOptimal Evaluation metrics")
        pprint(optimal_metrics)
        if model_description in MODELS_1D:
            cutoff = get_cutoff(classifier=classifier)
            print(f"\nCutoff is at: {cutoff}")
        print("-" * 80)

        # --- cross validation ---
        cv = ShuffleSplit(n_splits=200, test_size=test_size, random_state=RANDOM_SEED)
        scoring = {
            "roc_auc": make_scorer(roc_auc_score),
            "balanced_accuracy": make_scorer(balanced_accuracy_score),
            "f1": make_scorer(f1_score),
            "precision": make_scorer(precision_score, zero_division=0),
            "recall": make_scorer(recall_score),
            "npv": make_scorer(npv),
        }
        scores = cross_validate(classifier, X, y, cv=cv, scoring=scoring)

        cross_validation = {}
        for score in scoring:
            mean = round(float(np.mean(scores[f"test_{score}"])), ndigits=3)
            std = round(float(np.std(scores[f"test_{score}"])), ndigits=3)
            cross_validation[score] = rf"{optimal_metrics[score]} ({mean} $\pm$ {std})"
        comparison_table[model_description] = cross_validation

        # --- mean roc curves ---
        ip_fpr, mean_tpr, upper_tpr, lower_tpr = mean_roc_curve(
            classifier=classifier, cv=cv, X=X, y=y
        )
        # refit the classifier with train data
        classifier.fit(X_train, y_train)

        # --- optimal roc curves ---
        optimal_fpr, optimal_tpr, _optimal_thresholds = roc_curve(
            y_true=y, y_score=optimal_y_score
        )

        # ------------------ figures ------------------
        # supplemental roc figures
        if model_description in ["Data1A", "Data1B", "Data2"]:
            fig_roc_supp, ax_roc_supp = create_plot(roc=True)
            ax_roc_supp.plot(
                optimal_fpr,
                optimal_tpr,
                linestyle="--",
                color=color,
                label=f"{model_description} ROC (AUC={optimal_metrics['roc_auc']})",
            )
            ax_roc_supp.plot(
                ip_fpr, mean_tpr, color=color, label=r"CV ROC (Mean$\pm$SD)"
            )
            ax_roc_supp.fill_between(
                ip_fpr, upper_tpr, lower_tpr, color=color, alpha=0.3
            )
            label_plot(ax_roc_supp, model_description)
            ax_roc_supp.legend(
                fontsize=ICGSimulationExperiment.legend_font_size, loc="lower right"
            )
            fig_roc_supp.savefig(
                RESULTS_PATH_CLASSIFICATION / f"{model_description}_roc.svg"
            )

        # decision boundary plots
        if model_description in ["PBPK2", "Data2"]:
            plot_decision_boundary(
                df=df,
                X=X,
                classifier=classifier,
                model_description=model_description,
                features=features,
            )

        # metrics vs cutoff plots
        if model_description in MODELS_1D:
            plot_metrics_vs_cutoff(cv=cv, X=X, y=y, model_description=model_description)

        # roc plots
        ax_roc = ax_roc1 if model_description in MODELS_1D else ax_roc2
        ax_roc.plot(
            ip_fpr,
            mean_tpr,
            color=color,
            label=rf"{model_description} CV ROC (Mean$\pm$SD)",
        )
        ax_roc.fill_between(ip_fpr, upper_tpr, lower_tpr, color=color, alpha=0.3)
        ax_roc.plot(
            optimal_fpr,
            optimal_tpr,
            label=f"{model_description} ROC (AUC={optimal_metrics['roc_auc']})",
            linestyle="--",
            color=color,
        )
        ax_roc.legend(fontsize=ICGSimulationExperiment.legend_font_size)

    plt.show()
    for fig, name in zip([fig_roc1, fig_roc2], ["1D", "2D"], strict=True):
        fig.savefig(RESULTS_PATH_CLASSIFICATION / f"{name}_roc.svg")

    # show cross validation results
    print("-" * 80)
    print("Cross validation")
    print("-" * 80)
    comparison_df = pd.DataFrame(comparison_table)
    print(comparison_df)
    comparison_df.to_csv(
        RESULTS_PATH_CLASSIFICATION / "crossvalidation_results.tsv",
        sep="&",
        index=False,
    )


def main() -> None:
    """Run the classification of the survival on the classification dataset."""
    df = pd.read_csv(CLASSIFICATION_DATASET_PATH, sep="\t", index_col=False)
    data_prep(df=df)

    classification_models: dict[str, ClassificationModel] = {
        "Data1A": {"features": ["preop_r15"], "color": "tab:red"},
        "Data1B": {"features": ["postop_r15_data"], "color": "tab:green"},
        "PBPK1": {"features": ["postop_r15_model"], "color": "black"},
        "Data2": {"features": ["preop_r15", "resection_rate"], "color": "tab:orange"},
        "PBPK2": {"features": ["f_cirrhosis", "resection_rate"], "color": "black"},
    }

    classification(df=df, classification_models=classification_models)
    histogram(
        df=df,
        features=classification_models["PBPK1"]["features"],
        features_inset=classification_models["Data1A"]["features"],
    )


if __name__ == "__main__":
    main()
