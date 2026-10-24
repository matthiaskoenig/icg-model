"""Processing of the NHANES demographics data.

See https://www.cdc.gov/nchs/nhanes/ for more information.

Use of the continuous NHANES files (1999-2018), downloaded 2021-07-21.
The important subset used for analysis consists of the demographic information,
and the body measurements (BMX and DEMO subsets).

The nhanes.tsv provides the merged data set reduced to the following fields:

    SEQN: 'Respondent sequence number',
    RIAGENDR: 'Gender of the sample person (male/female)',
    RIDAGEYR: 'Best age in years of the sample person at time of HH screening.
               Individuals 85 and over are topcoded at 85 years of age.',
    RIDAGEMN: 'Best age in months.',
    RIDRETH1: 'Recode of reported race and ethnicity information',
    SEQN: 'Respondent sequence number',
    BMXWT: 'Weight (kg)',
    BMIWT: 'Weight Comment',
    BMXHT: 'Standing Height (cm)',
    BMIHT: 'Standing Height Comment',
    BMXBMI: 'Body Mass Index (kg/m**2)',
    BMXWAIST: 'Waist Circumference (cm)',
    BMIWAIST: 'Waist Circumference Comment'


nhanes.tsv provides concatenated data from R script.

The processed data (adults which are not obese, without the top coded ages) is
written to `NHANES_PROCESSED_PATH`, the individuals of the sampling, see
`sampling.sample_nhanes_individuals`.
"""

import pandas as pd

from pkdb_models.models.icg import NHANES_INPUT_PATH, NHANES_PROCESSED_PATH


def process_nhanes() -> pd.DataFrame:
    """Process and filter the NHANES data and write it to `NHANES_PROCESSED_PATH`.

    Returns:
        The processed data.
    """
    nhanes: pd.DataFrame = pd.read_csv(
        NHANES_INPUT_PATH, sep="\t", dtype={"SEQN": int}, index_col="SEQN"
    )
    # Remove top codings
    #   A: topcoded at 85 years of age
    #   B: topcoded at 85 years of age
    #   C: topcoded at 85 years of age
    #   D: topcoded at 85 years of age
    #   E: topcoded at 80 years of age
    #   F: topcoded at 80 years of age
    #   G: topcoded at 80 years of age
    #   H: topcoded at 80 years of age
    #   I: topcoded at 80 years of age
    #   J: topcoded at 80 years of age

    # - remove adolescents and children (age < 18)
    # - remove obese subjects (bmi >30)

    nhanes2 = nhanes[((nhanes.BMXBMI <= 30) & (nhanes.RIDAGEYR >= 18))]
    condition_topcode = (
        nhanes2["nhanes_id"].isin(["A", "B", "C", "D"]) & (nhanes2.RIDAGEYR < 85)
    ) | (
        nhanes2["nhanes_id"].isin(["E", "F", "G", "H", "I", "J"])
        & (nhanes2.RIDAGEYR < 80)
    )
    nhanes2 = nhanes2[condition_topcode]

    # https://www.cdc.gov/obesity/adult/defining.html
    # If your BMI is less than 18.5, it falls within the underweight range.
    # If your BMI is 18.5 to <25, it falls within the healthy weight range.
    # If your BMI is 25.0 to <30, it falls within the overweight range.
    # If your BMI is 30.0 or higher, it falls within the obesity range.
    nhanes2["BMICLASS"] = ""
    nhanes2.loc[nhanes.BMXBMI < 18.5, "BMICLASS"] = "underweight"
    nhanes2.loc[(nhanes.BMXBMI >= 18.5) & (nhanes.BMXBMI < 25), "BMICLASS"] = (
        "healthy weight"
    )
    nhanes2.loc[(nhanes.BMXBMI >= 25) & (nhanes.BMXBMI < 30), "BMICLASS"] = "overweight"

    nhanes2["AGECLASS"] = ""
    nhanes2.loc[(nhanes.RIDAGEYR <= 40), "AGECLASS"] = "18-40"
    nhanes2.loc[(nhanes.RIDAGEYR > 40) & (nhanes.RIDAGEYR <= 65), "AGECLASS"] = "41-65"
    nhanes2.loc[(nhanes.RIDAGEYR > 65), "AGECLASS"] = "66-84"

    # complete information for the following fields must exist
    subset = [
        "RIAGENDR",
        "RIDAGEYR",
        "BMXWT",
        "BMXHT",
        "BMXBMI",
        "RIDRETH1",
        "AGECLASS",
        "BMICLASS",
    ]
    nhanes2 = nhanes2.dropna(subset=subset)
    nhanes2 = nhanes2[subset]
    nhanes2.to_csv(NHANES_PROCESSED_PATH, sep="\t")
    return nhanes2


def main() -> None:
    """Process the NHANES data and show the processed data."""
    process_nhanes()

    # load the data
    nhanes = pd.read_csv(NHANES_PROCESSED_PATH, sep="\t")
    with pd.option_context("display.max_columns", None):
        print(nhanes.head(10))
    print("Number of subjects:", len(nhanes))


if __name__ == "__main__":
    main()
