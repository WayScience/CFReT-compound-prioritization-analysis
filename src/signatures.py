"""This module provides statistical tests to identify significant differences in morphology features
between two profiles (reference and experimental). It supports Welch’s t-test, Kolmogorov–Smirnov
test, and permutation test, using scipy and statsmodels.

The core function, get_signatures, compares the two profiles using a specified test and a list
of morphology features. It returns two lists of features: significant (on-morphology) and non-significant
- On-morphology signatures: significant features associated with the cellular state.
- Off-morphology signatures: non-significant features not associated with the cellular state
"""


import numpy as np
import polars as pl
from scipy.stats import ks_2samp, permutation_test
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.weightstats import ttest_ind


def __add_significance_label(
    pvals_df: pl.DataFrame, sig_threshold: float | None = 0.05
) -> pl.DataFrame:
    """Add significance labels to p-values based on a threshold.

    Adds a new column to the dataframe indicating whether each p-value is
    significant or not based on the provided threshold.

    Parameters
    ----------
    pvals_df : pl.DataFrame
        DataFrame containing p-values.
    sig_threshold : float, optional
        Significance threshold. Default is 0.05.

    Returns
    -------
    pl.DataFrame
        DataFrame with added significance labels.
    """
    # setting
    return pvals_df.with_columns(
        pl.when(pl.col("corrected_p_value") < sig_threshold)
        .then(True)
        .otherwise(False)
        .alias("is_significant")
    )


def p_val_correction(p_values: np.ndarray, method: str = "fdr_bh") -> np.ndarray:
    """
    Perform multiple testing correction on p-values.

    Parameters
    ----------
    p_values : np.ndarray
        Array of p-values to be corrected.
    method : str, optional
        Method for correction. Default is "fdr_bh".

    Returns
    -------
    np.ndarray
        Array of corrected p-values.
    """

    _, corrected_p_values, _, _ = multipletests(p_values, method=method)
    return corrected_p_values


def __split_morphology_features(pval_df: pl.DataFrame) -> tuple[list[str], list[str]]:
    """Split features into two groups based on their significance.

    This function separates features into two categories: those that are
    significant (based on the "is_significant" column) and those that are not.

    Parameters
    ----------
    pval_df : pl.DataFrame
        A DataFrame containing feature names, and it's p-values


    Returns
    -------
    tuple
        A tuple containing two lists:
        - The first list contains the names of significant features. (on-morphology signature)
        - The second list contains the names of non-significant features. (off-morphology signature)

    Raises
    ------
    TypeError
        If the input is not a polars DataFrame.
    """
    # type checking
    if not isinstance(pval_df, pl.DataFrame):
        raise TypeError("pval_df must be a DataFrame")

    # now separate the morphology features that are significant and non-significant
    on_morph_feats = pval_df.filter(pl.col("is_significant"))["features"].to_list()
    off_morph_feats = pval_df.filter(~pl.col("is_significant"))["features"].to_list()

    return on_morph_feats, off_morph_feats


def welchs_ttest(
    ref_profiles: pl.DataFrame,
    exp_profiles: pl.DataFrame,
    morph_feats: list[str],
    correction_method: str | None = "fdr_bh",
    sig_threshold: float | None = 0.05,
) -> pl.DataFrame:
    """
    Perform Welch's t-test for each feature in the provided profiles and return a DataFrame with p-values.

    Parameters:
    - ref_profile: polars.DataFrame
        Reference profile containing features to be tested.
    - exp_profile: polars.DataFrame
        Experimental profile containing features to be tested.
    - morph_feats: List[str]
        List of feature names to perform the statistical test on.
    - correction_method: Optional[str], default="fdr_bh"
        Method for multiple testing correction (e.g., "fdr_bh", "bonferroni").
    - sig_threshold: Optional[float], default=0.05
        Significance threshold for corrected p-values.

    Returns:
    - polars.DataFrame
        DataFrame containing features, raw p-values, corrected p-values, and significance labels.
    """
    # Dictionary to store p-values for each feature
    pvals = {}
    for morph_feat in morph_feats:
        try:
            # Perform Welch's t-test (two-sided, unequal variance)
            _, p_value, _ = ttest_ind(
                ref_profiles[morph_feat].to_numpy(),
                exp_profiles[morph_feat].to_numpy(),
                alternative="two-sided",
                usevar="unequal",
                value=0,
            )
        except ValueError as e:
            # Handle errors (e.g., insufficient data) and assign NaN for the feature
            print(f"Error in t-test for {morph_feat}: {e}")
            pvals[morph_feat] = np.nan
            continue

        # Store the computed p-value
        pvals[morph_feat] = p_value

    # Create a DataFrame to store features and their corresponding p-values
    pvals_df = pl.DataFrame(
        {
            "features": morph_feats,
            "pval": [pvals[morph_feat] for morph_feat in morph_feats],
        }
    )

    # Apply multiple testing correction and add corrected p-values
    return (
        pvals_df.with_columns(
            pl.Series(
                "corrected_p_value",
                p_val_correction(pvals_df["pval"].to_numpy(), method=correction_method),
            )
        )
        # Add significance labels based on corrected p-values
        .pipe(__add_significance_label, sig_threshold=sig_threshold)
    )


def permutation(
    ref_profiles: pl.DataFrame,
    exp_profiles: pl.DataFrame,
    morph_feats: list[str],
    n_resamples: int | None = 1000,
    correction_method: str | None = "fdr_bh",
    sig_threshold: float | None = 0.05,
    seed: int | None = 0,
):
    # setting internal statistical functions (since it's a required input for permutation_test)
    def diff_of_means(ref_vals, exp_vals):
        return np.mean(exp_vals) - np.mean(ref_vals)

    def diff_of_medians(ref_vals, exp_vals):
        return np.median(exp_vals) - np.median(ref_vals)

    # setting up dictionary to store functions
    pvals = {}

    # iterate through each feature and perform permutation test
    for morph_feat in morph_feats:
        # get the reference and experimental values
        ref_vals = ref_profiles[morph_feat].to_numpy()
        exp_vals = exp_profiles[morph_feat].to_numpy()

        # perform permutation test
        # if the permutation test fails, catch the exception and continue
        # sets pval to nan
        try:
            result = permutation_test(
                data=(ref_vals, exp_vals),
                statistic=diff_of_means,
                alternative="two-sided",
                n_resamples=n_resamples,
                random_state=seed,
            )
        except Exception as e:
            # handle the exception
            print(f"Error occurred for feature {morph_feat}: {e}")
            pvals[morph_feat] = np.nan
            continue

        # store p-value in dictionary
        pvals[morph_feat] = result.pvalue

    # convert p-values dictionary to a polars dataframe
    pvals_df = pl.DataFrame(
        {
            "features": pvals.keys(),
            "p_value": pvals.values(),
        }
    )

    # correct p-values using the specified method
    # correct p-values using Benjamini-Hochberg method
    return (
        pvals_df.with_columns(
            # calculate and add corrected p-values using the specified method
            pl.Series(
                "corrected_p_value",
                p_val_correction(
                    pvals_df["p_value"].to_numpy(), method=correction_method
                ),
            )
        )
        # from the output generated above:
        # Add significance label based on corrected p-values using the pipe() method
        .pipe(__add_significance_label, sig_threshold=sig_threshold)
    )


def ks_test(
    ref_profiles: pl.DataFrame,
    exp_profiles: pl.DataFrame,
    morph_feats: list[str],
    correction_method: str = "fdr_bh",
    sig_threshold: float | None = 0.05,
) -> pl.DataFrame:
    """Perform KS-test for each feature in the morphology profiles and identifies
    significant features.

    This function performs a Kolmogorov-Smirnov test for each feature in the morphology profiles
    and identifies significant features based on a specified p-value correction method and
    significance threshold. Returns a DataFrame containing feature names, p-values,
    corrected p-values,

    Parameters
    ----------
    ref_df : pl.DataFrame
        Reference DataFrame.
    exp_df : pl.DataFrame
        Experimental DataFrame.
    morph_feats : list[str]
        List of morphology feature names.
    correction_method : str, optional
        Method for p-value correction. Default is "fdr_bh".
    sig_threshold : float, optional
        Significance threshold for labeling features. Default is 0.05.

    Returns
    -------
    pl.DataFrame
        DataFrame containing feature names, p-values, corrected p-values, and significance labels.
    """

    # Perform KS-test for each column and directly create a DataFrame.
    # using a list comprehension to iterate over the morphology features
    # the list comprehension creates a list of dictionaries, each containing
    # the feature name and its corresponding p-value
    # Perform KS-test for each feature and store results in a list
    pvals = {}
    for morph_feat in morph_feats:
        # calculate the p-value using the KS-test
        # if the KS-test fails, catch the exception and continue
        # sets pval to nan
        try:
            p_value = ks_2samp(
                ref_profiles[morph_feat].to_numpy(),
                exp_profiles[morph_feat].to_numpy(),
                method="auto",
                nan_policy="omit",
            )[1]
        except Exception as e:
            # handle the exception
            print(f"Error occurred for feature {morph_feat}: {e}")
            pvals[morph_feat] = np.nan
            continue

        # store the p-value in the dictionary
        pvals[morph_feat] = p_value

    # Create a DataFrame from the results
    pvals_df = pl.DataFrame(
        {
            "features": morph_feats,
            "p_value": pvals.values(),
        }
    )

    # Apply p-value correction
    # adds a new column to the pvals_df with the corrected p-values
    # then adds a significance label based on the corrected p-values
    return (
        pvals_df.with_columns(
            # calculate and add corrected p-values using the specified method
            pl.Series(
                "corrected_p_value",
                p_val_correction(
                    pvals_df["p_value"].to_numpy(), method=correction_method
                ),
            )
        )
        # from the output generated above:
        # Add significance label based on corrected p-values using the pipe() method
        .pipe(__add_significance_label, sig_threshold=sig_threshold)
    )


def get_signatures(
    ref_profiles: pl.DataFrame,
    exp_profiles: pl.DataFrame,
    morph_feats: list[str],
    method: str | None = "ks_test",
    fdr_method: str | None = "fdr_bh",
    p_threshold: float | None = 0.05,
    permutation_resamples: int | None = 1000,
    seed: int | None = 0,
) -> tuple[list[str], list[str]]:
    """
    Identifies significant and non-significant features between two profiles.

    This function performs statistical tests to compare two profiles (reference and experimental)
    based on specified morphology features. It identifies significant features using the
    Kolmogorov-Smirnov (KS) test or other specified methods. The function applies p-value
    correction and labels features as significant or non-significant based on a given
    significance threshold.

    Parameters
    ----------
    ref_profile : pl.DataFrame
        Reference profile as a Polars DataFrame.
    exp_profile : pl.DataFrame
        Experimental profile as a Polars DataFrame.
    morph_feats : list[str]
        List of morphology feature names to compare.
    method : str, optional
        Statistical method to use for comparison. Default is "ks_test".
    p_threshold : Optional[float], optional
        Significance threshold for p-values. Default is 0.05.
    fdr_method : str, optional
        Method for p-value correction. Default is "fdr_bh".
    permutation_resamples : int, optional
        Number of resamples for permutation test. Default is 1000.

    Returns
    -------
    tuple
        A tuple containing two lists:
        - Significant features (on-morphology).
        - Non-significant features (off-morphology).

    Raises
    ------
    TypeError
        If input types are not as expected.
    ValueError
        If an invalid method is specified.
    """
    # set seed for reproducibility
    np.random.seed(seed)

    # checking if the method is valid
    supported_methods = ["ks_test", "permutation_test", "welchs_ttest"]
    if method not in supported_methods:
        raise ValueError(
            f"{method} is not a valid method. Supported methods are: {supported_methods}"
        )

    # type check
    if not isinstance(ref_profiles, pl.DataFrame):
        raise TypeError("ref_profile must be a DataFrame")
    if not isinstance(exp_profiles, pl.DataFrame):
        raise TypeError("exp_profile must be a DataFrame")
    if not isinstance(morph_feats, list):
        raise TypeError("feats must be a list")
    if not all(isinstance(feat, str) for feat in morph_feats):
        raise TypeError("All elements in feats must be strings")
    if not isinstance(p_threshold, (int, float)):
        raise TypeError("p_threshold must be an int or float")
    if not isinstance(fdr_method, str):
        raise TypeError("fdr_method must be a string")

    # selecting statistical test to determine the significance of the morphology features
    # and to create the on-morphology and off-morphology signatures
    if method == "ks_test":
        pvals_df = ks_test(
            ref_profiles=ref_profiles,
            exp_profiles=exp_profiles,
            morph_feats=morph_feats,
            correction_method=fdr_method,
            sig_threshold=p_threshold,
        )
    elif method == "permutation_test":
        pvals_df = permutation(
            ref_profiles=ref_profiles,
            exp_profiles=exp_profiles,
            morph_feats=morph_feats,
            n_resamples=permutation_resamples,
            correction_method=fdr_method,
            sig_threshold=p_threshold,
            seed=seed,
        )
    elif method == "welchs_ttest":
        pvals_df = welchs_ttest(
            ref_profiles=ref_profiles,
            exp_profiles=exp_profiles,
            morph_feats=morph_feats,
            correction_method=fdr_method,
            sig_threshold=p_threshold,
        )

    # Split the features into significant and non-significant based on the significance label
    # on_morphology = significant morphology features
    # off_morphology = non-significant morpholgoy features
    on_morphology_feats, off_morphology_feats = __split_morphology_features(
        pval_df=pvals_df
    )
    return on_morphology_feats, off_morphology_feats
