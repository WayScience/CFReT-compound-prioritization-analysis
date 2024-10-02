"""
module: signatures.py

Description
------------

This module contains functions to identify morphological features classified as "on" and "off"
morphological signatures.

On morphological signatures: indicate features that are significantly and distinctly
different from a reference (e.g., a negative control group). These features reflect a specific
change associated with the condition or treatment.

Off morphological signatures: represent features not associated with the given cellular
state and should remain unaffected when exposed to a compound. If these features change, they
may indicate off-target effects.

The module includes various statistical methods to determine the on/off morphological features.
The primary function, `get_morphology_signatures()`, returns a tuple of two lists representing the
"on" and "off" morphological signatures, respectively.
"""

from collections import defaultdict
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats


def weighted_ks_test(
    target: pd.DataFrame,
    reference: pd.DataFrame,
    p_tresh: Optional[float] = 0.05,
) -> dict:
    """Performs a weighted Kolmogorov-Smirnov (KS) test between the target and reference
    datasets for each morphological feature. This method is designed to handle
    unbalanced sample sizes by applying weights to the cumulative distribution functions (CDFs).

    The function calculates the KS statistic and p-value for each morphological feature,
    using Welch's t-test to compare the means of the weighted cumulative distributions
    between the two samples. It returns a dictionary storing the KS statistic,
    p-value, and whether the p-value falls below a specified threshold for each feature.

    Parameters
    ----------
    target : pd.DataFrame
        A DataFrame containing the morphological features of the target dataset (e.g.,
        treated cells). Each column represents a different feature, and each row represents
        a single observation (e.g., a cell).

    reference : pd.DataFrame
        A DataFrame containing the morphological features of the reference dataset (e.g.,
        control cells). Each column represents a different feature, and each row represents
        a single observation.

    p_tresh : Optional[float], default=0.05
        The significance threshold for the p-value. Features with a p-value below this threshold
        will be considered significant.

    Returns
    -------
    dict
        Contains the follow infromation:
        - 'feature_name': Name of the morphological feature
        - 'ks_stat': The KS statistic for the feature.
        - 'p_val': The p-value resulting from the Welch's t-test.
        - 'below_thresh': A boolean indicating whether the p-value is below the specified threshold.

    Notes
    -----
    - The test assumes that both the target and reference datasets have the same columns
        (features) in the same order.
    - The KS statistic measures the maximum distance between the empirical cumulative
        distribution functions of the two samples, while the Welch's t-test checks for
        significance by comparing the means of the weighted CDFs.
    - Weights are applied to each sample to account for the difference in sample sizes,
        ensuring that both contribute equally to the analysis.
    """

    # storing weighted ks stats and p values
    scores = defaultdict(lambda: None).fromkeys(
        ["feature_name", "ks_stat", "p_val", "below_thresh"]
    )

    for morphology_feature in target.columns.tolist():
        # step1: calculating weights for the target and reference samples ensures that each sample
        # contributes equally to the analysis, regardless of its size.
        target_weights = np.ones_like(target[morphology_feature]) / len(
            target[morphology_feature]
        )
        reference_weights = np.ones_like(reference[morphology_feature]) / len(
            reference[morphology_feature]
        )

        # Step 2: sort the values and weights
        sorted_reference_indices = np.argsort(reference[morphology_feature].to_numpy())
        sorted_target_indices = np.argsort(target[morphology_feature].to_numpy())

        sorted_reference_data = reference.iloc[sorted_reference_indices]
        sorted_target_data = target.iloc[sorted_target_indices]

        sorted_reference_weights = reference_weights[sorted_reference_indices]
        sorted_target_weights = target_weights[sorted_target_indices]

        # step 3: calculate weighted cumulative distribution functions (CDF)
        weighted_reference_cdf = np.cumsum(sorted_reference_weights) / np.sum(
            sorted_reference_indices
        )
        weighted_target_cdf = np.cumsum(sorted_target_weights) / np.sum(
            sorted_target_indices
        )

        # finding all the unique values
        all_values = np.unique(
            np.concatenate([sorted_reference_data, sorted_target_data])
        )

        # interpolate to get CDF values at these unique points
        target_cdf_at_values = np.interp(
            all_values, sorted_reference_data, weighted_reference_cdf, left=0, right=1
        )
        reference_cdf_at_values = np.interp(
            all_values, sorted_target_data, weighted_target_cdf, left=0, right=1
        )

        # calculating the maximum difference between the two weighted CDF's (KS stat)
        ks_stat = np.max(np.abs(target_cdf_at_values - reference_cdf_at_values))

        # using welch's t-test assuming no equal variances between groups, good for comparing
        # the means of unbalanced samples
        t_stat, p_val = stats.ttest_ind(
            reference_cdf_at_values, target_cdf_at_values, equal_var=False
        )

        # storing scores
        scores["feature_name"] = morphology_feature
        scores["ks_stat"] = ks_stat
        scores["p_val"] = p_val
        scores["below_thresh"] = p_val < p_tresh

    return scores


def get_morphology_signatures(
    target: pd.DataFrame,
    reference: pd.DataFrame,
    method: Optional[str] = "weighted_ks",
    p_cutoff: Optional[float] = 0.05,
) -> tuple[list[str], list[str]]:
    """A wrapper function to identify on and off morphological features by applying different
    statistical methods to compare the target and reference datasets.

    Parameters
    ----------
    target : pd.DataFrame
        A DataFrame containing the morphological features of the target dataset (e.g., treated cells).
        Each column represents a different feature, and each row represents a single observation (e.g., a cell).

    reference : pd.DataFrame
        A DataFrame containing the morphological features of the reference dataset (e.g., control cells).
        Each column represents a different feature, and each row represents a single observation.

    method : Optional[str], default="weighted_ks"
        The statistical method to use for selecting the on and off morphological features.
        Currently supported methods:
        - "weighted_ks": Uses a weighted Kolmogorov-Smirnov (KS) test to identify features.

    p_cutoff : Optional[float], default=0.05
        The significance threshold for the p-value. Features with p-values below this cutoff will be
        considered significantly different (on or off).
    """
    # checking if the method trying to use exists
    supported_methods = ["weighted_ks"]
    if method.lower() not in supported_methods:
        raise ValueError(f"{method} is not a supported method.")

    # checking if the feature space are identical
    feature_space_check = set(target.columns.tolist()) - set(reference.columns.tolist())
    if len(feature_space_check) != 0:
        raise ValueError("Feature spaces are not identical")

    # checking which method was seleted
    if method == "weighted_ks":
        result = weighted_ks_test(target=target, reference=reference, p_tresh=p_cutoff)
        return result
