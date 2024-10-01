"""
Module: utils.py

A collection of common utility functions for data processing,
as well as for saving, loading, and writing files.
"""

import pathlib
from typing import Optional

import pandas as pd
from pycytominer.cyto_utils import infer_cp_features


def create_results_dir() -> pathlib.Path:
    """Creates a results directory in the current path

    Returns
    -------
    pathlib.Path
        Path to results directory
    """

    # create results directory. if it does exist, do not raise error
    results_dir = pathlib.Path("./results").resolve()
    results_dir.mkdir(exist_ok=True)

    return results_dir


def split_meta_and_features(
    profile: pd.DataFrame,
    compartments=["Nuclei", "Cells", "Cytoplasm"],
    metadata_tag: Optional[bool] = False,
) -> tuple[list[str], list[str]]:
    """Splits metadata and feature column names

    Parameters
    ----------
    profile : pd.DataFrame
        image-based profile
    compartments : list, optional
        compartments used to generated image-based profiles, by default
        ["Nuclei", "Cells", "Cytoplasm"]
    metadata_tag : Optional[bool], optional
        indicating if the profiles have metadata columns tagged with 'Metadata_'
        , by default False

    Returns
    -------
    tuple[list[str], list[str]]
        Tuple containing metadata and feature column names

    Note
    ----
    This function was found in:
    https://github.com/axiomcura/predicting-cell-injury-compounds/blob/main/src/utils.py#L508
    """

    # identify features names
    features_cols = infer_cp_features(profile, compartments=compartments)

    # iteratively search metadata features and retain order if the Metadata tag is not added
    if metadata_tag is False:
        meta_cols = [
            colname
            for colname in profile.columns.tolist()
            if colname not in features_cols
        ]
    else:
        meta_cols = infer_cp_features(profile, metadata=metadata_tag)

    return (meta_cols, features_cols)
