#!/usr/bin/env python

# ## Downloading Data
# In this notebook, we will download the dataset required for our compound prioritization analysis. The dataset will be stored in the `../data` directory, which is located one level above the current notebook directory. This `../data` folder will contain all the necessary data and metadata needed for the analysis.
#
# Specifically:
# - The metadata and feature column names will be saved in a JSON file for easy reference.
# - The profiles will be stored in Parquet format to ensure efficient storage and processing.
#
# This setup ensures that all analytical notebooks have consistent access to the relevant data files.
#
# **NOTE**: In this analysis, we will use well-level aggregate data from the cell-injury dataset. This will be updated once the CFReT HCS data becomes available. The goal is to ensure that changing the data profiles should not necessitate altering the entire analysis pipeline.

# In[1]:


import pathlib

import pandas as pd

# In[2]:


# creating a data folder
data_dir = pathlib.Path("../data").resolve()
data_dir.mkdir(exist_ok=True)


# In[3]:


# dowloading labeled cell injury data from:
# https://github.com/axiomcura/predicting-cell-injury-compounds/blob/main/results/0.feature_selection/cell_injury_profile_fs.csv.gz
cell_injury_df = pd.read_csv(
    "https://github.com/axiomcura/predicting-cell-injury-compounds/raw/refs/heads/main/results/0.feature_selection/cell_injury_profile_fs.csv.gz",
    compression="gzip",
    low_memory=False,
)
cell_injury_df


# In[4]:


# Saving cell-injury well data as a parquet file
cell_injury_df.to_parquet(data_dir / "labeled_cell_injury_df.parquet", index=False)
