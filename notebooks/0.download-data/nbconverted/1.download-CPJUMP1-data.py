#!/usr/bin/env python

# # Download CPJUMP1 Data
#
# This notebook documents the workflow for downloading and processing data from the JUMP Cell Painting dataset, available in the Cell Painting Gallery.
#
# We focus on datasets where cells have been perturbed by overexpression of genes using open reading frame (ORF) vectors, which artificially increase the production of specific proteins.
#
# Key steps in this workflow:
# - Load configuration and metadata from a YAML file.
# - Filter experimental metadata to include only plates with ORF perturbations.
# - Download each plate's data as a CSV file, convert it to Parquet format, and save it in the `./data` directory.
#
# Each individual plate data file is saved as a `parquet` file in the `./data` folder. If a file with the same name already exists, it will be replaced with the newly downloaded data.

# In[1]:


import pathlib
import sys
import time

import polars
import tqdm

sys.path.append("../../")
from utils import io_utils

# In[2]:


# setting config path
config_path = pathlib.Path("../nb-configs.yaml").resolve(strict=True)

# setting results setting a data directory
data_dir = pathlib.Path("./data").resolve()
data_dir.mkdir(exist_ok=True)

# setting a path to save the experimental metadata
exp_metadata_path = (data_dir / "CPJUMP1-experimental-metadata.csv").resolve()


# In[3]:


# loading config file and setting experimental metadata URL
nb_configs = io_utils.load_configs(config_path)
CPJUMP1_exp_metadata_url = nb_configs["links"]["CPJUMP1-experimental-metadata-source"]

# read in the experimental metadata CSV file and only filter down to plays that
# have an ORF perturbation
exp_metadata = polars.read_csv(
    CPJUMP1_exp_metadata_url, separator="\t", has_header=True, encoding="utf-8"
)

# filtering the metadata to only includes plates that their perturbation types are orfs
exp_metadata = exp_metadata.filter(exp_metadata["Perturbation"].str.contains("orf"))

# save the experimental metadata as a csv file
exp_metadata.write_csv(exp_metadata_path)

# display
exp_metadata.head()


# In[ ]:


# setting CPJUMP1 source link, this points to the main directory where all the plate data
# is stored
header_link = nb_configs["links"]["CPJUMP1-source"]

# create a for loop with progress bar for downloading plate data
for plate in tqdm.tqdm(
    exp_metadata["Assay_Plate_Barcode"].to_list(), desc="Downloading plates"
):
    # constructing the plate data source URL
    plate_data_source = f"{header_link}/{plate}/{plate}_normalized_negcon.csv.gz"

    # reading the plate data from the source URL
    # if the plate cannot be downloaded and read, it will skip to the next plate
    try:
        orf_plate_df = polars.read_csv(plate_data_source, separator=",", has_header=True)
    except polars.errors.ReadError:
        print(f"Failed to download and read plate data for {plate}. Skipping...")
        continue

    # saving the plate data to a parquet file
    orf_plate_df.write_parquet(data_dir / f"{plate}_normalized_negcon.parquet")

    # sleep to avoid overwhelming the AWS hosting the data
    time.sleep(1)
