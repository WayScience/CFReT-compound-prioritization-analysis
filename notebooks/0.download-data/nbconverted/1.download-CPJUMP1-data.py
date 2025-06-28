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

# Setting parameters for the notebook:
#
# - `pert_type (str)`: Perturbation of interest

# In[2]:


pert_type = "orf"


# Setting input and output paths

# In[3]:


# setting config path
config_path = pathlib.Path("../nb-configs.yaml").resolve(strict=True)

# setting results setting a data directory
data_dir = pathlib.Path("./data").resolve()
data_dir.mkdir(exist_ok=True)

# creating a metadata directory
metadata_dir = (data_dir / "metadata").resolve()
metadata_dir.mkdir(exist_ok=True)

# creating a platemaps directory in the metadata directory
platemap_dir = (metadata_dir / "platemaps").resolve()
platemap_dir.mkdir(exist_ok=True)

# creating a profiles directory in the metadata directory
profiles_dir = (data_dir / "profiles").resolve()
profiles_dir.mkdir(exist_ok=True)


# Downloading the experimental metadata file and saving it into the `metadata` directory

# In[4]:


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
exp_metadata.write_csv(metadata_dir / "CPJUMP1-experimental-metadata.tsv", separator="\t")

# display
exp_metadata.head()


# To efficiently organize and download the CPJUMP1 data, we generate a dictionary (`batch_and_time_delay_dict`) that maps each experimental timepoint (e.g., Day0, Day1, Week2, Week4, DL) to its corresponding batch identifier.
#
# This mapping enables us to:
# - Group and save data files according to their acquisition timepoints after initial treatment.
# - Treat each delayed timepoint as a distinct batch for downstream processing.

# In[5]:


# extracting the unique batches and time delays from the experimental metadata
batch_and_time_delay = exp_metadata[["Batch", "Time_delay"]].unique(maintain_order=True)

# creating a dictionary to store the batch and time delay information
# {"time of delay": "Batch"} is the format of the dictionary
batch_and_time_delay_dict = {}
for row in batch_and_time_delay.rows(named=True):

    # extracting Batch and Time_delay from the row
    batch = row["Batch"]
    time_delay = row["Time_delay"]

    # ignore CPJUMP1_DL batch
    # there is no negative controled wells for this batch
    if batch.endswith("CPJUMP1_DL"):
        continue
    else:
        batch_and_time_delay_dict[time_delay] = batch

print("Delayed timepoints and their corresponding batches:")
batch_and_time_delay_dict


# Downloading Aggregated Profiles for Selected `pert_type`
#
# For each time point, all plates associated with the selected perturbation type (`pert_type`) will be downloaded. After downloading, the data from all plates for a given time point will be concatenated into a single file. Each output file will be prefixed with the corresponding time point to clearly indicate its contents.

# In[6]:


# setting CPJUMP1 source link, this points to the main directory where all the plate data
# is stored
header_link = nb_configs["links"]["CPJUMP1-profiles-source"]

# iterating over each time point
for delayed_time_point, batch_name in batch_and_time_delay_dict.items():

    # filter the experiential metadata to only include the plates that correspond
    # to the current time point
    exp_metadata_filtered = exp_metadata.filter(
        exp_metadata["Batch"] == batch_name
    )
    if exp_metadata_filtered.is_empty():
        raise ValueError(
            f"No plates found for batch {batch_name} at time point {delayed_time_point}. "
            "Please check the experimental metadata."
        )

    # get all the plate names for the current time points
    plate_names = exp_metadata_filtered["Assay_Plate_Barcode"].to_list()

    # iterate over each plate name and download the data
    loaded_profiles_df = []
    for plate_name in tqdm.tqdm(plate_names, desc=f"Downloading plates for {delayed_time_point}"):
        # set the profile src end point
        profile_end_point = f"{batch_name}/{plate_name}/{plate_name}_normalized_feature_select_negcon_batch.csv.gz"

        # set the full URL for the profile
        profile_url = f"{header_link}/{profile_end_point}"

        # downloading profile data from the URL link
        # if download fails, raise an error
        try:
            loaded_profile = polars.read_csv(
                profile_url,
                separator=",",
                has_header=True
            )

            # adding a small delay to avoid overwhelming the server
            time.sleep(0.7)
        except Exception as e:
            raise ValueError(f"Error downloading {plate_name} for {delayed_time_point}: {e}")

        # store the plate data in a list
        loaded_profiles_df.append(loaded_profile)

    # concatenate all the loaded profiles into a single DataFrame and save
    loaded_profiles_df = polars.concat(loaded_profiles_df, how="vertical")

    # save the concatenated profiles to a parquet file
    output_file = (profiles_dir / f"{batch_name}_profiles.parquet").resolve()
    loaded_profiles_df.write_parquet(output_file)


# Downloading the platemaps associated with the selected `pert_type`

# In[7]:


# downloading the plate maps
pert_platemap_url = nb_configs["links"]["CPJUMP1-platemaps-source"]
pert_platemap_name = f"JUMP-Target-1_{pert_type}_platemap.txt"

# constructing the full URL for the plate map
platemap_url = f"{pert_platemap_url}/{pert_platemap_name}"

# downloading the plate map
platemap_df = polars.read_csv(platemap_url, separator="\t")

# save it to the platemap directory
platemap_df.write_csv(platemap_dir / pert_platemap_name)

# display the first few rows of the plate map
platemap_df.head()
