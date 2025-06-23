
# This script executes the download and preprocessing of the data


# activate the conda environment
source activate sc-hit

# convert the notebook to a script
jupyter nbconvert --output-dir=nbconverted --to script *.ipynb

# execute the script
python nbconvert/0.download-data.py
