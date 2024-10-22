#!/bin/bash

# activate conda env
conda activate CFReT-compound-prio

# convert notebooks into python scripts
jupyter nbconvert --to python --output-dir=nbconverted/ *.ipynb

# run the scripts
python nbconverted/0.download-data.py
