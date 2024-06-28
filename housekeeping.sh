#!/usr/bin/env bash

# Lint all the python files

black ethical_sir.py
black scratch-demo-1.py
black scratch.py

# Keep the conda environment file up to date

conda env export > environment.yaml
