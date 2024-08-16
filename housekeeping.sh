#!/usr/bin/env bash

# Lint all the python files

black ethical_sir.py
black scratch-demo-1.py
black scratch.py
black create-grid-database.py

black ethics/model.py
black ethics/optimisation.py

# Keep the conda environment file up to date

conda env export > environment.yaml
