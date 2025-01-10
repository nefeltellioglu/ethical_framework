import json
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import sys
import os
import itertools
import ethics.model as em
import ethics.optimisation as eo
import numpy as np


if len(sys.argv) > 1:
    config_file = sys.argv[1]
else:
    config_file = "config/config-2024-10-14_manuscript.json"
    # config_file = "config/config-2024-10-28_limited_vaccine.json"
    # config_file = "config/config-2024-12-02_limited_low_R0.json"
assert os.path.exists(config_file)

# NOTE This assumes the configuration file is named with the format
# `config-YYYY-MM-DD-<some_name>.json`. The `config_date_name` is used
# as the name for the output directory.
config_date_name = os.path.basename(config_file).replace("config-", "").replace(".json", "")


print(70 * "=")
print("Ethics plotting example results")
print(f"Using config file: {config_file}")
print(70 * "=")
with open(config_file, "r") as f:
    CONFIG = json.load(f)

output_dir = f"out/{config_date_name}"
os.makedirs(output_dir, exist_ok=True)
plot_df = pd.read_csv(f"{output_dir}/ab-heatmap-data.csv")

# ====================================================================

plot_title = "Clinical burden with optimal vaccination"
colourbar_title = 'Days of hospitalisation'
colour_scheme = "viridis_r"

cb_df = plot_df.pivot(index="b", columns="a", values="cli_burden")

fig = plt.figure(figsize=(7, 8))
cb_heatmap = sns.heatmap(cb_df, cmap=colour_scheme,
                         cbar_kws={'label': colourbar_title,
                                   'location': 'bottom',
                                   'shrink': 0.5,},
                         linewidths=0.5,
                         square=True,
                         yticklabels=cb_df.index.values.round(2),
                                xticklabels=cb_df.columns.values.round(2))
for label in cb_heatmap.get_yticklabels():
    label.set_rotation(0)
cb_heatmap.set_xlabel("Equity in Infection Burden Multiplier (a)")
cb_heatmap.set_ylabel("Equity in Vaccination Burden Multiplier (b)")
# invert the y-axis so that the origin is in the bottom left
cb_heatmap.invert_yaxis()
cb_heatmap.set_title(plot_title, fontweight='bold', loc='left')

fig.savefig(f"{output_dir}/ab-heatmap-clinical-burden.png", bbox_inches='tight', dpi=300)
fig.savefig(f"{output_dir}/ab-heatmap-clinical-burden.svg", bbox_inches='tight')
