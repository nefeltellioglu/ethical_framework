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
    #config_file = "config/config-2024-10-14_manuscript.json"
     config_file = "config/config-2024-10-28_limited_vaccine.json"
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

output_dir = f"out/CZ_test_II/{config_date_name}"
os.makedirs(output_dir, exist_ok=True)
plot_df = pd.read_csv(f"{output_dir}/ab-heatmap-data.csv")

# ====================================================================

print("""
=======================================================
WARNING WARNING WARNING WARNING WARNING WARNING WARNING

This script generates a figure with an additional dummy
plot so that we can get a colourbar for the heatmaps.
You will need to post-process this in a vector graphics
editor to remove the dummy plot.
=======================================================
""")
colourbar_title = 'Vaccination (%)'
colour_scheme = "viridis_r"
fig, ax = plt.subplots(1, 3, figsize=(15, 4), gridspec_kw={'wspace': 0.2})

vac_1_df = 100*plot_df.pivot(index="b", columns="a", values="vac_1")
vac_2_df = 100*plot_df.pivot(index="b", columns="a", values="vac_2")

heat_map_1 = sns.heatmap(vac_1_df, cmap=colour_scheme,
                         ax=ax[0], cbar=False,
                         linewidth=0.5,
                         yticklabels=vac_1_df.index.round(1),
                         xticklabels=vac_1_df.columns.round(1))
heat_map_1.set_title("Group 1 vaccination",
                     fontweight='bold',
                     loc='left')
heat_map_1.set_xlabel("Equity in Infection Burden Multiplier (a)")
heat_map_1.set_ylabel("Equity in Vaccination Burden Multiplier (b)")
ax[0].invert_yaxis()
for label in ax[0].get_yticklabels():
    label.set_rotation(0)
for label in ax[0].get_xticklabels():
    label.set_rotation(0)

heat_map_2 = sns.heatmap(vac_2_df, cmap=colour_scheme,
                         ax=ax[1],
                         linewidth=0.5, cbar=False,
                         yticklabels=vac_2_df.index.round(1),
                         xticklabels=vac_2_df.columns.round(1))
heat_map_2.set_title("Group 2 vaccination",
                     fontweight='bold',
                     loc='left')
heat_map_2.set_xlabel("Equity in Infection Burden Multiplier (a)")
heat_map_2.set_ylabel("Equity in Vaccination Burden Multiplier (b)")
ax[1].invert_yaxis()
for label in ax[1].get_yticklabels():
    label.set_rotation(0)
for label in ax[1].get_xticklabels():
    label.set_rotation(0)

sns.heatmap(vac_2_df,
            cmap=colour_scheme,
            ax=ax[2],
            cbar_kws={'label': colourbar_title},
            cbar=True)
# remove the axis text for ax[2]
ax[2].set_yticks([])
ax[2].set_xticks([])

fig.savefig(f"{output_dir}/ab-heatmap-group-vaccination.png", bbox_inches='tight', dpi=300)
fig.savefig(f"{output_dir}/ab-heatmap-group-vaccination.svg", bbox_inches='tight')
