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

with open(CONFIG["database_file"], "rb") as f:
    db = pickle.load(f)
    assert len(db["model_parameters"]) == 1

# ==============================================================================

# plot_df.to_csv(f"{output_dir}/ab-heatmap-data.csv", index=False)
# read the plot_df from csv
plot_df = pd.read_csv(f"{output_dir}/ab-heatmap-data.csv")

figsize=(7,8)
labels = ('A', 'B', 'C', 'D', 'E', 'F')
j = 0
fig = plt.figure(figsize=figsize)
variables = ["cli_burden", "total_vacc_perc", "vac_1", "vac_2"]
titles = ["Total Clinical Burden",
          "Total vaccination (%)",
          "Vaccinations in group 1 (%)",
          "Vaccinations in group 2 (%)"]
colors = ["rocket_r", "viridis_r", "viridis_r", "viridis_r"]

for var, title, color in zip(variables, titles, colors):
    no = 221 + j
    ax = fig.add_subplot(no)
    ax.text(-0.3, 1.05, labels[j], transform=ax.transAxes,
      fontsize=12, fontweight='bold', va='top', ha='right')
    j += 1

    if var in ["cli_burden", "total_vacc_perc", "total_vacc",
               "inf_burden", "adv_burden"]:
        perc_var = var
    elif var in [ "inf_1", "inf_2","vac_1", "vac_2"]:
        perc_var = "%s_perc"%var
        plot_df[perc_var] = 100 * plot_df[var]

    else:
        perc_var = "%s_perc"%var
        plot_df[perc_var] = 100 * plot_df[var]/CONFIG["population_parameters"]["pop_size_%s"%(var.split("_")[1])]
    plot_df["a"] = [round(i,2) for i in plot_df["a"]]
    plot_df["b"] = [round(i,2) for i in plot_df["b"]]
    data = plot_df.pivot(index="b", columns="a", values = perc_var)

    ax = sns.heatmap(data, linewidth=0.5,
                     vmin=min(plot_df[perc_var]),
                     vmax=max(plot_df[perc_var]),
                     #yticklabels=['High','Medium','Low','No'],
                     cbar_kws={'label': title,"location":'bottom',
                               "pad":0.25},
                     cmap=color,annot=False, fmt='.0f')
    if var in [ "vac_1", "vac_2", "total_vacc_perc"]:
        plt.gca().collections[0].set_clim(0,100)
    cax = ax.figure.axes[-1]
    cax.tick_params(labelsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("Equity in Infection\nBurden Multiplier (a)")
    ax.set_ylabel("Equity in Vaccination\nBurden Multiplier (b)")
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.1)

fig.savefig(f"{output_dir}/ab-heatmap-vaccination-and-clinical-burden.png", bbox_inches='tight', dpi=300)
