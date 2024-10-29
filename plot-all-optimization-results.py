import json
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import ethics.model as em
import ethics.optimisation as eo
import seaborn as sns
from scipy.interpolate import griddata
import matplotlib.tri as tri
from matplotlib import cm
import sys
import os

if len(sys.argv) > 1:
    config_file = sys.argv[1]
else:
    # config_file = "config/config-2024-10-14_manuscript.json"
    config_file = "config/config-2024-10-28_limited_vaccine.json"
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


# ====================================================================
# We only want the simulation results that use an amount of vaccine
# less than the maximum amount specified in the configuration file so
# we filter the results in the database to only include those records.
# ====================================================================
if "vaccine_parameters" in CONFIG:
    max_vacc = CONFIG["vaccine_parameters"]['maximum_vacc_rollout']
else:
    max_vacc = float('inf')

valid_ics = [
    ic for ic in db["initial_conditions"]
    if ic["value"].total_number_vaccinated() <= max_vacc
]
valid_ic_ids = [ic["id"] for ic in valid_ics]
valid_configs = [
    c for c in db["configurations"]
    if c["initial_condition_id"] in valid_ic_ids
]
valid_config_ids = [c["id"] for c in valid_configs]
valid_outcomes = [
    o for o in db["outcomes"]
    if o["configuration_id"] in valid_config_ids
]

db = {
    "model_parameters": db["model_parameters"],
    "initial_conditions": valid_ics,
    "configurations": valid_configs,
    "outcomes": valid_outcomes,
    "burden_parameters": db["burden_parameters"]
}
# ====================================================================


# At the point where we need to make some plots!
model_param_id = 0
burden_param_id = 0
#ethical_a = 0.1
#ethical_b = 0.5

bp = [bp for bp in db["burden_parameters"] if bp["id"] == burden_param_id][0][
    "parameters"
]

configs = [
    c for c in db["configurations"] if c["model_parameters_id"] == model_param_id
]
config_ids = [c["id"] for c in configs]
ocs = [o for o in db["outcomes"] if o["configuration_id"] in config_ids]
step = CONFIG["grid_search_step"]["a_b_grid_step"]


plot_df = []

for ethical_a in np.arange(0.0, 1 , step):
    for ethical_b in np.arange(0.0, 1 - ethical_a , step):


        foo, bar = eo.optimal_initial_condition(
            ethical_a, ethical_b, model_param_id, burden_param_id, db, normalise=True
        )
        _optimal_ic = [ic for ic in db["initial_conditions"] if ic["id"] == foo]
        _optimal_config = [c for c in db["configurations"] if c["initial_condition_id"] == foo]
        _optimal_outcome = [
            o for o in db["outcomes"] if o["configuration_id"] == _optimal_config[0]["id"]
        ]
        best_vac_1 = _optimal_outcome[0]["outcome"].total_vac_1
        best_vac_2 = _optimal_outcome[0]["outcome"].total_vac_2
        best_inf_1 = _optimal_outcome[0]["outcome"].inf_1_no_vac
        best_inf_2 = _optimal_outcome[0]["outcome"].inf_2_no_vac


        plot_df.append(
            {   "a": ethical_a,
                "b": ethical_b,
                #"outcome_id": oc["id"],
                #"config_id": tmp_config_id,
                "vac_1": best_vac_1,
                "vac_2": best_vac_2,
                "inf_1": best_inf_1,
                "inf_2": best_inf_2,
                #"loss": _optimal_outcome,
            }
        )
plot_df = pd.DataFrame(plot_df)



# Use matplotlib to make a plot of the same data.
# Put a single large red dot at the best outcome.

plt.figure()
plt.scatter(plot_df["a"], plot_df["b"], c=plot_df["vac_1"])
cbar = plt.colorbar()
cbar.set_label("No vacc in group 1")
plt.xlabel("a")
plt.ylabel("b")
#plt.savefig(f"{output_dir}/vac_group_1_across_all.png", dpi=300)
#plt.savefig(f"{output_dir}/vac_group_1_across_all.svg", dpi=300)
plt.clf()



plt.figure()
plt.scatter(plot_df["a"], plot_df["b"],
            c=100 * plot_df["vac_1"]/CONFIG["population_parameters"]["pop_size_1"])
cbar = plt.colorbar()
cbar.set_label("Vaccinations in group 1 (%)")
plt.xlabel("a")
plt.ylabel("b")
#plt.savefig(f"{output_dir}/vac_group_1_across_all_perc.png", dpi=300)
#plt.savefig(f"{output_dir}/vac_group_1_across_all_perc.svg", dpi=300)

plt.clf()

plt.figure()
plt.scatter(plot_df["a"], plot_df["b"],
            c=100 * plot_df["vac_2"]/CONFIG["population_parameters"]["pop_size_2"])
cbar = plt.colorbar()
cbar.set_label("Vaccinations in group 2 (%)")
#plt.xlabel("a")
#plt.ylabel("b")
plt.xlabel("Equity in Vaccination Multiplier (a)")
plt.ylabel("Equity in Clinical Burden Multiplier (b)")
#plt.savefig(f"{output_dir}/vac_group_2_across_all_perc.png", dpi=300)
#plt.savefig(f"{output_dir}/vac_group_2_across_all_perc.svg", dpi=300)

plt.clf()

plt.figure()
plt.scatter(plot_df["a"], plot_df["b"],
            c=100 * plot_df["inf_2"]/CONFIG["population_parameters"]["pop_size_2"])
cbar = plt.colorbar()
cbar.set_label("Infections in group 2 (%)")
plt.xlabel("a")
plt.ylabel("b")
#plt.savefig(f"{output_dir}/inf_group_2_across_all_perc.png", dpi=300)
#plt.savefig(f"{output_dir}/inf_group_2_across_all_perc.svg", dpi=300)

plt.clf()

plt.figure()
plt.scatter(plot_df["a"], plot_df["b"],
            c=100 * plot_df["inf_1"]/CONFIG["population_parameters"]["pop_size_1"])
cbar = plt.colorbar()
cbar.set_label("Infections in group 1 (%)")
plt.xlabel("a")
plt.ylabel("b")
#plt.savefig(f"{output_dir}/inf_group_1_across_all_perc.png", dpi=300)
#plt.savefig(f"{output_dir}/inf_group_1_across_all_perc.svg", dpi=300)

plt.clf()


variables = ["inf_1", "inf_2", "vac_1", "vac_2"]
labels = ["Infections in group 1 (%)",
          "Infections in group 2 (%)",
          "Vaccinations in group 1 (%)",
          "Vaccinations in group 2 (%)"]
colors = ["Reds", "Reds", "Purples", "Purples"]
for var, label, color in zip(variables, labels, colors):
    perc_var = "%s_perc"%var
    plot_df[perc_var] = 100 * plot_df[var]/CONFIG["population_parameters"]["pop_size_%s"%(var.split("_")[1])]
    plot_df["a"] = [round(i,2) for i in plot_df["a"]]
    plot_df["b"] = [round(i,2) for i in plot_df["b"]]
    data = plot_df.pivot(index="b", columns="a", values = perc_var)
    plt.figure()
    ax = sns.heatmap(data, linewidth=0.5,
                     vmin=min(plot_df[perc_var]),
                     vmax=max(plot_df[perc_var]),
                     #yticklabels=['High','Medium','Low','No'],
                     cbar_kws={'label': label},
                     cmap=color,annot=False, fmt='.0f')
    ax.invert_yaxis()
    #plt.xlabel("a")
    #plt.ylabel("b")
    plt.xlabel("Equity in Vaccination Multiplier (a)")
    plt.ylabel("Equity in Clinical Burden Multiplier (b)")

    plt.savefig(f"{output_dir}/hm_%s_across_all_perc.png"%var, bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_dir}/hm_%s_across_all_perc.svg"%var, bbox_inches='tight', dpi=300)


variables = ["inf_1", "inf_2", "vac_1", "vac_2"]
labels = ["Infections in group 1 (%)",
          "Infections in group 2 (%)",
          "Vaccinations in group 1 (%)",
          "Vaccinations in group 2 (%)"]
colors = ["Reds", "Reds", "Purples", "Purples"]
for var, label, color in zip(variables, labels, colors):
    fig, ax1 = plt.subplots()

    perc_var = "%s_perc"%var
    plot_df[perc_var] = 100 * plot_df[var]/CONFIG["population_parameters"]["pop_size_%s"%(var.split("_")[1])]
    plot_df["a"] = [round(i,2) for i in plot_df["a"]]
    plot_df["b"] = [round(i,2) for i in plot_df["b"]]
    x = plot_df["a"]
    y = plot_df["b"]
    z = plot_df[perc_var]


    # Create grid values first.
    xi = np.arange(0, 1, step/2)
    yi = np.arange(0, 1, step/2)

    # Linearly interpolate the data (x, y) on a grid defined by (xi, yi).
    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)

    # Note that scipy.interpolate provides means to interpolate data on a grid
    # as well. The following would be an alternative to the four lines above:
    # from scipy.interpolate import griddata
    # zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='linear')

    #ax1.contour(xi, yi, zi, levels=14, linewidths=0.5, colors='k')
    cntr1 = ax1.contourf(xi, yi, zi, levels=14, cmap=color)

    #fig.colorbar(cntr1, ax=ax1)
    #ax1.plot(x, y, 'ko', ms=3)
    ax1.set(xlim=(0,1), ylim=(0,1))
    cbar = plt.colorbar(cntr1)
    cbar.set_label(label)

    plt.xlabel("Equity in Vaccination Multiplier (a)")
    plt.ylabel("Equity in Clinical Burden Multiplier (b)")

    plt.savefig(f"{output_dir}/cnt_%s_across_all_perc.png"%var, bbox_inches='tight', dpi=300)
    plt.savefig(f"{output_dir}/cnt_%s_across_all_perc.svg"%var, bbox_inches='tight', dpi=300)
