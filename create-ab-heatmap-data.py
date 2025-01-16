# This script will generate a CSV file that has all the data needed to
# generate the AB heatmaps.
#
# For example: ./out/2024-10-14_manuscript/ab-heatmap-vaccination-part-1.png
#
# This code is takes a non-trivial amount of time to run, so it is
# nice to be able to generate the CSV once and then plot different
# slices of it in other scripts.
#
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
    #config_file = "config/config-2024-10-14_manuscript_CZ_test_II.json"
    #config_file = "config/config-2024-10-14_manuscript.json"
     config_file = "config/config-2024-10-28_limited_vaccine.json"
    #config_file = "config/config-2024-12-02_limited_low_R0.json"
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

# specify output dir
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

#### ====================================================================

model_param_id = 0
burden_param_id = 0
bp = [bp for bp in db["burden_parameters"] if bp["id"] == burden_param_id][0][
    "parameters"
]
configs = [
    c for c in db["configurations"] if c["model_parameters_id"] == model_param_id
]
config_ids = [c["id"] for c in configs]
ocs = [o for o in db["outcomes"] if o["configuration_id"] in config_ids]
assert "a_b_grid_step" not in CONFIG["grid_search_step"], """

Remove 'a_b_grid_step' from the configuration file.
This is not supported.
Instead set 'a_b_grid_num_points'.

"""
num_points = CONFIG["grid_search_step"]["a_b_grid_num_points"]
grid_max = CONFIG["grid_search_step"]["a_b_grid_max"]
grid_min = CONFIG["grid_search_step"]["a_b_grid_min"]

plot_df = []

eth_a_vals = np.linspace(grid_min, grid_max, num_points)
eth_b_vals = np.linspace(grid_min, grid_max, num_points)
for eth_a, eth_b in itertools.product(eth_a_vals, eth_b_vals):
    if eth_a + eth_b > 1:
        continue

    print('optimising vaccination for w_EI = ' + str(eth_a) + ' w_EV = ' + str(eth_b))

    opt_ix, _ = eo.optimal_initial_condition(
        eth_a, eth_b, model_param_id, burden_param_id, db, normalise=True
    )
    _optimal_ic = [ic for ic in db["initial_conditions"] if ic["id"] == opt_ix]
    _optimal_config = [c for c in db["configurations"] if c["initial_condition_id"] == opt_ix]
    _optimal_outcome = [
        o for o in db["outcomes"] if o["configuration_id"] == _optimal_config[0]["id"]
    ]
    ic_ab = _optimal_ic[0]["value"]
    oc_ab = _optimal_outcome[0]["outcome"]
    vac_1 = oc_ab.total_vac_1 / ic_ab.pop_size(1)
    vac_2 = oc_ab.total_vac_2 / ic_ab.pop_size(2)
    inf_1 = (oc_ab.inf_1_no_vac + oc_ab.inf_1_vu) / ic_ab.pop_size(1)
    inf_2 = (oc_ab.inf_2_no_vac + oc_ab.inf_2_vu) / ic_ab.pop_size(2)
    clinical_burden = em.total_clinical_burden(oc_ab, bp)
    infections_burden = em.total_burden_infections(oc_ab, bp)
    adverse_burden = em.total_burden_adverse(oc_ab, bp)
    total_vaccination = em.total_vaccinations(oc_ab)
    total_vaccination_perc = 100 * total_vaccination / (ic_ab.pop_size(1) + ic_ab.pop_size(2))
    plot_df.append(
        {   "a": eth_a,
            "b": eth_b,
            "vac_1": vac_1,
            "vac_2": vac_2,
            "inf_1": inf_1,
            "inf_2": inf_2,
            "cli_burden": clinical_burden,
            "inf_burden": infections_burden,
            "adv_burden": adverse_burden,
            "total_vacc": total_vaccination,
            "total_vacc_perc": total_vaccination_perc,
         }
    )

plot_df = pd.DataFrame(plot_df)
plot_df.to_csv(f"{output_dir}/ab-heatmap-data.csv", index=False)
