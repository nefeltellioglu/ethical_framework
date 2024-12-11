import json
import pickle
import matplotlib.pyplot as plt
import sys
import os

from ethics.model import (
    OptParams,
    BurdenParams,
    SIRParams,
    SIRSolution,
    SIROutcome,
    sir_vacc
)

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

configs = db["configurations"]
m_params = db["model_parameters"]
ics = db["initial_conditions"]

tmp = {"total_vacc_1": [], "total_vacc_2": [], "total_inf_1": [], "total_inf_2": []}
for oc in db["outcomes"]:
    oc_obj = oc["outcome"]
    tmp["total_vacc_1"].append(oc_obj.total_vac_1)
    tmp["total_vacc_2"].append(oc_obj.total_vac_2)
    tmp["total_inf_1"].append(oc_obj.inf_1_no_vac + oc_obj.inf_1_vu + oc_obj.inf_1_vp)
    tmp["total_inf_2"].append(oc_obj.inf_2_no_vac + oc_obj.inf_2_vu + oc_obj.inf_2_vp)


plt.figure()
plt.scatter(tmp["total_vacc_1"], tmp["total_inf_1"], c=tmp["total_vacc_2"])
cbar = plt.colorbar()
cbar.set_label("Total Vaccinations in Group 2")
plt.xlabel("Total Vaccinations in Group 1")
plt.ylabel("Total Infections in Group 1")
plt.savefig(f"{output_dir}/vacc-vs-inf-group-1.png")
plt.savefig(f"{output_dir}/vacc-vs-inf-group-1.svg")
plt.clf()


plt.figure()
plt.scatter(tmp["total_vacc_2"], tmp["total_inf_2"], c=tmp["total_vacc_1"])
cbar = plt.colorbar()
cbar.set_label("Total Vaccinations in Group 1")
plt.xlabel("Total Vaccinations in Group 2")
plt.ylabel("Total Infections in Group 2")
plt.savefig(f"{output_dir}/vacc-vs-inf-group-2.png")
plt.savefig(f"{output_dir}/vacc-vs-inf-group-2.svg")
plt.clf()


plt.figure()
plt.scatter([100 * i/CONFIG["population_parameters"]["pop_size_1"] for i in tmp["total_vacc_1"]],
            [100 * i/CONFIG["population_parameters"]["pop_size_1"] for i in tmp["total_inf_1"]],
            c=[100 * i/CONFIG["population_parameters"]["pop_size_2"] for i in tmp["total_vacc_2"]])
cbar = plt.colorbar()
cbar.set_label("Total Vaccinations in Group 2 (%)")
plt.xlabel("Total Vaccinations in Group 1 (%)")
plt.ylabel("Total Infections in Group 1 (%)")
plt.savefig(f"{output_dir}/vacc-vs-inf-group-1-perc.png")
plt.savefig(f"{output_dir}/vacc-vs-inf-group-1-perc.svg")
plt.clf()


plt.figure()
plt.scatter([100 * i/CONFIG["population_parameters"]["pop_size_2"] for i in tmp["total_vacc_2"]],
            [100 * i/CONFIG["population_parameters"]["pop_size_2"] for i in tmp["total_inf_2"]],
            c=[100 * i/CONFIG["population_parameters"]["pop_size_1"] for i in tmp["total_vacc_1"]])
cbar = plt.colorbar()
cbar.set_label("Total Vaccinations in Group 1 (%)")
plt.xlabel("Total Vaccinations in Group 2 (%)")
plt.ylabel("Total Infections in Group 2 (%)")
plt.savefig(f"{output_dir}/vacc-vs-inf-group-2-perc.png")
plt.savefig(f"{output_dir}/vacc-vs-inf-group-2-perc.svg")
plt.clf()
