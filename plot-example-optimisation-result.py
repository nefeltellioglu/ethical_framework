import json
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

import ethics.model as em
import ethics.optimisation as eo


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


model_param_id = 0
burden_param_id = 0
ethical_a = 0.0
ethical_b = 0.95

bp = [bp for bp in db["burden_parameters"] if bp["id"] == burden_param_id][0][
    "parameters"
]

configs = [
    c for c in db["configurations"] if c["model_parameters_id"] == model_param_id
]
config_ids = [c["id"] for c in configs]
ocs = [o for o in db["outcomes"] if o["configuration_id"] in config_ids]

plot_df = []
for oc in ocs:
    tmp_config_id = oc["configuration_id"]

    tmp_ic = [
        ic
        for ic in db["initial_conditions"]
        if ic["id"]
        == [c for c in configs if c["id"] == tmp_config_id][0]["initial_condition_id"]
    ][0]["value"]


    tmp_loss_tcb, tmp_loss_ecb, tmp_loss_evb = em.loss_terms(oc["outcome"], tmp_ic, bp)
    tmp_loss = (
        (1 - ethical_a - ethical_b) * tmp_loss_tcb
        + ethical_a * tmp_loss_ecb
        + ethical_b * tmp_loss_evb
    )
    plot_df.append(
        {
            "outcome_id": oc["id"],
            "config_id": tmp_config_id,
            "vac_1": oc["outcome"].total_vac_1,
            "vac_2": oc["outcome"].total_vac_2,
            "loss": tmp_loss,
        }
    )
plot_df = pd.DataFrame(plot_df)

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

# Use matplotlib to make a plot of the same data.
# Put a single large red dot at the best outcome.

plt.figure()
plt.scatter(plot_df["vac_1"], plot_df["vac_2"], c=plot_df["loss"])
cbar = plt.colorbar()
plt.scatter(best_vac_1, best_vac_2, s=500, c="red", marker="o")
cbar.set_label("Loss")
plt.xlabel("Total Vaccinations in Group 1")
plt.ylabel("Total Vaccinations in Group 2")
plt.savefig(f"{output_dir}/example-optimisation-results.png")
plt.savefig(f"{output_dir}/example-optimisation-results.svg")
plt.clf()





plt.figure()
plt.scatter(100 * plot_df["vac_1"]/CONFIG["population_parameters"]["pop_size_1"],
            100 * plot_df["vac_2"]/CONFIG["population_parameters"]["pop_size_2"],
            c=plot_df["loss"])
cbar = plt.colorbar()
plt.scatter(100 * best_vac_1/CONFIG["population_parameters"]["pop_size_1"],
            100 * best_vac_2/CONFIG["population_parameters"]["pop_size_2"],
            s=500, c="red", marker="o")
cbar.set_label("Loss")
plt.xlabel("Total Vaccinations in Group 1 (%)")
plt.ylabel("Total Vaccinations in Group 2 (%)")
plt.savefig(f"{output_dir}/example-optimisation-results-perc.png")
plt.savefig(f"{output_dir}/example-optimisation-results-perc.svg")

plt.clf()
