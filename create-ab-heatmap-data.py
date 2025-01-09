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
#define function to calculate total burden from SIROutcome object
# burden from adverse vaccination reactions (group 1)
def burden_adverse_group_1(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    return (dbp.days_hosp_vacc_1
            * sir.total_vac_1
            * dbp.prop_hosp_vacc_1)

# burden from adverse vaccination reactions (group 2)
def burden_adverse_group_2(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    return (dbp.days_hosp_vacc_2
            * sir.total_vac_2
            * dbp.prop_hosp_vacc_2)

# count of infections
def count_infections_group_1(sir: em.SIROutcome) -> float:
    return(sir.inf_1_no_vac + sir.inf_1_vu)

def count_infections_group_2(sir: em.SIROutcome) -> float:
    return(sir.inf_2_no_vac + sir.inf_2_vu)

def count_vaccinations_group_1(sir: em.SIROutcome) -> float:
    return sir.total_vac_1

def count_vaccinations_group_2(sir: em.SIROutcome) -> float:
    return sir.total_vac_2

#burden from infections in unvaccinated people (group 1)
def burden_infections_group_1_noVacc(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    return (dbp.prop_hosp_inf_1 *
            dbp.days_hosp_inf_1 * sir.inf_1_no_vac)

# burden from infections in unvaccinated people (group 2)
def burden_infections_group_2_noVacc(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    return (dbp.prop_hosp_inf_2 *
            dbp.days_hosp_inf_2 * sir.inf_2_no_vac)

#burden from infections in vaccinated people (group 1)
def burden_infections_group_1_Vacc(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    return (dbp.prop_hosp_inf_1 *
            (1 - dbp.vacc_protection_from_disease_1) *
            dbp.days_hosp_inf_1 *
            sir.inf_1_vu )

#burden from infections in vaccinated people (group 2)
def burden_infections_group_2_Vacc(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    return (dbp.prop_hosp_inf_2 *
            (1 - dbp.vacc_protection_from_disease_2) *
            dbp.days_hosp_inf_2 *
            sir.inf_2_vu )

# total infection burden group 1
def total_burden_infections_group_1(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    tot= (burden_infections_group_1_noVacc(sir, dbp) +
            burden_infections_group_1_Vacc(sir, dbp))
    return (tot)

def total_burden_infections_group_2(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    tot= (burden_infections_group_2_noVacc(sir, dbp) +
            burden_infections_group_2_Vacc(sir, dbp))
    return (tot)


def total_vaccinations(sir: em.SIROutcome) -> float:
    return (sir.total_vac_1 + sir.total_vac_2)

# aggregate burden components
def total_burden_infections(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    tot_1 = total_burden_infections_group_1(sir, dbp)
    tot_2 = total_burden_infections_group_2(sir, dbp)
    return (tot_1 + tot_2)

def total_burden_adverse(sir: em.SIROutcome, dbp: em.BurdenParams) -> float:
    return(burden_adverse_group_1(sir, dbp) +
           burden_adverse_group_2(sir, dbp))

def total_clinical_burden(sir:em.SIROutcome, dbp: em.BurdenParams) -> float:
    return (total_burden_infections(sir, dbp) +
            total_burden_adverse(sir, dbp))


# per-capita burdens:
def pop_1(ic: em.SIRInitialCondition) -> int:
    return ic.pop_size(1)

def pop_2(ic: em.SIRInitialCondition) -> int:
    return ic.pop_size(2)

def total_pop(ic: em.SIRInitialCondition) -> int:
    return (pop_1(ic) + pop_2(ic))

def adverse_per_capita_1(sir: em.SIROutcome,
                         dbp: em.BurdenParams,
                         ic: em.SIRInitialCondition) -> float:
    return(burden_adverse_group_1(sir, dbp) / pop_1(ic))

def adverse_per_capita_2(sir: em.SIROutcome,
                         dbp: em.BurdenParams,
                         ic: em.SIRInitialCondition) -> float:
    return(burden_adverse_group_2(sir, dbp) / pop_2(ic))

def infection_burden_per_capita_1(sir: em.SIROutcome,
                                  dbp: em.BurdenParams,
                                  ic: em.SIRInitialCondition) -> float:
    return(total_burden_infections_group_1(sir, dbp) / pop_1(ic))

def infection_burden_per_capita_2(sir: em.SIROutcome,
                                  dbp: em.BurdenParams,
                                  ic: em.SIRInitialCondition) -> float:
    return(total_burden_infections_group_2(sir, dbp) / pop_2(ic))




#### ====================================================================











































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
step = CONFIG["grid_search_step"]["a_b_grid_step"]
grid_max = CONFIG["grid_search_step"]["a_b_grid_max"]
grid_min = CONFIG["grid_search_step"]["a_b_grid_min"]


plot_df = []

for ethical_a in np.arange(grid_min, grid_max, step):
    for ethical_b in np.arange(grid_min, grid_max - ethical_a , step):

        foo, _ = eo.optimal_initial_condition(
            ethical_a, ethical_b, model_param_id, burden_param_id, db, normalise=True
        )
        _optimal_ic = [ic for ic in db["initial_conditions"] if ic["id"] == foo]
        _optimal_config = [c for c in db["configurations"] if c["initial_condition_id"] == foo]
        _optimal_outcome = [
            o for o in db["outcomes"] if o["configuration_id"] == _optimal_config[0]["id"]
        ]
        ic_ab = _optimal_ic[0]["value"]
        oc_ab = _optimal_outcome[0]["outcome"]

        # this is where outcomes are logged
        vac_1 = count_vaccinations_group_1(oc_ab) / pop_1(ic_ab)
        vac_2 = count_vaccinations_group_2(oc_ab) / pop_2(ic_ab)
        inf_1 = count_infections_group_1(oc_ab) / pop_1(ic_ab)
        inf_2 = count_infections_group_2(oc_ab) / pop_2(ic_ab)

        clinical_burden = total_clinical_burden(oc_ab, bp)
        infections_burden = total_burden_infections(oc_ab, bp)
        adverse_burden = total_burden_adverse(oc_ab, bp)
        total_vaccination = total_vaccinations(oc_ab)
        total_vaccination_perc = 100 * total_vaccination / (pop_1(ic_ab) + pop_2(ic_ab))


        plot_df.append(
            {   "a": ethical_a,
                "b": ethical_b,
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
