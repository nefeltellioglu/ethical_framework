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
import itertools

from ethics.model import (
    OptParams,
    BurdenParams,
    SIRParams,
    SIRInitialCondition,
    SIRSolution,
    SIROutcome,
    sir_vacc,
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


# ====================================================================
#calculate beta values
#contact_per_capita_ij are the rescaled values
contact_per_capita_11=CONFIG["model_parameters"]["contact_per_capita_11"]
contact_per_capita_12=CONFIG["model_parameters"]["contact_per_capita_12"]
contact_per_capita_21=CONFIG["model_parameters"]["contact_per_capita_21"]
contact_per_capita_22=CONFIG["model_parameters"]["contact_per_capita_22"]
gamma=CONFIG["model_parameters"]["gamma"]
R0=CONFIG["model_parameters"]["R0"]


#calculation of beta from R0 and contact_per_capita multipliers
beta = R0 * 2 * gamma / (contact_per_capita_11 + contact_per_capita_22 +
                         (contact_per_capita_11**2 
                          - 2 * contact_per_capita_22 * contact_per_capita_11
                          + contact_per_capita_22 ** 2
                          + 4 * contact_per_capita_12 * contact_per_capita_22
                          )**(0.5))


CONFIG["model_parameters"]["beta_11"] = beta * contact_per_capita_11
CONFIG["model_parameters"]["beta_12"] = beta * contact_per_capita_12
CONFIG["model_parameters"]["beta_21"] = beta * contact_per_capita_21
CONFIG["model_parameters"]["beta_22"] = beta * contact_per_capita_22

########################




model_parameters = [
    {
        "id": 0,
        "parameters": SIRParams(
            beta_11=CONFIG["model_parameters"]["beta_11"],
            beta_12=CONFIG["model_parameters"]["beta_12"],
            beta_21=CONFIG["model_parameters"]["beta_21"],
            beta_22=CONFIG["model_parameters"]["beta_22"],
            gamma=CONFIG["model_parameters"]["gamma"],
        ),
    }
]

_num_model_parameters = len(model_parameters)
assert _num_model_parameters == 1

pop_size_1 = CONFIG["population_parameters"]["pop_size_1"]
pop_size_2 = CONFIG["population_parameters"]["pop_size_2"]
vac_protection_from_inf = CONFIG["vacc_protection_from_infection"]

initial_conditions = []
ic_ix = 0

selected_vaccinations = [(int(pop_size_1 * 0.5), int(pop_size_2 * 0.5)), 
                         (int(pop_size_1 * 0.2), int(pop_size_2 * 0.8))]
for (num_vac_1, num_vac_2) in selected_vaccinations:
        # Print out the initial condition being added to the list
        print(
            f"Adding initial condition {ic_ix} with {num_vac_1} vaccinated in population 1 and {num_vac_2} vaccinated in population 2."
        )
        s0_1_vp = int(num_vac_1 * vac_protection_from_inf)
        s0_2_vp = int(num_vac_2 * vac_protection_from_inf)
        initial_conditions.append(
            {
                "id": ic_ix,
                "value": SIRInitialCondition(
                    s0_1=pop_size_1 - num_vac_1 - 1,
                    s0_2=pop_size_2 - num_vac_2 - 1,
                    i0_1=1,
                    i0_2=1,
                    r0_1=0,
                    r0_2=0,
                    s0_1_vp= s0_1_vp,
                    s0_2_vp= s0_2_vp,
                    s0_1_vu= num_vac_1 - s0_1_vp,
                    s0_2_vu= num_vac_2 - s0_2_vp,
                    i0_1_vu=0,
                    i0_2_vu=0,
                    r0_1_vu=0,
                    r0_2_vu=0,
                ),
            }
        )
        ic_ix += 1
_num_initial_conditions = len(initial_conditions)
assert _num_initial_conditions == (len(selected_vaccinations))

configurations = [
    {"id": c_ix, "model_parameters_id": mp["id"], "initial_condition_id": ic["id"]}
    for c_ix, (mp, ic) in enumerate(
        itertools.product(model_parameters, initial_conditions)
    )
]
_num_configurations = len(configurations)
assert _num_configurations == _num_model_parameters * _num_initial_conditions


def _compute_sol(config) -> SIRSolution:
    model_params = next(
        (
            mp["parameters"]
            for mp in model_parameters
            if mp["id"] == config["model_parameters_id"]
        ),
        None,
    )
    ic = next(
        (
            ic["value"]
            for ic in initial_conditions
            if ic["id"] == config["initial_condition_id"]
        ),
        None,
    )
    return sir_vacc(params=model_params, sir_0=ic, ts=np.linspace(0, 100, 100))[0]


solutions = [_compute_sol(c) for c in configurations]
_num_solutions = len(solutions)
assert _num_solutions == _num_configurations

# NOTE the seed is always 0 because we are not using any randomness in
# the model
outcomes = [
    {
        "id": o_ix,
        "configuration_id": c["id"],
        "seed": 0,
        "outcome": SIROutcome(
            inf_1_no_vac=sol.r1[-1],
            inf_1_vu=sol.r1_vu[-1],
            inf_1_vp=0,
            total_vac_1=sol.s1_vu[0] + sol.s1_vp[0],
            inf_2_no_vac=sol.r2[-1],
            inf_2_vu=sol.r2_vu[-1],
            inf_2_vp=0,
            total_vac_2=sol.s2_vu[0] + sol.s2_vp[0],
        ),
    }
    for o_ix, (c, sol) in enumerate(zip(configurations, solutions))
]
_num_outcomes = len(outcomes)
assert _num_outcomes == _num_solutions


burden_parameters = [
    {
        "id": 0,
        "parameters": BurdenParams(
            prop_hosp_inf_1=CONFIG["burden_parameters"]["prop_hosp_inf_1"],
            prop_hosp_inf_2=CONFIG["burden_parameters"]["prop_hosp_inf_2"],
            days_hosp_inf_1=CONFIG["burden_parameters"]["days_hosp_inf_1"],
            days_hosp_inf_2=CONFIG["burden_parameters"]["days_hosp_inf_2"],
            prop_hosp_vacc_1=CONFIG["burden_parameters"]["prop_hosp_vacc_1"],
            prop_hosp_vacc_2=CONFIG["burden_parameters"]["prop_hosp_vacc_2"],
            days_hosp_vacc_1=CONFIG["burden_parameters"]["days_hosp_vacc_1"],
            days_hosp_vacc_2=CONFIG["burden_parameters"]["days_hosp_vacc_2"],
            vacc_protection_from_disease_1=CONFIG["burden_parameters"][
                "vacc_protection_from_disease_1"
            ],
            vacc_protection_from_disease_2=CONFIG["burden_parameters"][
                "vacc_protection_from_disease_2"
            ],
        ),
    }
]
_num_burden_parameters = len(burden_parameters)
assert _num_burden_parameters == 1

db = {
    "model_parameters": model_parameters,
    "initial_conditions": initial_conditions,
    "configurations": configurations,
    "outcomes": outcomes,
    "burden_parameters": burden_parameters,
}

output_file = CONFIG["database_file"]
with open(output_file, "wb") as f:
    pickle.dump(db, f)
