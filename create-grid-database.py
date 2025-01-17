# Create a database of the solutions for our SIR model with varying
# levels of vaccination.
#
# The "database" created by this script is saved in a pickle file as a
# dictionary of lists (representing a set of tables). The schema for
# the database that is produced is in the README.md file. There is a
# graphical representation of the schema in database-schema.png.
#
import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
from dataclasses import dataclass
import itertools
import json
import numpy as np
import pandas as pd
import pickle
import os
import sys
import ethics.model as em


if len(sys.argv) > 1:
    config_file = sys.argv[1]
else:
    config_file = "config/config-2024-10-14_manuscript_CZ_test_II.json"
    #config_file = "config/config-2024-10-14_manuscript.json"

    #config_file = "config/config-2024-10-28_limited_vaccine.json"
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

# ------------------------------------------------------------------------------
# calculate beta values since the contact_per_capita_ij in the
# configuration file are the rescaled values.
# ------------------------------------------------------------------------------
contact_per_capita_11=CONFIG["model_parameters"]["contact_per_capita_11"]
contact_per_capita_12=CONFIG["model_parameters"]["contact_per_capita_12"]
contact_per_capita_21=CONFIG["model_parameters"]["contact_per_capita_21"]
contact_per_capita_22=CONFIG["model_parameters"]["contact_per_capita_22"]
gamma=CONFIG["model_parameters"]["gamma"]
R0=CONFIG["model_parameters"]["R0"]

beta = R0 * 2 * gamma / (contact_per_capita_11 + contact_per_capita_22 +
                         (contact_per_capita_11**2
                          - 2 * contact_per_capita_22 * contact_per_capita_11
                          + contact_per_capita_22 ** 2
                          + 4 * contact_per_capita_12 * contact_per_capita_22
                          )**(0.5))

model_parameters = [
    {
        "id": 0,
        "parameters": em.SIRParams(
            beta_11=beta * contact_per_capita_11,
            beta_12=beta * contact_per_capita_12,
            beta_21=beta * contact_per_capita_21,
            beta_22=beta * contact_per_capita_22,
            gamma=CONFIG["model_parameters"]["gamma"],
        ),
    }
]

_num_model_parameters = len(model_parameters)
assert _num_model_parameters == 1

vac_protection_from_inf = CONFIG["vacc_protection_from_infection"]
pop_size_1 = CONFIG["population_parameters"]["pop_size_1"]
pop_size_2 = CONFIG["population_parameters"]["pop_size_2"]

# --------------------------------------------------------------------
# Construct an iterator which contains all the possible combinations
# of the number of vaccinated individuals in each population and make
# a couple of assertions to ensure that the 0 and all vaccinated
# options are included. This is because we want to make sure out
# database covers a good range of possible options.
pop_grid_step_size = CONFIG["grid_search_step"]["grid_step"]

v1_step = int(np.ceil(pop_grid_step_size * pop_size_1))
_vac_vals_pop_1 = list(range(0, pop_size_1, v1_step))
if _vac_vals_pop_1[0] != 0 :
    _vac_vals_pop_1[0] = 0
if _vac_vals_pop_1[-1] != pop_size_1:
    _vac_vals_pop_1[-1] = pop_size_1

v2_step = int(np.ceil(pop_grid_step_size * pop_size_2))
_vac_vals_pop_2 = list(range(0, pop_size_2, v2_step))
if _vac_vals_pop_2[0] != 0 :
    _vac_vals_pop_2[0] = 0
if _vac_vals_pop_2[-1] != pop_size_2:
    _vac_vals_pop_2[-1] = pop_size_2

assert 0 in _vac_vals_pop_1, "Must have 0 vaccinated in population 1 as an option."
assert pop_size_1 in _vac_vals_pop_1, "Must have all vaccinated in population 1 as an option."
assert 0 in _vac_vals_pop_2, "Must have 0 vaccinated in population 2 as an option."
assert pop_size_2 in _vac_vals_pop_2, "Must have all vaccinated in population 2 as an option."
vac_12_combinations = itertools.product(_vac_vals_pop_1, _vac_vals_pop_2)
# --------------------------------------------------------------------

initial_conditions = []
for ic_ix, (num_vac_1, num_vac_2) in enumerate(vac_12_combinations):
    print(
        f"Adding initial condition {ic_ix} with {num_vac_1} vaccinated in population 1 and {num_vac_2} vaccinated in population 2."
    )
    s0_1_vp = int(num_vac_1 * vac_protection_from_inf)
    s0_2_vp = int(num_vac_2 * vac_protection_from_inf)
    initial_conditions.append(
        {
            "id": ic_ix,
            "value": em.SIRInitialCondition(
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
_num_initial_conditions = len(initial_conditions)

configurations = [
    {"id": c_ix, "model_parameters_id": mp["id"], "initial_condition_id": ic["id"]}
    for c_ix, (mp, ic) in enumerate(
        itertools.product(model_parameters, initial_conditions)
    )
]
_num_configurations = len(configurations)
assert _num_configurations == _num_model_parameters * _num_initial_conditions


def _compute_sol(config) -> em.SIRSolution:
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
    sim_times = np.linspace(
        CONFIG["time_parameters"]["start_time"],
        CONFIG["time_parameters"]["finish_time"],
        num = CONFIG["time_parameters"]["num_points"]
    )
    return em.sir_vacc(params=model_params, sir_0=ic, ts=sim_times)[0]


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
        "outcome": em.SIROutcome(
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


# --------------------------------------------------------------------
# Check that the total number of infected individuals at the end of
# the simulation is 0. If this is not the case, then the simulation
# has *probably* not run for long enough. In this case you should try
# extending the final time in the configuration file.
for sol in solutions:
    _total_i = round(sol.i1[-1], 0) + round(sol.i1_vu[-1],0) + \
            round(sol.i2[-1],0) + round(sol.i2_vu[-1],0)
    if round(_total_i, 0) > 0:
        raise ValueError(f"""
        Total infected at the end of the simulation is {_total_i} and should round to 0.
        This might be because you haven't run the simulation for long enough.
        Try increasing the value of the final time in the configuration file.
        """)
# --------------------------------------------------------------------


burden_parameters = [
    {
        "id": 0,
        "parameters": em.BurdenParams(
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
