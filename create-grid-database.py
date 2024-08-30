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

from ethics.model import (
    OptParams,
    BurdenParams,
    SIRParams,
    SIRInitialCondition,
    SIRSolution,
    SIROutcome,
    optimal_initial_conditions,
    loss_clinical_burden,
    loss_equity_of_burden,
    loss_equity_of_vaccination,
    sir_vacc,
    sir_vacc_SSA,
)


if len(sys.argv) > 1:
    config_file = sys.argv[1]
else:
    config_file = "config/config-2024-08-07.json"
assert os.path.exists(config_file)

with open(config_file, "r") as f:
    CONFIG = json.load(f)


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

initial_conditions = []
ic_ix = 0
for num_vac_1 in range(0, pop_size_1, 100):
    for num_vac_2 in range(0, pop_size_2, 100):
        # Print out the initial condition being added to the list
        print(
            f"Adding initial condition {ic_ix} with {num_vac_1} vaccinated in population 1 and {num_vac_2} vaccinated in population 2."
        )
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
                    s0_1_vp=num_vac_1,
                    s0_2_vp=num_vac_2,
                    s0_1_vu=0,
                    s0_2_vu=0,
                    i0_1_vu=0,
                    i0_2_vu=0,
                    r0_1_vu=0,
                    r0_2_vu=0,
                ),
            }
        )
        ic_ix += 1
_num_initial_conditions = len(initial_conditions)
assert _num_initial_conditions == (len(range(0, pop_size_1, 100)) * len(range(0, pop_size_2, 100)))

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
    for o_ix, (c, sol) in enumerate(
        zip(configurations, solutions)
    )
]
_num_outcomes = len(outcomes)
assert _num_outcomes == _num_solutions


burden_parameters = [
    {
        "id": 0,
        "parameters": BurdenParams(
            perc_hosp_inf=CONFIG["burden_parameters"]["perc_hosp_inf"],
            days_hosp_inf_1=CONFIG["burden_parameters"]["days_hosp_inf_1"],
            days_hosp_inf_2=CONFIG["burden_parameters"]["days_hosp_inf_2"],
            perc_hosp_vacc_1=CONFIG["burden_parameters"]["perc_hosp_vacc_1"],
            perc_hosp_vacc_2=CONFIG["burden_parameters"]["perc_hosp_vacc_2"],
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
