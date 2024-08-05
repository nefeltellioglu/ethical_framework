import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
from dataclasses import dataclass
import itertools
import numpy as np
import pandas as pd
import pickle

from ethical_sir_fancy import (
    OptParams,
    BurdenParams,
    SIRParams,
    SIRInitialConditions,
    SIRSolution,
    optimal_initial_conditions,
    loss_clinical_burden,
    loss_equity_of_burden,
    loss_equity_of_vaccination,
    sir_vacc,
    sir_vacc_SSA,
)


model_parameters = [
    {
        "id": 0,
        "parameters": SIRParams(
            beta_11=0.19188309203678316,
            beta_12=0.16789770553218525,
            beta_21=0.16789770553218525,
            beta_22=0.14391231902758736,
            gamma=0.125,
        ),
    }
]


pop_size_1 = 1000
pop_size_2 = 500

initial_conditions = []
ic_ix = 0
for num_vac_1 in range(0, pop_size_1, 200):
    for num_vac_2 in range(0, pop_size_2, 200):
        # Print out the initial condition being added to the list
        print(f"Adding initial condition {ic_ix} with {num_vac_1} vaccinated in population 1 and {num_vac_2} vaccinated in population 2.")
        initial_conditions.append(
            {
                "id": ic_ix,
                "value": SIRInitialConditions(
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


configurations = [
    {"id": c_ix, "model_parameters_id": mp["id"], "initial_condition_id": ic["id"]}
    for c_ix, (mp, ic) in enumerate(
        itertools.product(model_parameters, initial_conditions)
    )
]


def _compute_sol(config) -> SIRSolution:
    model_params = next((mp["parameters"] for mp in model_parameters if mp["id"] == config["model_parameters_id"]), None)
    ic = next((ic["value"] for ic in initial_conditions if ic["id"] == config["initial_condition_id"]), None)
    return sir_vacc(params = model_params, sir_0 = ic, ts = np.linspace(0, 100, 100))[0]

solutions = [_compute_sol(c) for c in configurations]

num_seeds = 2
outcomes = [
    {
        "id": o_ix,
        "configuration_id": c["id"],
        "seed": seed,
        "inf_1_no_vac": sol.r1[-1],
        "inf_1_vu": sol.r1_vu[-1],
        "inf_1_vp": 0,
        "vac_1": sol.s1_vu[0] + sol.s1_vp[0],
        "inf_2_no_vac": sol.r2[-1],
        "inf_2_vu": sol.r2_vu[-1],
        "inf_2_vp": 0,
        "vac_2": sol.s2_vu[0] + sol.s2_vp[0],
    }
    for o_ix, (c, sol, seed) in enumerate(
        itertools.product(configurations, solutions, range(num_seeds))
    )
]

burden_parameters = [
    {
        "id": 0,
        "parameters": BurdenParams(
            perc_hosp_inf=0.02,
            days_hosp_inf_1=3.075,
            days_hosp_inf_2=7.60,
            perc_hosp_vacc_1=0.004,
            perc_hosp_vacc_2=0.002,
            days_hosp_vacc_1=6.0,
            days_hosp_vacc_2=6.0,
            vacc_protection_1=0.531,
            vacc_protection_2=0.531,
            vacc_protection_dis_1=0.825,
            vacc_protection_dis_2=0.825,
        ),
    }
]


db = {
    "model_parameters": model_parameters,
    "initial_conditions": initial_conditions,
    "configurations": configurations,
    "outcomes": outcomes,
    "burden_parameters": burden_parameters,
}

output_file = "grid_database.pkl"
with open(output_file, "wb") as f:
    pickle.dump(db, f)
