import json
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import itertools 
import ethics.model as em
import ethics.optimisation as eo
import numpy as np

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
    config_file = "config/config-2024-10-14_manuscript.json"
    #config_file = "config/config-2024-10-28_limited_vaccine.json"
    config_file = "config/config-2024-12-02_limited_low_R0.json"
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




# ====================================================================
#find best vacc parameters
# ====================================================================

ethical_a_b_list = [(0.05, 0.05), (0.85, 0.05), (0.05, 0.85)]

times = 500
if "low" in config_file: 
   times = int(times * 20)
elif "high" in config_file: 
   times = int(times * 0.5)
   
selected_vaccinations = []
for (ethical_a, ethical_b) in ethical_a_b_list:
    #ethical_a = 0.1
    #ethical_b = 0.1
    
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
    
    selected_vaccinations.append((int(best_vac_1), int(best_vac_2)))



# ====================================================================
#calculate beta values
#contact_per_capita_ij are the rescaled values
# ====================================================================

contact_per_capita_11=CONFIG["model_parameters"]["contact_per_capita_11"]
contact_per_capita_12=CONFIG["model_parameters"]["contact_per_capita_12"]
contact_per_capita_21=CONFIG["model_parameters"]["contact_per_capita_21"]
contact_per_capita_22=CONFIG["model_parameters"]["contact_per_capita_22"]
gamma=CONFIG["model_parameters"]["gamma"]
R0=CONFIG["model_parameters"]["R0"]


#calculation of beta from R0 and contact_per_capita multipliers
beta = R0 * 2 * gamma / (contact_per_capita_11 + contact_per_capita_22 +
                         (contact_per_capita_11** 2 
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

#selected_vaccinations = [(int(pop_size_1 * 0.5), int(pop_size_2 * 0.5)), 
#                         (int(pop_size_1 * 0.0), int(pop_size_2 * 0.0))]
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
    return sir_vacc(params=model_params, sir_0=ic, ts=np.linspace(0, times, times + 1))[0]


solutions = [_compute_sol(c) for c in configurations]
_num_solutions = len(solutions)
assert _num_solutions == _num_configurations


# ====================================================================
#plot three trajectories
# ====================================================================

# NOTE the seed is always 0 because we are not using any randomness in
# the model

figsize=(10,2.5)#no in x axis, no in yaxis
fig, axs = plt.subplots(1, 3, figsize=figsize)
labels = ('A', 'B', 'C', 'D', 'E', 'F')
i = 0
    

for o_ix, (c, sol, (a,b)) in enumerate(zip(configurations, solutions, ethical_a_b_list)):
    
    ax = axs[i]
    ax.text(-0.25, 1.15, labels[i], transform=ax.transAxes,
      fontsize=12, fontweight='bold', va='top', ha='right')
    i += 1

    
    total_s1 = 100 * (sol.s1 + sol.s1_vp + sol.s1_vu) / pop_size_1
    total_s2 = 100 * (sol.s2 + sol.s2_vp + sol.s2_vu) / pop_size_2
    
    total_i1 = 100 * (sol.i1 + sol.i1_vu) / pop_size_1
    total_i2 = 100 * (sol.i2 + sol.i2_vu) / pop_size_2

    ax.plot(sol.times, total_s1,color = "tab:blue",label = "S1"
             )
    ax.plot(sol.times, total_s2,color = "tab:green",label = "S2"
             )
    ax.plot(sol.times, total_i1,color = "tab:red",label = "I1"
             )
    ax.plot(sol.times, total_i2,color = "tab:orange",label = "I2"
             )
    
    #ax.plot(sol.times, 100 * sol.s1_vp/ pop_size_1 ,color = "tab:purple",label = "test"
    #         )
    #ax.plot(sol.times, 100 * (pop_size_1 - (sol.s1 + sol.s1_vp + sol.s1_vu)
    #           - sol.r1 - sol.r1_vu)/ pop_size_1 ,color = "tab:purple",
    #        label = "test I1")
    
    vacc = [int(100 * x/y) for (x, y) in zip(selected_vaccinations[o_ix], (pop_size_1, pop_size_2))]
    
    ax.set_title(f'a = {a}, b = {b}', fontweight="bold"#, size = 8
                 )
    
    textstr = '\n'.join((
        'Vaccinations',
    f'Group 1: {vacc[0]}%',
    f'Group 2: {vacc[1]}%',
    ))

    props = dict(#boxstyle='round', 
                 facecolor='white', 
                 alpha=0.5)
    if (total_s1[-1] > 60) or (total_s2[-1] > 60):
        ax.text(0.975, 0.5, textstr, transform=ax.transAxes, #fontsize=14,
            verticalalignment='top', horizontalalignment='right', bbox=props)
    else:
        # place a text box in upper left in axes coords
        ax.text(0.975, 0.97, textstr, transform=ax.transAxes, #fontsize=14,
            verticalalignment='top', horizontalalignment='right', bbox=props)
    
    
    #ax.set_title(f'{vacc[0]}% Group 1 & {vacc[1]}% Group 2 vaccinated', fontweight="bold", size = 8)
    ax.set_xlabel('Time')
    ax.set_ylabel('Prevalence (% per group)')
    
    ax.legend().set_visible(False) 

ax.legend().set_visible(True) 
ax.legend(loc= "lower center", bbox_to_anchor=(-0.9,-0.55), ncol= 2)
plt.subplots_adjust(left=0.1,
            bottom=0.1, 
            right=0.9, 
            top=0.9, 
            wspace=0.4, 
            hspace=0.4)

fig.savefig(f"{output_dir}/tractectories.png", bbox_inches='tight', dpi=300)
#fig.savefig(f"{output_dir}/tractectories.svg", bbox_inches='tight', dpi=300)


# ====================================================================
#plot where opt vaccinations stays in the all vaccinations
# ====================================================================


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

# Use matplotlib to make a plot of the same data.
# Put a single large red dot at the best outcome.



figsize=(4,3)#no in x axis, no in yaxis
plt.figure( figsize=figsize)
labels = ('A', 'B', 'C', 'D', 'E', 'F')
    
plt.scatter(100 * plot_df["vac_1"]/CONFIG["population_parameters"]["pop_size_1"],
            100 * plot_df["vac_2"]/CONFIG["population_parameters"]["pop_size_2"],
            c=plot_df["loss"])
cbar = plt.colorbar()
cbar.set_label("Loss")
plt.xlabel("Total Vaccinations\nin Group 1 (%)")
plt.ylabel("Total Vaccinations\nin Group 2 (%)")

for o_ix, (c, sol, (a,b)) in enumerate(zip(configurations, solutions, ethical_a_b_list)):
    vacc = [int(100 * x/y) for (x, y) in zip(selected_vaccinations[o_ix], (pop_size_1, pop_size_2))]
    props = dict(boxstyle='round', 
                 facecolor='white', 
                 alpha=1)
    idx = [i for  i,v in enumerate(selected_vaccinations) 
           if v == selected_vaccinations[o_ix]]
    if len(idx) == 1:
        label = labels[idx[0]]
        coordx, coordy = vacc[0],vacc[1]
    else: 
        label = labels[o_ix]#labels[idx[-1]]
        coordx, coordy = vacc[0],vacc[1]
        #label = " - ".join([labels[i] for i in idx])
        if label != labels[idx[0]]:
            coordx, coordy = vacc[0] -0.3,vacc[1]
        
        
    #plt.annotate(label,
    #        xy=(vacc[0]/100,vacc[1]/100), xycoords='data',
    #        xytext=(vacc[0]/100,vacc[0]/100), textcoords='data',
    #        arrowprops=dict(facecolor='black', shrink=0.05),
    #        horizontalalignment='right', verticalalignment='top')
    plt.text(coordx, coordy,label ,bbox=props,c="black",weight="bold",#transform=ax.transAxes,
                #100 * best_vac_2/CONFIG["population_parameters"]["pop_size_2"],
                #s=500,  marker="o"
                )
plt.savefig(f"{output_dir}/example-optimisation-results-perc.png", 
            bbox_inches='tight', dpi=300)
#plt.savefig(f"{output_dir}/example-optimisation-results-perc.svg")

plt.clf()