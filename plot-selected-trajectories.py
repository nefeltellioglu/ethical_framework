import json
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
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

from ethics.optimisation import normalisation, get_extreme_burdens



if len(sys.argv) > 1:
    config_file = sys.argv[1]
else:
    #config_file = "config/config-2024-10-14_manuscript.json"
    config_file = "config/config-2024-10-28_limited_vaccine.json"
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

fig, axs = plt.subplots(1, 3, figsize=(10, 2.5))
subplot_labels = ['A', 'B', 'C']

# TODO These colour codes should not be hard-coded, but they are just
# colorbrewer2 8-class Dark2 values so they shouldn't be too
# mysterious for now.
green_hex = "#1b9e77"
orange_hex = "#d95f02"

for ix, (c, sol, (a,b)) in enumerate(zip(configurations, solutions, ethical_a_b_list)):
    ax = axs[ix]
    ax.text(-0.25, 1.15, subplot_labels[ix], transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='right')

    total_s1 = 100 * (sol.s1 + sol.s1_vp + sol.s1_vu) / pop_size_1
    total_s2 = 100 * (sol.s2 + sol.s2_vp + sol.s2_vu) / pop_size_2
    total_i1 = 100 * (sol.i1 + sol.i1_vu) / pop_size_1
    total_i2 = 100 * (sol.i2 + sol.i2_vu) / pop_size_2

    ax.plot(sol.times, total_s1, color = green_hex, linestyle="solid", label = "Group 1 Susceptibles")
    ax.plot(sol.times, total_s2, color = green_hex, linestyle="dashed", label = "Group 2 Susceptibles")
    ax.plot(sol.times, total_i1, color = orange_hex, linestyle="solid", label = "Group 1 Infecteds")
    ax.plot(sol.times, total_i2, color = orange_hex, linestyle="dashed", label = "Group 2 Infecteds")

    vacc = [int(100 * x/y) for (x, y) in zip(selected_vaccinations[ix], (pop_size_1, pop_size_2))]

    ax.set_title(f'a = {a}, b = {b}', fontweight="bold"#, size = 8
                 )

    textstr = '\n'.join((
        'Optimal vaccination',
    f'Group 1: {vacc[0]}%',
    f'Group 2: {vacc[1]}%',
    ))

    props = {"facecolor":"white", "alpha":1.0}
    if (total_s1[-1] > 60) or (total_s2[-1] > 60):
        ax.text(0.315, 0.5, textstr, transform=ax.transAxes, #fontsize=14,
            verticalalignment='top', horizontalalignment='left', bbox=props)
    else:
        # place a text box in upper left in axes coords
        ax.text(0.315, 0.97, textstr, transform=ax.transAxes, #fontsize=14,
            verticalalignment='top', horizontalalignment='left', bbox=props)

    ax.set_xlabel('Day')
    ax.set_ylabel('Group percentage (%)')
    ax.legend().set_visible(False)

ax.legend().set_visible(True)
ax.legend(loc= "lower center", bbox_to_anchor=(-0.9,-0.55), ncol= 2)
plt.subplots_adjust(left=0.1,
            bottom=0.1,
            right=0.9,
            top=0.9,
            wspace=0.4,
            hspace=0.4)

fig.savefig(f"{output_dir}/trajectories.png", bbox_inches='tight', dpi=300)
fig.savefig(f"{output_dir}/trajectories.svg", bbox_inches='tight')


# ====================================================================
#plot where opt vaccinations stays in the all vaccinations
# ====================================================================


plot_df_list = []

selected_vaccinations = []
for (ethical_a, ethical_b) in ethical_a_b_list:
    #ethical_a = 0.1
    #ethical_b = 0.1
    plot_df = []

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
    extreme_burdens = get_extreme_burdens(model_param_id, burden_param_id, db)
    
    for oc in ocs:
        tmp_config_id = oc["configuration_id"]

        tmp_ic = [
            ic
            for ic in db["initial_conditions"]
            if ic["id"]
            == [c for c in configs if c["id"] == tmp_config_id][0]["initial_condition_id"]
        ][0]["value"]
        

        inf_1 = count_infections_group_1(oc["outcome"]) / pop_1(tmp_ic)
        inf_2 = count_infections_group_2(oc["outcome"]) / pop_2(tmp_ic)
        
        clinical_burden = total_clinical_burden(oc["outcome"], bp)
        infections_burden = total_burden_infections(oc["outcome"], bp)
        adverse_burden = total_burden_adverse(oc["outcome"], bp)
        
        total_vaccination = total_vaccinations(oc["outcome"])

        p_adverse_1 = adverse_per_capita_1(oc["outcome"], bp, tmp_ic)
        p_adverse_2 = adverse_per_capita_2(oc["outcome"], bp, tmp_ic)

        p_infection_burden_1 = infection_burden_per_capita_1(oc["outcome"], bp, tmp_ic)
        p_infection_burden_2 = infection_burden_per_capita_2(oc["outcome"], bp, tmp_ic)
        
        (tmp_loss_tcb, tmp_loss_ecb, tmp_loss_evb,
        tmp_loss_tcb_nonnorm, tmp_loss_ecb_nonnorm, tmp_loss_evb_nonnorm) = normalisation( 
            ethical_a,
            ethical_b,
            oc["outcome"], 
            tmp_ic, 
            bp, 
            extreme_burdens)
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
                "inf_1": inf_1,
                "inf_2": inf_2,
                "cli_burden": clinical_burden,
                "inf_burden": infections_burden,
                "adv_burden": adverse_burden,
                "total_vacc": total_vaccination,
                "p_adverse_1": p_adverse_1,
                "p_adverse_2": p_adverse_2,
                "p_infburd_1": p_infection_burden_1,
                "p_infburd_2": p_infection_burden_2,
                "loss": tmp_loss,
                "loss_tcb": tmp_loss_tcb,
                "loss_eib": tmp_loss_ecb,
                "loss_evb": tmp_loss_evb,
                "tcb": tmp_loss_tcb_nonnorm, 
                "teib": tmp_loss_ecb_nonnorm,
                "tevb": tmp_loss_evb_nonnorm,
            }
        )
    
    plot_df = pd.DataFrame(plot_df)
    plot_df_list.append(plot_df)

# Use matplotlib to make a plot of the same data.
# Put a single large red dot at the best outcome.


plot_df = plot_df_list[0]
myvars_list = [["cli_burden", "inf_burden", "adv_burden"],
          ["loss_tcb", "loss_eib", "loss_evb"],
          ["tcb", "teib", "tevb"]]
mylabels_list = [["Total clinical burden", "Infection burden", "Vaccination burden"],
          [f"$L_{{CB}}$", f"$L_{{EI}}$", f"$L_{{EV}}$"],
          ["Loss in clinical burden", 
           "Loss in equity in\ninfection burden", 
           "Loss in equity in\nvaccination burden"]]
fnames = ["Burden", "Normalized_individual_loss", "Individual_loss"]
#myvars = ["cli_burden","teib"]
#myvars = ["loss"]
for myvars, mylabels, fname in zip(myvars_list, mylabels_list,fnames):
    fig, axs = plt.subplots(1, 3, figsize=(10,4))
    subplot_labels1 = ['a', 'b', 'c']
    subplot_labels = ['A', 'B', 'C']
    
    
    for ix, (myvar, mylabel) in enumerate(zip(myvars, mylabels)):
        ax = axs[ix]
        ax.text(-0.25, 1.05, subplot_labels1[ix], transform=ax.transAxes,
                fontsize=12, fontweight='bold', va='top', ha='right')
        cax1 = ax.scatter(100 * plot_df["vac_1"]/CONFIG["population_parameters"]["pop_size_1"],
                    100 * plot_df["vac_2"]/CONFIG["population_parameters"]["pop_size_2"],
                    c=plot_df[myvar], cmap = "viridis_r"#, norm=mpl.colors.LogNorm()
                    )
        x_max = max(100 * plot_df["vac_1"]/CONFIG["population_parameters"]["pop_size_1"])
        cbar = fig.colorbar(cax1, ax=ax, location ='bottom',
                  pad = 0.23)
        
        #cbar = plt.colorbar(ax = ax)
        """if mylabel in ["Loss in equity in\ninfection burden", 
                       "Loss in equity in\nvaccination burden"]:
            cbar.set_label(mylabel,rotation=270, labelpad = 20, )
        else:
            cbar.set_label(mylabel,rotation=270, labelpad = 10, )"""
        cbar.set_label(mylabel)
        ax.set_xlabel("Total Vaccinations\nin Group 1 (%)")
        ax.set_ylabel("Total Vaccinations\nin Group 2 (%)")
        
        
        for o_ix, (c, sol, (a,b)) in enumerate(zip(configurations, solutions, ethical_a_b_list)):
            vacc = [int(100 * x/y) for (x, y) in zip(selected_vaccinations[o_ix], (pop_size_1, pop_size_2))]
            props = dict(boxstyle='round', 
                         facecolor='white', 
                         alpha=1)
            idx = [i for  i,v in enumerate(selected_vaccinations) 
                   if v == selected_vaccinations[o_ix]]
            if len(idx) == 1:
                label = subplot_labels[idx[0]]
                coordx, coordy = vacc[0],vacc[1]
            else: 
                label = subplot_labels[o_ix]#labels[idx[-1]]
                coordx, coordy = vacc[0],vacc[1]
                label = " - ".join([subplot_labels[i] for i in idx])
                if label != subplot_labels[idx[0]]:
                    coordx, coordy = vacc[0] -1/10 * x_max,vacc[1]
            
            
            ax.text(coordx, coordy,label ,bbox=props,c="black",weight="bold",#transform=ax.transAxes,
                        #100 * best_vac_2/CONFIG["population_parameters"]["pop_size_2"],
                        #s=500,  marker="o"
                        )
    
    plt.subplots_adjust(left=0.1,
                bottom=0.1, 
                right=0.9, 
                top=0.9, 
                wspace=0.4, 
                hspace=0.6)
    plt.savefig(f"{output_dir}/example-optimisation-%s.png"%fname, 
                bbox_inches='tight', dpi=300)
    #plt.savefig(f"{output_dir}/example-optimisation-results-perc.svg")
    #plt.show()
    #plt.clf()
    
fname = "Aggregated_loss"
myvar = "loss"
mylabel = "$L$"
fig, axs = plt.subplots(1, 3, figsize=(10, 4))
subplot_labels = ['A', 'B', 'C']


for ix, plot_df in enumerate(plot_df_list):
    ax = axs[ix]
    ax.text(-0.25, 1.05, subplot_labels[ix], transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='right')
    cax1 = ax.scatter(100 * plot_df["vac_1"]/CONFIG["population_parameters"]["pop_size_1"],
                100 * plot_df["vac_2"]/CONFIG["population_parameters"]["pop_size_2"],
                c=plot_df[myvar], cmap = "viridis_r"#, norm=mpl.colors.LogNorm()
                )
    x_max = max(100 * plot_df["vac_1"]/CONFIG["population_parameters"]["pop_size_1"])
    cbar = fig.colorbar(cax1, ax=ax, location ='bottom',
              pad = 0.23)
    #cbar = plt.colorbar(ax = ax)
    cbar.set_label(mylabel)
    ax.set_xlabel("Total Vaccinations\nin Group 1 (%)")
    ax.set_ylabel("Total Vaccinations\nin Group 2 (%)")
    
    vacc = [int(100 * x/y) for (x, y) in zip(selected_vaccinations[ix], (pop_size_1, pop_size_2))]
    
    ax.text(vacc[0],vacc[1],subplot_labels[ix],bbox=props,c="black",weight="bold",#transform=ax.transAxes,
                #100 * best_vac_2/CONFIG["population_parameters"]["pop_size_2"],
                #s=500,  marker="o"
                )

plt.subplots_adjust(left=0.1,
            bottom=0.1, 
            right=0.9, 
            top=0.9, 
            wspace=0.4, 
            hspace=0.6)
plt.savefig(f"{output_dir}/example-optimisation-%s.png"%fname, 
            bbox_inches='tight', dpi=300)
#plt.savefig(f"{output_dir}/example-optimisation-results-perc.svg")
#plt.show()