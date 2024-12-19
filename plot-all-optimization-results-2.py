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

if len(sys.argv) > 1:
    config_file = sys.argv[1]
else:
    # config_file = "config/config-2024-10-14_manuscript.json"
    config_file = "config/config-2024-10-28_limited_vaccine.json"
    config_file = "config/config-2024-10-14_manuscript.json"
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

# At the point where we need to make some plots!
model_param_id = 0
burden_param_id = 0
#ethical_a = 0.1
#ethical_b = 0.5

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

        foo, bar = eo.optimal_initial_condition(
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
        if vac_2 >= 1:
            print("here")
        inf_1 = count_infections_group_1(oc_ab) / pop_1(ic_ab)
        inf_2 = count_infections_group_2(oc_ab) / pop_2(ic_ab)
        
        clinical_burden = total_clinical_burden(oc_ab, bp)
        infections_burden = total_burden_infections(oc_ab, bp)
        adverse_burden = total_burden_adverse(oc_ab, bp)
        
        total_vaccination = total_vaccinations(oc_ab)

        p_adverse_1 = adverse_per_capita_1(oc_ab, bp, ic_ab)
        p_adverse_2 = adverse_per_capita_2(oc_ab, bp, ic_ab)

        p_infection_burden_1 = infection_burden_per_capita_1(oc_ab, bp, ic_ab)
        p_infection_burden_2 = infection_burden_per_capita_2(oc_ab, bp, ic_ab)



        # add outcomes of interest to plotting dataframe. 
        plot_df.append(
            {   "a": ethical_a,
                "b": ethical_b,
                #"outcome_id": oc["id"],
                #"config_id": tmp_config_id,
                "vac_1": vac_1,
                "vac_2": vac_2,
                "inf_1": inf_1,
                "inf_2": inf_2,
                "cli_burden": clinical_burden,
                "inf_burden": infections_burden,
                "adv_burden": adverse_burden,
                "total_vacc": total_vaccination,
                "p_adverse_1": p_adverse_1,
                "p_adverse_2": p_adverse_2,
                "p_infburd_1": p_infection_burden_1,
                "p_infburd_2": p_infection_burden_2
                #"loss": _opti mal_outcome,
            }
        )
plot_df = pd.DataFrame(plot_df)




figsize=(10,4)#no in x axis, no in yaxis
labels = ('A', 'B', 'C', 'D', 'E', 'F')
j = 0
fig = plt.figure(figsize=figsize)
variables = ["cli_burden", "inf_burden", "adv_burden"]
titles = ["Total Clinical Burden",
          "Total Infection Burden",
          "Total Adverse Burden"]
colors = ["Reds", "Reds", "Reds"]

for var, title, color in zip(variables, titles, colors):
    no = 131 + j
    ax = fig.add_subplot(no)
    ax.text(-0.3, 1.05, labels[j], transform=ax.transAxes,
      fontsize=12, fontweight='bold', va='top', ha='right')
    j += 1
    
    if var in ["cli_burden", "total_vacc", "inf_burden", "adv_burden"]:
        perc_var = var
    else:
        perc_var = "%s_perc"%var
        plot_df[perc_var] = 100 * plot_df[var]/CONFIG["population_parameters"]["pop_size_%s"%(var.split("_")[1])]
    plot_df["a"] = [round(i,2) for i in plot_df["a"]]
    plot_df["b"] = [round(i,2) for i in plot_df["b"]]
    data = plot_df.pivot(index="b", columns="a", values = perc_var)
    #plt.figure()
    
    ax = sns.heatmap(data, linewidth=0.5,
                     vmin=min(plot_df[perc_var]),
                     vmax=max(plot_df[perc_var]),
                     #yticklabels=['High','Medium','Low','No'],
                     cbar_kws={'label': title,"location":'bottom',
                               "pad":0.25},
                     cmap=color,annot=False, fmt='.0f')
    
    ax.invert_yaxis()
    #plt.xlabel("a")
    #plt.ylabel("b")
    #ax.set_title(titles[j-1])
    ax.set_xlabel("Equity in Infection\nBurden Multiplier (a)")
    ax.set_ylabel("Equity in Vaccination\nBurden Multiplier (b)")
plt.subplots_adjust(left=0.1,
            bottom=0.1, 
            right=0.9, 
            top=0.9, 
            wspace=0.4, 
            hspace=0.6)

fig.savefig(f"{output_dir}/hm_all_burdens.png", bbox_inches='tight', dpi=300)

figsize=(7,8)#no in x axis, no in yaxis
labels = ('A', 'B', 'C', 'D', 'E', 'F')
j = 0
fig = plt.figure(figsize=figsize)
variables = ["inf_1", "inf_2", "vac_1", "vac_2"]
titles = ["Infections in group 1 (%)",
          "Infections in group 2 (%)",
          "Vaccinations in group 1 (%)",
          "Vaccinations in group 2 (%)"]
colors = ["Reds", "Reds", "Purples", "Purples"]
#sns.set(font_scale = 0.9)
for var, title, color in zip(variables, titles, colors):
    no = 221 + j
    ax = fig.add_subplot(no)
    ax.text(-0.3, 1.05, labels[j], transform=ax.transAxes,
      fontsize=12, fontweight='bold', va='top', ha='right')
    j += 1
    
    if var in ["cli_burden", "total_vacc", "inf_burden", "adv_burden"]:
        perc_var = var
    elif var in [ "inf_1", "inf_2","vac_1", "vac_2"]:
        perc_var = "%s_perc"%var
        plot_df[perc_var] = 100 * plot_df[var]
    
    else:
        perc_var = "%s_perc"%var
        plot_df[perc_var] = 100 * plot_df[var]/CONFIG["population_parameters"]["pop_size_%s"%(var.split("_")[1])]
    plot_df["a"] = [round(i,2) for i in plot_df["a"]]
    plot_df["b"] = [round(i,2) for i in plot_df["b"]]
    data = plot_df.pivot(index="b", columns="a", values = perc_var)
    #plt.figure()
    
    ax = sns.heatmap(data, linewidth=0.5,
                     vmin=min(plot_df[perc_var]),
                     vmax=max(plot_df[perc_var]),
                     #yticklabels=['High','Medium','Low','No'],
                     cbar_kws={'label': title,"location":'bottom',
                               "pad":0.25},
                     cmap=color,annot=False, fmt='.0f')
    #plt.gca().collections[0].set_clim(min(plot_df[perc_var]),max(plot_df[perc_var]))
    if var in [ "vac_1", "vac_2"]:
        plt.gca().collections[0].set_clim(0,100)
    #plt.setp(ax._legend.legendHandles, fontsize=16)
    cax = ax.figure.axes[-1]
    cax.tick_params(labelsize=10)
    #if var == "vac_2":
    #    cax.tick_params(labelsize=6)
    #ax.legend(fontsize = 10)
    ax.invert_yaxis()
    #plt.xlabel("a")
    #plt.ylabel("b")
    #ax.set_title(titles[j-1])
    ax.set_xlabel("Equity in Infection\nBurden Multiplier (a)")
    ax.set_ylabel("Equity in Vaccination\nBurden Multiplier (b)")
plt.subplots_adjust(left=0.1,
            bottom=0.1, 
            right=0.9, 
            top=0.9, 
            wspace=0.4, 
            hspace=0.1)

fig.savefig(f"{output_dir}/hm_inf_vacc.png", bbox_inches='tight', dpi=300)


variables = ["inf_1", "inf_2", "vac_1", "vac_2", "cli_burden", "total_vacc"]
labels = ["Infections in group 1 (%)",
          "Infections in group 2 (%)",
          "Vaccinations in group 1 (%)",
          "Vaccinations in group 2 (%)",
          "Total Clinical Burden",
          "Total Number of Vaccinated Individuals"]
colors = ["Reds", "Reds", "Purples", "Purples", "Reds", "Purples"]
for var, label, color in zip(variables, labels, colors):
    if var in ["cli_burden", "total_vacc"]:
        perc_var = var
    elif var in [ "inf_1", "inf_2","vac_1", "vac_2"]:
        perc_var = "%s_perc"%var
        plot_df[perc_var] = 100 * plot_df[var]
    else:
        perc_var = "%s_perc"%var
        plot_df[perc_var] = 100 * plot_df[var]
    plot_df["a"] = [round(i,2) for i in plot_df["a"]]
    plot_df["b"] = [round(i,2) for i in plot_df["b"]]
    data = plot_df.pivot(index="b", columns="a", values = perc_var)
    plt.figure()
    ax = sns.heatmap(data, linewidth=0.5,
                     vmin=min(plot_df[perc_var]),
                     vmax=max(plot_df[perc_var]),
                     #yticklabels=['High','Medium','Low','No'],
                     cbar_kws={'label': label},
                     cmap=color,annot=False, fmt='.0f')
    if var in [ "vac_1", "vac_2"]:
        plt.gca().collections[0].set_clim(0,100)
    ax.invert_yaxis()
    #plt.xlabel("a")
    #plt.ylabel("b")
    plt.xlabel("Equity in Infection Burden Multiplier (a)")
    plt.ylabel("Equity in Vaccination Burden Multiplier (b)")
    if var == "total_vacc":
        plt.savefig(f"{output_dir}/hm_%s_across_all.png"%perc_var, bbox_inches='tight', dpi=300)
        
    #plt.savefig(f"{output_dir}/hm_%s_across_all.png"%perc_var, bbox_inches='tight', dpi=300)
    #plt.savefig(f"{output_dir}/hm_%s_across_all.svg"%perc_var, bbox_inches='tight', dpi=300)

