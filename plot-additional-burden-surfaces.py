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
import seaborn as sns

from matplotlib.colors import LogNorm


if len(sys.argv) > 1:
    config_file = sys.argv[1]
else:
    #config_file = "config/config-2024-10-14_manuscript_CZ_test_II.json"
    #config_file = "config/config-2024-10-14_manuscript.json"
     #config_file = "config/config-2024-10-28_limited_vaccine.json"
    # config_file = "config/config-2025-03-12_unlimited_low_R0.json"
    #config_file = "config/config-2024-12-02_unlimited_low_R0.json"
    # config_file = "config/config-2025-03-21_test_unlimited_low_R0.json"
    config_file = "config/config-2025-03-12_unlimited_low_R0.json"
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

# at this point, db contains only a placeholder id and the set of model 
# parameters specified in the config file. 


# ====================================================================
# We only want the simulation results that use an amount of vaccine
# less than the maximum amount specified in the configuration file so
# we filter the results in the database to only include those records.
# If a vaccination limit is not given in the configuration, it assumes
# that there is an unlimited supply.
# ====================================================================

print(70 * "=")
print("filtering vaccination combos")
print(70 * "=")

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

# db now contains the following lists: 
# - model parameters : from the config file
# - initial conditions : for each vaccination strategy in grid search (proportions vaccinated)
# - configurations : specifying the indices of the configurations 
#       which do not exceed the max number of vaccines
# - outcomes specitying : how the numbers of infected individuals from the SIR simulations 
# - burden parameters : a single-element list specifying parameters 
#       used for calculating clinical burdens, given an SIR outcome (prop. hospitalised etc.). 



# --------------------------------------------------------------------
# In the following figure generation we assume that there is a single
# set of model parametrs and burden parameters. We can then avoid a
# lot of computation below by extracting some data based on this.

all_model_param_ids = [p['id'] for p in db['model_parameters']]
assert len(all_model_param_ids) == 1, "There should only be one model parameter set in the database."

unique_model_param_id = all_model_param_ids[0]

unique_model_param = [
    mp for mp in db["model_parameters"] if
    mp["id"] == unique_model_param_id][0]["parameters"]

configs = [c for c in db["configurations"] if c["model_parameters_id"] == unique_model_param_id]

all_burden_param_ids = [p['id'] for p in db['burden_parameters']]
assert len(all_burden_param_ids) == 1, "There should only be one burden parameter set in the database."
unique_burden_param_id = all_burden_param_ids[0]
unique_burden_param = [
    bp for bp in db["burden_parameters"] if
    bp["id"] == unique_burden_param_id][0]["parameters"]
# --------------------------------------------------------------------

# ====================================================================
# Find best the best vaccination strategy for each of the given
# combinations of ethical parameters a and b so that we know which
# trajectories to plot later on.
# ====================================================================

print(70 * "=")
print("locating optimal strategies for (a, b) examples")
print(70 * "=")

ethical_a_b_list = [(0.0, 0.0), (0.99, 0.0), (0.0, 0.99)]

opt_vacc_strat = {}

for (ethical_a, ethical_b) in ethical_a_b_list:

    config_ids = [c["id"] for c in configs]

    ocs = [o for o in db["outcomes"] if o["configuration_id"] in config_ids]

    best_ic_id, _ = eo.optimal_initial_condition(
        ethical_a, ethical_b, unique_model_param_id, unique_burden_param_id, db, normalise=True
    )

    _optimal_config = [c for c in db["configurations"] if c["initial_condition_id"] == best_ic_id]
    
    _optimal_outcome = [
        o for o in db["outcomes"] if o["configuration_id"] == _optimal_config[0]["id"]
    ]

    best_vac_1 = _optimal_outcome[0]["outcome"].total_vac_1

    best_vac_2 = _optimal_outcome[0]["outcome"].total_vac_2

    opt_vacc_strat[(ethical_a, ethical_b)] = (int(best_vac_1), int(best_vac_2))


pop_size_1 = CONFIG["population_parameters"]["pop_size_1"]
pop_size_2 = CONFIG["population_parameters"]["pop_size_2"]
vac_protection_from_inf = CONFIG["vacc_protection_from_infection"]

initial_conditions = {}
for ethical_a_b in ethical_a_b_list:
    num_vac_1, num_vac_2 = opt_vacc_strat[ethical_a_b]
    print(
        f"""
        For ethical (a, b) =  {ethical_a_b}, adding
        initial condition with {num_vac_1} vaccinated in population 1
        and {num_vac_2} vaccinated in population 2.
        """
    )
    s0_1_vp = int(num_vac_1 * vac_protection_from_inf)
    s0_2_vp = int(num_vac_2 * vac_protection_from_inf)
    
    initial_conditions[ethical_a_b] = em.SIRInitialCondition(
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
    )



# ====================================================================
# Plot the loss surfaces for each of the different burden components along
# with their component specific optimal and the global optimal.
# ====================================================================

print(70 * "=")
print("plotting loss surfaces for burden components")
print(70 * "=")


# We start by constructing a mapping from the initial conditions to
# the outcomes given these initial conditions so that it is easy to
# get the data into a plottable format.
ics_ids = {ic["value"].number_vaccinated_by_group(): ic["id"]
           for ic in db["initial_conditions"]}

ics_objs = {ic["value"].number_vaccinated_by_group(): ic["value"]
            for ic in db["initial_conditions"]}

_num_ics_ids = len(ics_ids.values())
ocs = {o["configuration_id"]: o["outcome"] for o in db["outcomes"]}

_num_ocs = len(ocs.values())
cfs = {c["initial_condition_id"]: c["id"] for c in db["configurations"]}

_num_cfs = len(cfs.values())

assert _num_ics_ids == _num_ocs
assert _num_ocs == _num_cfs
assert _num_cfs == _num_ics_ids

g1_vac_nums, g2_vac_nums = zip(*list(ics_ids.keys()))
uniq_sorted = lambda x: sorted(list(set(x)))
g1_vac_nums = uniq_sorted(g1_vac_nums)
g2_vac_nums = uniq_sorted(g2_vac_nums)

# for seaborne heatmap, dataframe
num_rows = len(g1_vac_nums) * len(g2_vac_nums)
v1 = np.zeros((num_rows))
v2 = np.zeros((num_rows))

infections_tot = np.zeros((num_rows))
cb_infections_tot = np.zeros((num_rows))
cb_vaccinations_tot = np.zeros((num_rows))

infections_g1 = np.zeros((num_rows))
cb_infections_g1 = np.zeros((num_rows))
cb_vaccinations_g1 =  np.zeros((num_rows))

infections_g2 = np.zeros((num_rows))
cb_infections_g2 = np.zeros((num_rows))
cb_vaccinations_g2 =  np.zeros((num_rows))

# TODO: these need to be taken from the valid configurations. 
iter = -1

# populations are constant: 
pop_1 = em.pop_1(ics_objs[(0, 0)])
pop_2 = em.pop_2(ics_objs[(0, 0)])

for ix, g1_vac_num in enumerate(g1_vac_nums):
    for jx, g2_vac_num in enumerate(g2_vac_nums):

        iter += 1

        v1[iter] = g1_vac_num / pop_1
        v2[iter] = g2_vac_num / pop_2

        ic_key = (g1_vac_num, g2_vac_num)

        if ic_key in ics_objs:

            ic_obj = ics_objs[ic_key]
            ic_id = ics_ids[ic_key]
            cf_id = cfs[ic_id]
            oc = ocs[cf_id]


            infections_tot[iter] = em.count_infections_group_1(oc) + em.count_infections_group_2(oc)
            cb_infections_tot[iter] = em.total_burden_infections(oc, unique_burden_param)
            cb_vaccinations_tot[iter] = em.total_burden_adverse(oc, unique_burden_param)

            infections_g1[iter] = em.count_infections_group_1(oc)
            cb_infections_g1[iter] = em.total_burden_infections_group_1(oc, unique_burden_param)
            cb_vaccinations_g1[iter] = em.burden_adverse_group_1(oc, unique_burden_param)

            infections_g2[iter] = em.count_infections_group_2(oc)
            cb_infections_g2[iter] = em.total_burden_infections_group_2(oc, unique_burden_param)
            cb_vaccinations_g2[iter] =em.burden_adverse_group_2(oc, unique_burden_param)


        else: 

            infections_tot[iter] = float("nan")
            cb_infections_tot[iter] = float("nan")
            cb_vaccinations_tot[iter] = float("nan")

            infections_g1[iter] = float("nan")
            cb_infections_g1[iter] = float("nan")
            cb_vaccinations_g1[iter] = float("nan")

            infections_g2[iter] = float("nan")
            cb_infections_g2[iter] = float("nan")
            cb_vaccinations_g2[iter] = float("nan")



burden_components = {'v1':v1, 
                     'v2':v2, 
                     'infections_tot':infections_tot, 
                     'cb_infections_tot':cb_infections_tot,
                     'cb_vaccinations_tot':cb_vaccinations_tot,

                     'infections_g1':infections_g1,
                     'cb_infections_g1':cb_infections_g1,
                     'cb_vaccinations_g1':cb_vaccinations_g1,

                     'infections_g2':infections_g2,
                     'cb_infections_g2':cb_infections_g2,
                     'cb_vaccinations_g2':cb_vaccinations_g2}

burden_components_df = pd.DataFrame(burden_components)

mtx_infections_tot = burden_components_df.pivot(index="v1", columns="v2", values = "infections_tot")
mtx_cb_infections_tot = burden_components_df.pivot(index="v1", columns="v2", values = "cb_infections_tot")
mtx_cb_vaccinations_tot = burden_components_df.pivot(index="v1", columns="v2", values = "cb_vaccinations_tot")

mtx_infections_g1 = burden_components_df.pivot(index="v1", columns="v2", values = "infections_g1")
mtx_cb_infections_g1 = burden_components_df.pivot(index="v1", columns="v2", values = "cb_infections_g1")
mtx_cb_vaccinations_g1 = burden_components_df.pivot(index="v1", columns="v2", values = "cb_vaccinations_g1")

mtx_infections_g2 = burden_components_df.pivot(index="v1", columns="v2", values = "infections_g2")
mtx_cb_infections_g2 = burden_components_df.pivot(index="v1", columns="v2", values = "cb_infections_g2")
mtx_cb_vaccinations_g2 = burden_components_df.pivot(index="v1", columns="v2", values = "cb_vaccinations_g2")

# --------------------------------------------------------------------
# Draw the actual figure
# --------------------------------------------------------------------
x_ann_shift = 5
y_ann_shift = 5

def setup_axes(my_ax, g2_vac_nums, g1_vac_nums):
    n_ticks = 5
    x_tick_space = int(round(len(g2_vac_nums)/n_ticks))
    x_ticks_thinned = range(0, len(g2_vac_nums), x_tick_space)

    x_labels_thinned = [round(g2_vac_nums[i]/pop_2 * 100) for i in x_ticks_thinned]

    my_ax.set_xticks(x_ticks_thinned, labels=x_labels_thinned, rotation=45, size = 16)
    my_ax.set_xlabel("% vaccinated (group 2, 70+)",size=20)

    y_tick_space = int(round(len(g1_vac_nums)/n_ticks))
    y_ticks_thinned = range(0, len(g1_vac_nums), y_tick_space)

    y_labels_thinned = [round(g1_vac_nums[i]/pop_1 * 100) for i in y_ticks_thinned] 

    my_ax.set_yticks(y_ticks_thinned, labels=y_labels_thinned, size = 16)
    my_ax.set_ylabel("% vaccinated (group 1, 0-69)", size=20)
    my_ax.set_aspect(len(g2_vac_nums) / len(g1_vac_nums))
    #my_ax.figure.colorbar(im, ax=my_ax)
    #colorbar font size: 
    my_ax.figure.axes[-1].xaxis.label.set_size(25)
    my_ax.figure.axes[-1].tick_params(labelsize = 16)
    my_ax.figure.axes[-1].xaxis.set_label_coords(0.5, 3.0)

def annotate_vacc_opt_choice(my_ax, opt_vacc_strat):
    for (name, (g1_v, g2_v)) in zip(["A", "B", "C"], opt_vacc_strat.values()):
        x_index = g2_vac_nums.index(g2_v)
        y_index = g1_vac_nums.index(g1_v)
        my_ax.annotate(name + " optimal", (x_index - x_ann_shift, y_index - y_ann_shift), color="red")
        my_ax.scatter(x_index, y_index, color="red", s=100, marker="o")

# ....................................................................
# fig, ax = plt.subplots(3, 3, figsize=(15, 20))

# ax_inf_tot = ax[0][0]
# ax_cb_inf_tot = ax[0][1]
# ax_cb_vac_tot = ax[0][2]

# ax_inf_g1 = ax[1][0]
# ax_cb_inf_g1 = ax[1][1]
# ax_cb_vac_g1 = ax[1][2]

# ax_inf_g2 = ax[2][0]
# ax_cb_inf_g2 = ax[2][1]
# ax_cb_vac_g2 = ax[2][2]

#-------------------
fig, ax = plt.subplots(1, 3, figsize=(20, 10))

ax_inf_tot = ax[0]
ax_inf_g1 = ax[1]
ax_inf_g2 = ax[2]


# ....................................................................
im = sns.heatmap(mtx_infections_tot, cmap = "viridis_r",  ax=ax_inf_tot, 
                 cbar_kws=dict(location='bottom'), norm=LogNorm())
cntr = ax_inf_tot.contour(mtx_infections_tot, levels = (10, 100, 1000, 10000), colors = 'black')
cntr_crit = ax_inf_tot.contour(mtx_infections_tot, [15], colors = 'red')
ax_inf_tot.clabel(cntr, cntr.levels)
ax_inf_tot.clabel(cntr_crit, cntr_crit.levels)
im.invert_yaxis()
ax_inf_tot.set_title("Total infections", fontweight="bold")
setup_axes(ax_inf_tot, g2_vac_nums, g1_vac_nums)
annotate_vacc_opt_choice(ax_inf_tot, opt_vacc_strat)
# ....................................................................
# im = sns.heatmap(mtx_cb_infections_tot, cmap = "viridis_r",  ax=ax_cb_inf_tot, 
#                  cbar_kws=dict(location='bottom'), norm=LogNorm())
# im.invert_yaxis()
# ax_cb_inf_tot.set_title("Tot. inf. burden", fontweight="bold")
# setup_axes(ax_cb_inf_tot, g2_vac_nums, g1_vac_nums)
# annotate_vacc_opt_choice(ax_cb_inf_tot, opt_vacc_strat)
# ....................................................................
# im = sns.heatmap(mtx_cb_vaccinations_tot, cmap = "viridis_r",  ax=ax_cb_vac_tot, 
#                  cbar_kws=dict(location='bottom'), norm=LogNorm())
# im.invert_yaxis()
# ax_cb_vac_tot.set_title("Tot. vacc. burden", fontweight="bold")
# setup_axes(ax_cb_vac_tot, g2_vac_nums, g1_vac_nums)
# annotate_vacc_opt_choice(ax_cb_vac_tot, opt_vacc_strat)


# ....................................................................
im = sns.heatmap(mtx_infections_g1, cmap = "viridis_r",  ax=ax_inf_g1, 
                 cbar_kws=dict(location='bottom'), norm=LogNorm())
im.invert_yaxis()
ax_inf_g1.set_title("Grp 1 infections", fontweight="bold")
setup_axes(ax_inf_g1, g2_vac_nums, g1_vac_nums)
annotate_vacc_opt_choice(ax_inf_g1, opt_vacc_strat)
# ....................................................................
# im = sns.heatmap(mtx_cb_infections_g1, cmap = "viridis_r",  ax=ax_cb_inf_g1, 
#                  cbar_kws=dict(location='bottom'), norm=LogNorm())
# im.invert_yaxis()
# ax_cb_inf_g1.set_title("Grp 1 inf. burden", fontweight="bold")
# setup_axes(ax_cb_inf_g1, g2_vac_nums, g1_vac_nums)
# annotate_vacc_opt_choice(ax_cb_inf_g1, opt_vacc_strat)
# ....................................................................
# im = sns.heatmap(mtx_cb_vaccinations_g1, cmap = "viridis_r",  ax=ax_cb_vac_g1, 
#                  cbar_kws=dict(location='bottom'), norm=LogNorm())
# im.invert_yaxis()
# ax_cb_vac_g1.set_title("Grp 1 vacc. burden", fontweight="bold")
# setup_axes(ax_cb_vac_g1, g2_vac_nums, g1_vac_nums)
# annotate_vacc_opt_choice(ax_cb_vac_g1, opt_vacc_strat)


# ....................................................................
im = sns.heatmap(mtx_infections_g2, cmap = "viridis_r",  ax=ax_inf_g2, 
                 cbar_kws=dict(location='bottom'), norm=LogNorm())
im.invert_yaxis()
ax_inf_g2.set_title("Grp 2 infections", fontweight="bold")
setup_axes(ax_inf_g2, g2_vac_nums, g1_vac_nums)
annotate_vacc_opt_choice(ax_inf_g2, opt_vacc_strat)
# ....................................................................
# im = sns.heatmap(mtx_cb_infections_g2, cmap = "viridis_r",  ax=ax_cb_inf_g2, 
#                  cbar_kws=dict(location='bottom'), norm=LogNorm())
# im.invert_yaxis()
# ax_cb_inf_g2.set_title("Grp 2 inf. burden", fontweight="bold")
# setup_axes(ax_cb_inf_g2, g2_vac_nums, g1_vac_nums)
# annotate_vacc_opt_choice(ax_cb_inf_g2, opt_vacc_strat)
# ....................................................................
# im = sns.heatmap(mtx_cb_vaccinations_g2, cmap = "viridis_r",  ax=ax_cb_vac_g2, 
#                  cbar_kws=dict(location='bottom'), norm=LogNorm())
# im.invert_yaxis()
# ax_cb_vac_g2.set_title("Grp 2 vacc. burden", fontweight="bold")
# setup_axes(ax_cb_vac_g2, g2_vac_nums, g1_vac_nums)
# annotate_vacc_opt_choice(ax_cb_vac_g2, opt_vacc_strat)



# ....................................................................
fig.tight_layout()
fig.savefig(f"{output_dir}/glamorous-burden_surfaces.png", bbox_inches='tight', dpi=300)
fig.savefig(f"{output_dir}/glamorous-burden_surfaces.svg", bbox_inches='tight')
# --------------------------------------------------------------------









