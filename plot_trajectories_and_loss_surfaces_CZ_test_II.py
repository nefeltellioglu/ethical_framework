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


if len(sys.argv) > 1:
    config_file = sys.argv[1]
else:
    config_file = "config/config-2024-10-14_manuscript.json"
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

output_dir = f"out/CZ_test_II/glam/{config_date_name}"
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


# TODO We should really remove this hard-coded value.
times = 500
if "low" in config_file:
   times = int(times * 20)
elif "high" in config_file:
   times = int(times * 0.5)

# for plotting trajectories 
solutions = {ethical_a_b: em.sir_vacc(params=unique_model_param, 
                          sir_0=initial_conditions[ethical_a_b], 
                          ts=np.linspace(0, times, times + 1))[0]
             for ethical_a_b in ethical_a_b_list}

# ====================================================================
# Plot the optimal trajectories for each ethical configuration so we
# can see what the results look like.
# ====================================================================

print(70 * "=")
print("plotting optimal trajectories for (a, b) examples")
print(70 * "=")

fig, axs = plt.subplots(1, 3, figsize=(10, 2.5))
subplot_labels = ['(a)', '(b)', '(c)']

# TODO These colour codes should not be hard-coded, but they are just
# colorbrewer2 8-class Dark2 values so they shouldn't be too
# mysterious for now.
green_hex = "#1b9e77"
orange_hex = "#d95f02"

for ix, ethical_a_b in enumerate(ethical_a_b_list):
    a, b = ethical_a_b
    ax = axs[ix]
    ax.text(-0.25, 1.15, subplot_labels[ix], transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='right')
    
    sol = solutions[ethical_a_b]

    total_s1 = 100 * (sol.s1 + sol.s1_vp + sol.s1_vu) / pop_size_1
    total_s2 = 100 * (sol.s2 + sol.s2_vp + sol.s2_vu) / pop_size_2
    total_i1 = 100 * (sol.i1 + sol.i1_vu) / pop_size_1
    total_i2 = 100 * (sol.i2 + sol.i2_vu) / pop_size_2

    ax.plot(sol.times, total_s1, color = green_hex, linestyle="solid", label = "Group 1 susceptible")
    ax.plot(sol.times, total_s2, color = green_hex, linestyle="dashed", label = "Group 2 susceptible")
    ax.plot(sol.times, total_i1, color = orange_hex, linestyle="solid", label = "Group 1 infectious")
    ax.plot(sol.times, total_i2, color = orange_hex, linestyle="dashed", label = "Group 2 infectious")

    vacc = [int(100 * x/y) for (x, y) in zip(opt_vacc_strat[ethical_a_b], (pop_size_1, pop_size_2))]

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

fig.savefig(f"{output_dir}/glamorous-trajectories.png", bbox_inches='tight', dpi=300)
fig.savefig(f"{output_dir}/glamorous-trajectories.svg", bbox_inches='tight')


# ====================================================================
# Plot the loss surfaces for each of the different components along
# with their component specific optimal and the global optimal.
# ====================================================================

print(70 * "=")
print("plotting loss surfaces for individual components")
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

# for seaborne heatmap, use a matrix. 

loss_mtx_cb = np.zeros((len(g1_vac_nums), len(g2_vac_nums)))
loss_mtx_ei = np.zeros((len(g1_vac_nums), len(g2_vac_nums)))
loss_mtx_ev = np.zeros((len(g1_vac_nums), len(g2_vac_nums)))
# TODO: these need to be taken from the valid configurations. 
for ix, g1_vac_num in enumerate(g1_vac_nums):
    for jx, g2_vac_num in enumerate(g2_vac_nums):
        ic_key = (g1_vac_num, g2_vac_num)

        if ic_key in ics_objs:

            ic_obj = ics_objs[ic_key]
            ic_id = ics_ids[ic_key]
            cf_id = cfs[ic_id]
            oc = ocs[cf_id]

            l_cb, l_ei, l_ev = em.loss_terms(oc, ic_obj, unique_burden_param)

            loss_mtx_cb[ix, jx] = l_cb
            loss_mtx_ei[ix, jx] = l_ei
            loss_mtx_ev[ix, jx] = l_ev
        else: 
            loss_mtx_cb[ix, jx] = float("nan")
            loss_mtx_ei[ix, jx] = float("nan")
            loss_mtx_ev[ix, jx] = float("nan")


# --------------------------------------------------------------------
# Draw the actual figure
# --------------------------------------------------------------------
x_ann_shift = 5
y_ann_shift = 5

def setup_axes(my_ax, g2_vac_nums, g1_vac_nums):
    my_ax.set_xticks(range(len(g2_vac_nums)), labels=g2_vac_nums, rotation=45)
    my_ax.set_xlabel("Group 2 vaccinations")
    y_ticks_thinned = range(0, len(g1_vac_nums), 7)
    y_labels_thinned = [g1_vac_nums[i] for i in y_ticks_thinned]
    my_ax.set_yticks(y_ticks_thinned, labels=y_labels_thinned)
    my_ax.set_ylabel("Group 1 vaccinations")
    my_ax.set_aspect(len(g2_vac_nums) / len(g1_vac_nums))
    #my_ax.figure.colorbar(im, ax=my_ax)

def annotate_vacc_opt_choice(my_ax, opt_vacc_strat):
    for (name, (g1_v, g2_v)) in zip(["A", "B", "C"], opt_vacc_strat.values()):
        x_index = g2_vac_nums.index(g2_v)
        y_index = g1_vac_nums.index(g1_v)
        my_ax.annotate(name + " optimal", (x_index - x_ann_shift, y_index - y_ann_shift), color="red")
        my_ax.scatter(x_index, y_index, color="red", s=100, marker="o")

# ....................................................................
fig, ax = plt.subplots(1, 3, figsize=(15, 9))
ax_cb = ax[0]
ax_ei = ax[1]
ax_ev = ax[2]
# ....................................................................
#im = ax_cb.imshow(loss_mtx_cb, cmap="viridis_r", aspect="auto", origin="lower")
im = sns.heatmap(loss_mtx_cb, cmap = "viridis_r",  ax=ax_cb, cbar_kws=dict(location='bottom'))
ax_cb.set_title("Total clinical burden", fontweight="bold")
setup_axes(ax_cb, g2_vac_nums, g1_vac_nums)
annotate_vacc_opt_choice(ax_cb, opt_vacc_strat)
#annotate_global_opt(ax_cb, loss_mtx_cb)
# ....................................................................
#im = ax_ei.imshow(loss_mtx_ei, cmap="viridis_r", aspect="auto", origin="lower")
im = sns.heatmap(loss_mtx_ei, cmap = "viridis_r",  ax=ax_ei, cbar_kws=dict(location='bottom'))
ax_ei.set_title("Inequity of infection burden", fontweight="bold")
setup_axes(ax_ei, g2_vac_nums, g1_vac_nums)
annotate_vacc_opt_choice(ax_ei, opt_vacc_strat)
#annotate_global_opt(ax_ei, loss_mtx_ei)
# ....................................................................
#im = ax_ev.imshow(loss_mtx_ev, cmap="viridis_r", aspect="auto", origin="lower")
im = sns.heatmap(loss_mtx_ev, cmap = "viridis_r",  ax=ax_ev, cbar_kws=dict(location='bottom'))
ax_ev.set_title("Inequity of vaccination burden", fontweight="bold")
setup_axes(ax_ev, g2_vac_nums, g1_vac_nums)
annotate_vacc_opt_choice(ax_ev, opt_vacc_strat)
#annotate_global_opt(ax_ev, loss_mtx_ev)
# ....................................................................
fig.tight_layout()
fig.savefig(f"{output_dir}/glamorous-loss_surfaces.png", bbox_inches='tight', dpi=300)
fig.savefig(f"{output_dir}/glamorous-loss_surfaces.svg", bbox_inches='tight')
# --------------------------------------------------------------------


# ====================================================================
# Plot aggregate loss surface for each (a, b) combination
# with the (a,b) specific optimal (v1, v2) combination annotated
# ====================================================================

print(70 * "=")
print("plotting global loss surfaces")
print(70 * "=")


# use the previously-constructed mapping from the initial conditions to
# the outcomes given these initial conditions so that it is easy to
# get the data into a plottable format.

# for seaborne heatmap, use a matrix. 

loss_mtxs_ab = {}

for i_ab, (a, b) in enumerate(ethical_a_b_list):

    loss_mtxs_ab[(a, b)] = np.zeros((len(g1_vac_nums), len(g2_vac_nums)))

    extreme_burdens = eo.get_extreme_burdens(0, 0, db)

    for ix, g1_vac_num in enumerate(g1_vac_nums):
        for jx, g2_vac_num in enumerate(g2_vac_nums):
            ic_key = (g1_vac_num, g2_vac_num)

            if ic_key in ics_objs:
                ic_obj = ics_objs[ic_key]
                ic_id = ics_ids[ic_key]
                cf_id = cfs[ic_id]
                oc = ocs[cf_id]

                l_cb, l_ei, l_ev, l_1, l_2, l_3 = eo.normalisation(a, b, oc, ic_obj, unique_burden_param, extreme_burdens)

                loss = em.global_loss(l_cb, l_ei, l_ev, a, b)
                #loss = (1 - a - b) * l_cb + a * l_ei + b * l_ev 

                loss_mtxs_ab[(a, b)][ix, jx] = loss 
            else:
                loss_mtxs_ab[(a, b)][ix, jx] = float("nan") 



# TODO: matrices of global loss terms over (v1, v2) for each (a, b) combo. 

# --------------------------------------------------------------------
# Draw the actual figure
# --------------------------------------------------------------------
x_ann_shift = 5
y_ann_shift = 5



def annotate_global_opt(my_ax, loss_mtx):
    min_idx = np.argmin(loss_mtx)
    min_idx_g1, min_idx_g2 = np.unravel_index(min_idx, loss_mtx_cb.shape)
    print(f"Global optimal at ({g2_vac_nums[min_idx_g2]}, {g1_vac_nums[min_idx_g1]})")
    tmp = min_idx_g1 - 2 * y_ann_shift
    y_ann_ix = tmp if tmp > 0 else min_idx_g1 + 2 * y_ann_shift
    tmp = min_idx_g2 - x_ann_shift
    x_ann_ix = tmp if tmp > 0 else min_idx_g2
    #my_ax.annotate("Global optimal", (x_ann_ix, y_ann_ix), color="blue")
    my_ax.scatter(min_idx_g2, min_idx_g1, color="blue", s=100, marker="o")



# ....................................................................
fig, ax = plt.subplots(1, 3, figsize=(15, 9))
ax_cb = ax[0]
ax_ei = ax[1]
ax_ev = ax[2]
# ....................................................................
#im = ax_cb.imshow(loss_mtx_cb, cmap="viridis_r", aspect="auto", origin="lower")
loss_mtx = loss_mtxs_ab[ethical_a_b_list[0]]
im = sns.heatmap(loss_mtx, cmap = "viridis_r",  ax=ax_cb, cbar_kws=dict(location='bottom'))
ax_cb.set_title("L", fontweight="bold")
setup_axes(ax_cb, g2_vac_nums, g1_vac_nums)
annotate_global_opt(ax_cb, loss_mtx)
# ....................................................................
#im = ax_ei.imshow(loss_mtx_ei, cmap="viridis_r", aspect="auto", origin="lower")
loss_mtx = loss_mtxs_ab[ethical_a_b_list[1]]
im = sns.heatmap(loss_mtx, cmap = "viridis_r",  ax=ax_ei, cbar_kws=dict(location='bottom'))
ax_ei.set_title("L", fontweight="bold")
setup_axes(ax_ei, g2_vac_nums, g1_vac_nums)
annotate_global_opt(ax_ei, loss_mtx)
# ....................................................................
#im = ax_ev.imshow(loss_mtx_ev, cmap="viridis_r", aspect="auto", origin="lower")
loss_mtx = loss_mtxs_ab[ethical_a_b_list[2]]
im = sns.heatmap(loss_mtx, cmap = "viridis_r",  ax=ax_ev, cbar_kws=dict(location='bottom'))
ax_ev.set_title("L", fontweight="bold")
setup_axes(ax_ev, g2_vac_nums, g1_vac_nums)
annotate_global_opt(ax_ev, loss_mtx)
# ....................................................................
fig.tight_layout()
fig.savefig(f"{output_dir}/glamorous-loss_surfaces_global.png", bbox_inches='tight', dpi=300)
fig.savefig(f"{output_dir}/glamorous-loss_surfaces_global.svg", bbox_inches='tight')
# --------------------------------------------------------------------


