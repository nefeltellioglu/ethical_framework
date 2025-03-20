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



p_vac_1 = 0.3
p_vac_2 = 1.0


if len(sys.argv) > 1:
    config_file = sys.argv[1]
else:
    #config_file = "config/config-2024-10-14_manuscript_CZ_test_II.json"
    #config_file = "config/config-2024-10-14_manuscript.json"
     #config_file = "config/config-2024-10-28_limited_vaccine.json"
    config_file = "config/config-2025-03-12_unlimited_low_R0.json"
assert os.path.exists(config_file)

# NOTE This assumes the configuration file is named with the format
# `config-YYYY-MM-DD-<some_name>.json`. The `config_date_name` is used
# as the name for the output directory.
config_date_name = os.path.basename(config_file).replace("config-", "").replace(".json", "")


print(70 * "=")
print("Plotting example results")
print(f"Using config file: {config_file}")
print(70 * "=")
with open(config_file, "r") as f:
    CONFIG = json.load(f)

output_dir = f"out/additional_trajectories/{config_date_name}"
os.makedirs(output_dir, exist_ok=True)

with open(CONFIG["database_file"], "rb") as f:
    db = pickle.load(f)
    assert len(db["model_parameters"]) == 1

# at this point, db contains only a placeholder id and the set of model 
# parameters specified in the config file. 

all_model_param_ids = [p['id'] for p in db['model_parameters']]
assert len(all_model_param_ids) == 1, "There should only be one model parameter set in the database."

unique_model_param_id = all_model_param_ids[0]

unique_model_param = [
    mp for mp in db["model_parameters"] if
    mp["id"] == unique_model_param_id][0]["parameters"]


pop_size_1 = CONFIG["population_parameters"]["pop_size_1"]
pop_size_2 = CONFIG["population_parameters"]["pop_size_2"]
vac_protection_from_inf = CONFIG["vacc_protection_from_infection"]
# --------------------------------------------------------------------

# ====================================================================
# specify an arbitrary vaccination strategy to plot
# ====================================================================
num_vac_1 = p_vac_1 * pop_size_1 
num_vac_2 = p_vac_2 * pop_size_2

v1_label = int(p_vac_1* 100)
v2_label = int(p_vac_2 * 100)

v_label = f'-v1-{v1_label}-v2-{v2_label}'

print(
    f"""
    plotting trajectories with {num_vac_1} vaccinated in population 1
    and {num_vac_2} vaccinated in population 2.
    """
)
s0_1_vp = int(num_vac_1 * vac_protection_from_inf)
s0_2_vp = int(num_vac_2 * vac_protection_from_inf)

initial_conditions = em.SIRInitialCondition(
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
times = 100

# for plotting trajectories 
sol = em.sir_vacc(params=unique_model_param, 
                          sir_0=initial_conditions, 
                          ts=np.linspace(0, times, times + 1))[0]

# TODO These colour codes should not be hard-coded, but they are just
# colorbrewer2 8-class Dark2 values so they shouldn't be too
# mysterious for now.
green_hex = "#1b9e77"
orange_hex = "#d95f02"

#-------------------------------
# plot infections vs time for group 2: 

fig, ax = plt.subplots(1, 1, figsize=(5, 5))

# total_s1 = 100 * (sol.s1 + sol.s1_vp + sol.s1_vu) / pop_size_1
# total_s2 = 100 * (sol.s2 + sol.s2_vp + sol.s2_vu) / pop_size_2
# total_i1 = 100 * (sol.i1 + sol.i1_vu) / pop_size_1
# total_i2 = 100 * (sol.i2 + sol.i2_vu) / pop_size_2

total_s1 = (sol.s1 + sol.s1_vp + sol.s1_vu) 
total_s2 = (sol.s2 + sol.s2_vp + sol.s2_vu)
total_i1 = (sol.i1 + sol.i1_vu)
total_i2 = (sol.i2 + sol.i2_vu)

#ax.plot(sol.times, total_s1, color = green_hex, linestyle="solid", label = "Susceptible (Group 1)")
#ax.plot(sol.times, total_s2, color = green_hex, linestyle="dashed", label = "Susceptible (Group 2)")
#ax.plot(sol.times, total_i1, color = orange_hex, linestyle="solid", label = "Infectious (Group 1)")
ax.plot(sol.times, total_i2, color = orange_hex, linestyle="dashed", label = "Infectious (Group 2)")

textstr = '\n'.join((
    'Vaccination',
f'Group 1: {v1_label}%, Group 2: {v2_label}%',
))

props = {"facecolor":"white", "alpha":1.0}

# ax.text(0.315, 0.5, textstr, transform=ax.transAxes, #fontsize=14,
#     verticalalignment='top', horizontalalignment='left', bbox=props)

ax.set_title(textstr) 

ax.set_xlabel('Day')
ax.set_ylabel('number')

ax.legend().set_visible(True)

plt.subplots_adjust(left=0.1,
            bottom=0.1,
            right=0.9,
            top=0.9,
            wspace=0.4,
            hspace=0.4)

fig.savefig(f"{output_dir}/glamorous-trajectories-I2{v_label}.png", bbox_inches='tight', dpi=300)
#fig.savefig(f"{output_dir}/glamorous-trajectories-I2{v_label}.svg", bbox_inches='tight')

#----------------------------

fig, ax = plt.subplots(1, 1, figsize=(5, 5))

# total_s1 = 100 * (sol.s1 + sol.s1_vp + sol.s1_vu) / pop_size_1
# total_s2 = 100 * (sol.s2 + sol.s2_vp + sol.s2_vu) / pop_size_2
# total_i1 = 100 * (sol.i1 + sol.i1_vu) / pop_size_1
# total_i2 = 100 * (sol.i2 + sol.i2_vu) / pop_size_2

total_s1 = (sol.s1 + sol.s1_vp + sol.s1_vu) 
total_s2 = (sol.s2 + sol.s2_vp + sol.s2_vu)
total_i1 = (sol.i1 + sol.i1_vu)
total_i2 = (sol.i2 + sol.i2_vu)


#ax.plot(sol.times, total_s1, color = green_hex, linestyle="solid", label = "Susceptible (Group 1)")
#ax.plot(sol.times, total_s2, color = green_hex, linestyle="dashed", label = "Susceptible (Group 2)")
ax.plot(sol.times, total_i1, color = orange_hex, linestyle="solid", label = "Infectious (Group 1)")
#ax.plot(sol.times, total_i2, color = orange_hex, linestyle="dashed", label = "Infectious (Group 2)")

textstr = '\n'.join((
    'Vaccination',
f'Group 1: {v1_label}%, Group 2: {v2_label}%',
))

props = {"facecolor":"white", "alpha":1.0}

# ax.text(0.315, 0.5, textstr, transform=ax.transAxes, #fontsize=14,
#     verticalalignment='top', horizontalalignment='left', bbox=props)

ax.set_title(textstr) 
ax.set_xlabel('Day')
ax.set_ylabel('number')

ax.legend().set_visible(True)
plt.subplots_adjust(left=0.1,
            bottom=0.1,
            right=0.9,
            top=0.9,
            wspace=0.4,
            hspace=0.4)

fig.savefig(f"{output_dir}/glamorous-trajectories-I1{v_label}.png", bbox_inches='tight', dpi=300)
#fig.savefig(f"{output_dir}/glamorous-trajectories-I1{v_label}.svg", bbox_inches='tight')

