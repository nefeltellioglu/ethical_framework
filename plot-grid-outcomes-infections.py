import json
import pickle
import matplotlib.pyplot as plt

from ethics.model import (
    OptParams,
    BurdenParams,
    SIRParams,
    SIRSolution,
    SIROutcome,
    optimal_initial_conditions,
    loss_clinical_burden,
    loss_equity_of_burden,
    loss_equity_of_vaccination,
    sir_vacc,
    sir_vacc_SSA,
)

with open(#"config/config-2024-10-14_manuscript.json",
          "config/config-2024-10-28_limited_vaccine.json",
          "r") as f:
    CONFIG = json.load(f)


input_file = CONFIG["database_file"]

with open(input_file, "rb") as f:
    db = pickle.load(f)

assert len(db["model_parameters"]) == 1

if "vaccine_parameters" in CONFIG:
    max_vacc = CONFIG["vaccine_parameters"]['maximum_vacc_rollout']
    
    keys =['initial_conditions', 'configurations', 'outcomes']
    result = []
    #result.append(keys)
    for i in range(len(db[keys[2]])):
        output = [db[x][i] for x in keys]
        result.append(output)
    
    db_limited_vacc = {key: db[key] for key in 
                       ["model_parameters", "burden_parameters"]}
    for tmp_r in result:
        tmp_vacc = tmp_r[0]["value"].s0_1_vp + tmp_r[0]["value"].s0_1_vu + \
                   tmp_r[0]["value"].s0_2_vp + tmp_r[0]["value"].s0_2_vu + \
                   tmp_r[0]["value"].i0_1_vu + tmp_r[0]["value"].i0_2_vu + \
                   tmp_r[0]["value"].r0_1_vu + tmp_r[0]["value"].r0_2_vu 
        if tmp_vacc <= max_vacc:
            if 'initial_conditions' in db_limited_vacc.keys():
              db_limited_vacc['initial_conditions'].append(tmp_r[0])
              db_limited_vacc['configurations'].append(tmp_r[1])
              db_limited_vacc['outcomes'].append(tmp_r[2])
              
            else:
                db_limited_vacc['initial_conditions'] = [tmp_r[0]]
                db_limited_vacc['configurations'] = [tmp_r[1]]
                db_limited_vacc['outcomes'] = [tmp_r[2]]
    db = db_limited_vacc.copy()

configs = db["configurations"]
m_params = db["model_parameters"]
ics = db["initial_conditions"]

tmp = {"total_vacc_1": [], "total_vacc_2": [], "total_inf_1": [], "total_inf_2": []}
for oc in db["outcomes"]:
    oc_obj = oc["outcome"]
    tmp["total_vacc_1"].append(oc_obj.total_vac_1)
    tmp["total_vacc_2"].append(oc_obj.total_vac_2)
    tmp["total_inf_1"].append(oc_obj.inf_1_no_vac + oc_obj.inf_1_vu + oc_obj.inf_1_vp)
    tmp["total_inf_2"].append(oc_obj.inf_2_no_vac + oc_obj.inf_2_vu + oc_obj.inf_2_vp)


plt.figure()
plt.scatter(tmp["total_vacc_1"], tmp["total_inf_1"], c=tmp["total_vacc_2"])
cbar = plt.colorbar()
cbar.set_label("Total Vaccinations in Group 2")
plt.xlabel("Total Vaccinations in Group 1")
plt.ylabel("Total Infections in Group 1")
plt.savefig("out/vacc-vs-inf-group-1.png")
plt.savefig("out/vacc-vs-inf-group-1.svg")
plt.clf()


plt.figure()
plt.scatter(tmp["total_vacc_2"], tmp["total_inf_2"], c=tmp["total_vacc_1"])
cbar = plt.colorbar()
cbar.set_label("Total Vaccinations in Group 1")
plt.xlabel("Total Vaccinations in Group 2")
plt.ylabel("Total Infections in Group 2")
plt.savefig("out/vacc-vs-inf-group-2.png")
plt.savefig("out/vacc-vs-inf-group-2.svg")
plt.clf()


plt.figure()
plt.scatter([100 * i/CONFIG["population_parameters"]["pop_size_1"] for i in tmp["total_vacc_1"]], 
            [100 * i/CONFIG["population_parameters"]["pop_size_1"] for i in tmp["total_inf_1"]], 
            c=[100 * i/CONFIG["population_parameters"]["pop_size_2"] for i in tmp["total_vacc_2"]])
cbar = plt.colorbar()
cbar.set_label("Total Vaccinations in Group 2 (%)")
plt.xlabel("Total Vaccinations in Group 1 (%)")
plt.ylabel("Total Infections in Group 1 (%)")
plt.savefig("out/vacc-vs-inf-group-1-perc.png")
plt.savefig("out/vacc-vs-inf-group-1-perc.svg")
plt.clf()


plt.figure()
plt.scatter([100 * i/CONFIG["population_parameters"]["pop_size_2"] for i in tmp["total_vacc_2"]], 
            [100 * i/CONFIG["population_parameters"]["pop_size_2"] for i in tmp["total_inf_2"]], 
            c=[100 * i/CONFIG["population_parameters"]["pop_size_1"] for i in tmp["total_vacc_1"]])
cbar = plt.colorbar()
cbar.set_label("Total Vaccinations in Group 1 (%)")
plt.xlabel("Total Vaccinations in Group 2 (%)")
plt.ylabel("Total Infections in Group 2 (%)")
plt.savefig("out/vacc-vs-inf-group-2-perc.png")
plt.savefig("out/vacc-vs-inf-group-2-perc.svg")
plt.clf()
