import json
import pickle
import pandas as pd
import matplotlib.pyplot as plt

import ethics.model as em
import ethics.optimisation as eo

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
            
# At the point where we need to make some plots!
model_param_id = 0
burden_param_id = 0
ethical_a = 0.0
ethical_b = 0.95

bp = [bp for bp in db["burden_parameters"] if bp["id"] == burden_param_id][0][
    "parameters"
]

configs = [
    c for c in db["configurations"] if c["model_parameters_id"] == model_param_id
]
config_ids = [c["id"] for c in configs]
ocs = [o for o in db["outcomes"] if o["configuration_id"] in config_ids]

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

# Use matplotlib to make a plot of the same data.
# Put a single large red dot at the best outcome.

plt.figure()
plt.scatter(plot_df["vac_1"], plot_df["vac_2"], c=plot_df["loss"])
cbar = plt.colorbar()
plt.scatter(best_vac_1, best_vac_2, s=500, c="red", marker="o")
cbar.set_label("Loss")
plt.xlabel("Total Vaccinations in Group 1")
plt.ylabel("Total Vaccinations in Group 2")
plt.savefig("out/example-optimisation-result_a_%s_b_%s.png"%(ethical_a, ethical_b))
plt.savefig("out/example-optimisation-result_a_%s_b_%s.svg"%(ethical_a, ethical_b))
plt.clf()





plt.figure()
plt.scatter(100 * plot_df["vac_1"]/CONFIG["population_parameters"]["pop_size_1"], 
            100 * plot_df["vac_2"]/CONFIG["population_parameters"]["pop_size_2"], 
            c=plot_df["loss"])
cbar = plt.colorbar()
plt.scatter(100 * best_vac_1/CONFIG["population_parameters"]["pop_size_1"], 
            100 * best_vac_2/CONFIG["population_parameters"]["pop_size_2"],
            s=500, c="red", marker="o")
cbar.set_label("Loss")
plt.xlabel("Total Vaccinations in Group 1 (%)")
plt.ylabel("Total Vaccinations in Group 2 (%)")
plt.savefig("out/example-optimisation-result-perc_a_%s_b_%s.png"%(ethical_a, ethical_b))
plt.savefig("out/example-optimisation-result-perc_a_%s_b_%s.svg"%(ethical_a, ethical_b))

plt.clf()


