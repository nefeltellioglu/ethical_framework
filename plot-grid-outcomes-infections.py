import json
import pickle
import matplotlib.pyplot as plt

from ethics.model import (
    OptParams,
    BurdenParams,
    SIRParams,
    SIRInitialConditions,
    SIRSolution,
    SIROutcome,
    optimal_initial_conditions,
    loss_clinical_burden,
    loss_equity_of_burden,
    loss_equity_of_vaccination,
    sir_vacc,
    sir_vacc_SSA,
)

with open("config/config-2024-08-07.json", "r") as f:
    CONFIG = json.load(f)


input_file = CONFIG["database_file"]

with open(input_file, "rb") as f:
    db = pickle.load(f)

assert len(db["model_parameters"]) == 1

configs = db["configurations"]
m_params = db["model_parameters"]
ics = db["initial_conditions"]

tmp = {"total_vacc_1": [],
       "total_vacc_2": [],
       "total_inf_1": [],
       "total_inf_2": []}
for oc in db["outcomes"]:
    oc_obj = oc["outcome"]
    tmp["total_vacc_1"].append(oc_obj.total_vac_1)
    tmp["total_vacc_2"].append(oc_obj.total_vac_2)
    tmp["total_inf_1"].append(oc_obj.inf_1_no_vac + oc_obj.inf_1_vac_unprtct + oc_obj.inf_1_vac_prtct)
    tmp["total_inf_2"].append(oc_obj.inf_2_no_vac + oc_obj.inf_2_vac_unprtct + oc_obj.inf_2_vac_prtct)


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
