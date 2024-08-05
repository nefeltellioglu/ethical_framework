import pickle
import matplotlib.pyplot as plt


input_file = "out/grid_database.pkl"
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
    tmp["total_vacc_1"].append(oc["vac_1"])
    tmp["total_vacc_2"].append(oc["vac_2"])
    tmp["total_inf_1"].append(oc["inf_1_no_vac"] + oc["inf_1_vu"] + oc["inf_1_vp"])
    tmp["total_inf_2"].append(oc["inf_2_no_vac"] + oc["inf_2_vu"] + oc["inf_2_vp"])


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
