import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
from dataclasses import dataclass
import numpy as np
import pandas as pd

from ethics.model import OptParams, BurdenParams, SIRParams, SIRInitialConditions, SIRSolution, optimal_initial_conditions, loss_clinical_burden, loss_equity_of_burden, loss_equity_of_vaccination, sir_vacc, sir_vacc_SSA



# ================================

# TODO Params still need to be fixed to include the values from conmat.

#R0 for Delta strain [13]
R0 = 2.7

#Vacc efficacy values [17]
##########
#TODO: Check these assumptions:
#Hospitilised - adverse - given vaccination: days_hosp_vacc * perc_hosp_vacc
#Hospitilised - infection - no vaccination: days_hosp_inf * perc_hosp_inf
#Hospitilised - infection - given vaccination: days_hosp_inf * perc_hosp_inf * (1 - vacc_protection_dis)
#########

vacc_protection_inf = 0.531
vacc_protection_dis = 0.825

#TODO: add references and double check rest of the burden params
disease_burden_params = BurdenParams(
                        perc_hosp_inf= 0.02,
                        days_hosp_inf_1= 3.075,
                        days_hosp_inf_2= 7.60,
                        perc_hosp_vacc_1= 0.004,
                        perc_hosp_vacc_2= 0.002,
                        days_hosp_vacc_1= 6.0,
                        days_hosp_vacc_2= 6.0,
                        vacc_protection_from_disease_1 = vacc_protection_dis,
                        vacc_protection_from_disease_2 = vacc_protection_dis)



#recovery rate delta ~8 days UK[15]
gamma = 0.125 

######
#TODO:Make sure that Bij and contact ij corresponds to the same values!!!
#In the model Bij means "i infects j"    
#Infected i group infects Susceptible j group
#TODO: contact matrix needs to be updated
######
contact_per_capita = np.array([0.40, 0.35, 0.35, 0.30])



#given a contact between I and S, multiply beta by 
#the infectiousness of that I group
#For simplicity, currently consider infectiousness mult as 1
infectiousness_multipliers = np.array([1, 1, 1, 1])



beta_multipliers = contact_per_capita * infectiousness_multipliers

#calculation of beta from R0 and Beta multipliers
beta = R0 * 2 * gamma / (beta_multipliers[0] + beta_multipliers[3] +
                         (beta_multipliers[0]**2 
                          - 2 * beta_multipliers[3] * beta_multipliers[0]
                          + beta_multipliers[3] ** 2
                          + 4 * beta_multipliers[1] * beta_multipliers[2]
                          )**(0.5))

params = SIRParams(beta * beta_multipliers[0], 
                   beta * beta_multipliers[1],
                   beta * beta_multipliers[2],
                   beta * beta_multipliers[3],
                   gamma)

opt_params = OptParams(
    model_type = "Tau-Hybrid",
    #model_type = "SSA",
    no_runs = 100, 
    initial_vacc_1 = 0.5,
    initial_vacc_2 = 0.5,
    stat_type = "mean")

opt_params = OptParams(
    model_type = "ODE",
    no_runs = 1, 
    initial_vacc_1 = 0.5,
    initial_vacc_2 = 0.5,
    stat_type = "mean")
    
#population sizes and sim time horizon
pop_size_1 = 890
pop_size_2 = 110
ts = np.arange(0, 100, 1 / 24)

# tmp_ic = optimal_initial_conditions(params, pop_size_1, pop_size_2, 0, 0)
# tmp_sol = sir_vacc(params, tmp_ic["opt_init_cond"], ts)

foo = []
# Iterate from 0.1 up to 1 in steps of 0.1 then do the same with b conditional on a
for a in np.arange(0.1, 1, 0.1):
    for b in np.arange(0.1, 1 - a, 0.1):
        print("Running a=%s, b=%s"%(a, b))
        tmp_ic = optimal_initial_conditions(params, disease_burden_params,
                                  opt_params, ts,
                                            pop_size_1, pop_size_2,
                                            vacc_protection_inf, vacc_protection_inf,
                                            a, b)
        if opt_params.model_type == "ODE":
            tmp_sol = sir_vacc(params,  
                               tmp_ic["opt_init_cond"], ts)[0]
            
            foo.append(
                {
                    "a": a,
                    "b": b,
                    "total_infections_1": tmp_sol.total_infections()["inf_in_1"],
                    "total_infections_2": tmp_sol.total_infections()["inf_in_2"],
                    "total_vaccinated_1": tmp_sol.total_vaccinated()["vacc_1"],
                    "total_vaccinated_2": tmp_sol.total_vaccinated()["vacc_2"],
                    "loss_clinical_burden": loss_clinical_burden([tmp_sol], 
                                             disease_burden_params)[0],
                    "loss_equity_of_burden": loss_equity_of_burden([tmp_sol], 
                                              disease_burden_params)[0],
                    "loss_equity_of_vaccination": loss_equity_of_vaccination(
                                                [tmp_sol], disease_burden_params)[0],
                }
            )
            
        elif opt_params.model_type in ["SSA", "Tau-Hybrid"]:
            
            tmp_sols = sir_vacc_SSA(params,  
                               tmp_ic["opt_init_cond"], opt_params, ts)
            
            total_infections_1 = []
            total_infections_2 = []
            total_vaccinated_1 = []
            total_vaccinated_2 = []
            for tmp_sol in tmp_sols:
                total_infections_1.append(tmp_sol.total_infections()["inf_in_1"])
                total_infections_2.append(tmp_sol.total_infections()["inf_in_2"])
                total_vaccinated_1.append(tmp_sol.total_vaccinated()["vacc_1"])
                total_vaccinated_2.append(tmp_sol.total_vaccinated()["vacc_2"])
                
            if opt_params.stat_type == "mean":
                
                foo.append(
                    {
                        "a": a,
                        "b": b,
                        "total_infections_1": np.mean(total_infections_1),
                        "total_infections_2": np.mean(total_infections_2),
                        "total_vaccinated_1": np.mean(total_vaccinated_1),
                        "total_vaccinated_2": np.mean(total_vaccinated_2),
                        "loss_clinical_burden": np.mean(loss_clinical_burden(tmp_sols, 
                                                     disease_burden_params)),
                        "loss_equity_of_burden": np.mean(loss_equity_of_burden(tmp_sols, 
                                                     disease_burden_params)),
                        "loss_equity_of_vaccination": np.mean(loss_equity_of_vaccination(
                                                    tmp_sols, disease_burden_params)),
                    }
                )
            elif opt_params.stat_type == "median":
                
                foo.append(
                    {
                        "a": a,
                        "b": b,
                        "total_infections_1": np.median(total_infections_1),
                        "total_infections_2": np.median(total_infections_2),
                        "total_vaccinated_1": np.median(total_vaccinated_1),
                        "total_vaccinated_2": np.median(total_vaccinated_2),
                        "loss_clinical_burden": np.median(loss_clinical_burden(tmp_sols, 
                                                     disease_burden_params)),
                        "loss_equity_of_burden": np.median(loss_equity_of_burden(tmp_sols, 
                                                     disease_burden_params)),
                        "loss_equity_of_vaccination": np.median(loss_equity_of_vaccination(
                                                    tmp_sols, disease_burden_params)),
                    }
                )
                


# Saving a copy of the results in this way is useful as a way to track
# that we haven't broken anything in the code in the future.
df = pd.DataFrame(foo)
df.to_csv("scratch-fancy-%s.csv"%opt_params.model_type)

# Do a matplotlib plot.
# Do a heatmap where the `a` column is on the x-axis and the `b` column is on the y-axis and the `loss_clinical_burden` is the colour.
s_size = 200
fig_demo_1_heatmap = "scratch-fancy-%s-heatmap-CB.png"%opt_params.model_type
plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["loss_clinical_burden"], cmap="viridis", s=s_size)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Clinical Burden Loss")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_1_heatmap)



# Do a heatmap where the `a` column is on the x-axis and the `b` column is on the y-axis and the `loss_equity_of_burden` is the colour.

fig_demo_2_heatmap = "scratch-fancy-%s-heatmap-EB.png"%opt_params.model_type
plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["loss_equity_of_burden"], cmap="viridis", s=s_size)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Equity of Burden Loss")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_2_heatmap)

# Do a heatmap where the `a` column is on the x-axis and the `b` column is on the y-axis and the `loss_equity_of_vaccination` is the colour.

fig_demo_3_heatmap = "scratch-fancy-%s-heatmap-EV.png"%opt_params.model_type
plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["loss_equity_of_vaccination"], cmap="viridis", s=s_size)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Equity of Vaccination Loss")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_3_heatmap)

# Do a heatmap where the `a` column is on the x-axis and the `b` column is on the y-axis and the `loss_clinical_burden + loss_equity_of_burden + loss_equity_of_vaccination` is the colour.

fig_demo_4_heatmap = "scratch-fancy-%s-heatmap-ALL.png"%opt_params.model_type
plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["loss_clinical_burden"] + df["loss_equity_of_burden"] + df["loss_equity_of_vaccination"], cmap="viridis", s=s_size)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Combined Loss")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_4_heatmap)

# Do a heatmap where the `a` column is on the x-axis and the `b` column is on the y-axis and the `total_vaccine_<X>` is the colour.

fig_demo_5_1_heatmap = "scratch-fancy-%s-heatmap-TV1.png"%opt_params.model_type
fig_demo_5_2_heatmap = "scratch-fancy-%s-heatmap-TV2.png"%opt_params.model_type

plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["total_vaccinated_1"], cmap="viridis", s=s_size)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Total Vaccinations 1")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_5_1_heatmap)

plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["total_vaccinated_2"], cmap="viridis", s=s_size)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Total Vaccinations 2")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_5_2_heatmap)

# Do a heatmap where the `a` column is on the x-axis and the `b` column is on the y-axis and the `total_infections_<X>` is the colour.

fig_demo_6_1_heatmap = "scratch-fancy-%s-heatmap-TI1.png"%opt_params.model_type
fig_demo_6_2_heatmap = "scratch-fancy-%s-heatmap-TI2.png"%opt_params.model_type

plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["total_infections_1"], cmap="viridis", s=s_size)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Total Infections 1")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_6_1_heatmap)

plt.figure(figsize=(12, 8))
plt.scatter(df["a"], df["b"], c=df["total_infections_2"], cmap="viridis", s=s_size)
plt.xlabel("a")
plt.ylabel("b")
plt.title("Total Infections 2")
plt.colorbar()
plt.grid()
plt.savefig(fig_demo_6_2_heatmap)
