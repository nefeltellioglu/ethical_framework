import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
from dataclasses import dataclass
import numpy as np

from ethical_sir import SIRParams, SIRInitialConditions, SIRSolution, optimal_initial_conditions, loss_clinical_burden, loss_equity_of_burden, loss_equity_of_vaccination, sir_vacc, initial_cond_from_vacc


# ================================

params = SIRParams(0.40, 0.35, 0.35, 0.30, 0.2)
vacc_prop_1 = 0.1
vacc_prop_2 = 0.1
pop_size_1 = 500
pop_size_2 = 500
init_cond = initial_cond_from_vacc(vacc_prop_1, vacc_prop_2, pop_size_1, pop_size_2)
ts = np.arange(0, 100, 1 / 24)
result = sir_vacc(params, init_cond, ts)


print("======================================================================")
print("Demo Results:")
print("Total Infections: " + str(result.total_infections()))
print("Total Vaccinated: " + str(result.total_vaccinated()))
print("Clinical Burden: " + str(loss_clinical_burden(result)))
print("======================================================================")
print("Results with optimal vaccination proportions for ethics: (1,0,0):")
optimal_init_cond = optimal_initial_conditions(params, ts, pop_size_1, pop_size_2, 0, 0)
result_optimal = sir_vacc(params, optimal_init_cond["opt_init_cond"], ts)
print("Objective Value: " + str(optimal_init_cond["obejctive_value"]))
print("Total Infections: " + str(result_optimal.total_infections()))
print("Total Vaccinated: " + str(result_optimal.total_vaccinated()))
print("Loss -- clinical burden: " + str(loss_clinical_burden(result_optimal)))
print("Loss -- equity of burden: " + str(loss_equity_of_burden(result_optimal)))
print(
    "Loss -- equity of vaccination: " + str(loss_equity_of_vaccination(result_optimal))
)
print("======================================================================")
print("Results with optimal vaccination proportions for ethics: (0.5,0.5,0.0):")
optimal_init_cond = optimal_initial_conditions(params, ts, pop_size_1, pop_size_2, 0.5, 0.0)
result_optimal = sir_vacc(params, optimal_init_cond["opt_init_cond"], ts)
print("Total Infections: " + str(result_optimal.total_infections()))
print("Total Vaccinated: " + str(result_optimal.total_vaccinated()))
print("Loss -- clinical burden: " + str(loss_clinical_burden(result_optimal)))
print("Loss -- equity of burden: " + str(loss_equity_of_burden(result_optimal)))
print(
    "Loss -- equity of vaccination: " + str(loss_equity_of_vaccination(result_optimal))
)
print("======================================================================")
print("Results with optimal vaccination proportions for ethics: (0.5,0.0,0.5):")
optimal_init_cond = optimal_initial_conditions(params, ts, pop_size_1, pop_size_2, 0.0, 0.5)
result_optimal = sir_vacc(params, optimal_init_cond["opt_init_cond"], ts)
print("Total Infections: " + str(result_optimal.total_infections()))
print("Total Vaccinated: " + str(result_optimal.total_vaccinated()))
print("Loss -- clinical burden: " + str(loss_clinical_burden(result_optimal)))
print("Loss -- equity of burden: " + str(loss_equity_of_burden(result_optimal)))
print(
    "Loss -- equity of vaccination: " + str(loss_equity_of_vaccination(result_optimal))
)
print("======================================================================")
print("Results with optimal vaccination proportions for ethics: (0.6,0.2,0.2):")
optimal_init_cond = optimal_initial_conditions(params, ts, pop_size_1, pop_size_2, 0.2, 0.2)
result_optimal = sir_vacc(params, optimal_init_cond["opt_init_cond"], ts)
print("Total Infections: " + str(result_optimal.total_infections()))
print("Total Vaccinated: " + str(result_optimal.total_vaccinated()))
print("Loss -- clinical burden: " + str(loss_clinical_burden(result_optimal)))
print("Loss -- equity of burden: " + str(loss_equity_of_burden(result_optimal)))
print(
    "Loss -- equity of vaccination: " + str(loss_equity_of_vaccination(result_optimal))
)
print("======================================================================")
print("Results with optimal vaccination proportions for ethics: (0.4,0.3,0.3):")
optimal_init_cond = optimal_initial_conditions(params, ts, pop_size_1, pop_size_2, 0.3, 0.3)
result_optimal = sir_vacc(params, optimal_init_cond["opt_init_cond"], ts)
print("Total Infections: " + str(result_optimal.total_infections()))
print("Total Vaccinated: " + str(result_optimal.total_vaccinated()))
print("Loss -- clinical burden: " + str(loss_clinical_burden(result_optimal)))
print("Loss -- equity of burden: " + str(loss_equity_of_burden(result_optimal)))
print(
    "Loss -- equity of vaccination: " + str(loss_equity_of_vaccination(result_optimal))
)
