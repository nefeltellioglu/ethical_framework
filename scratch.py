import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
from dataclasses import dataclass
import numpy as np

from ethical_sir import SIRParams, SIRInitialConditions, SIRSolution


# TODO Make sure that this returns an integer for the initial
# conditions so it works with the stochastic model.
def initial_cond_from_vacc(
    vacc_prop_1: float, vacc_prop_2: float, pop_size_1: int, pop_size_2: int
) -> SIRInitialConditions:
    return SIRInitialConditions(
        pop_size_1 * (1 - vacc_prop_1) - 1,
        pop_size_2 * (1 - vacc_prop_2) - 1,
        1,
        1,
        pop_size_1 * vacc_prop_1,
        pop_size_2 * vacc_prop_2,
    )


# TODO This needs to include the cost of vaccination in the same units
# as clinical burden.
def loss_clinical_burden(sir_sol: SIRSolution) -> float:
    ttl_infs = sir_sol.total_infections()
    ttl_vacc = sir_sol.total_vaccinated()
    return ttl_infs["total"] + 0.5 * ttl_vacc["total"]


# TODO This needs to be updated use the stochastic burden of the cost
# of an infection. Recall the $ C_{I}^{i} $ is random per infection
# (zero-inflated), although this should be calculated in the stochastic
# simulation.
def loss_equity_of_burden(sir_sol: SIRSolution) -> float:
    ttl_infs = sir_sol.total_infections()
    ttl_pop = sir_sol.total_population()
    exp_burden_1 = ttl_infs["total"] * (ttl_pop["pop_1"] / ttl_pop["total"])
    exp_burden_2 = ttl_infs["total"] * (ttl_pop["pop_2"] / ttl_pop["total"])
    obs_burden_1 = ttl_infs["inf_in_1"]
    obs_burden_2 = ttl_infs["inf_in_2"]
    return abs(exp_burden_1 - obs_burden_1) + abs(exp_burden_2 - obs_burden_2)


# TODO This needs to make use of the cost of vaccination (which is
# stored in the SIRSolution object). Recall that $ C_{V}^{i} $ is
# random per vaccination cost (zero-inflated) and should be calculated
# in the stochastic simulation.
def loss_equity_of_vaccination(sir_sol: SIRSolution) -> float:
    ttl_vacc = sir_sol.total_vaccinated()
    ttl_pop = sir_sol.total_population()
    exp_vacc_1 = ttl_vacc["total"] * (ttl_pop["pop_1"] / ttl_pop["total"])
    exp_vacc_2 = ttl_vacc["total"] * (ttl_pop["pop_2"] / ttl_pop["total"])
    obs_vacc_1 = ttl_vacc["vacc_1"]
    obs_vacc_2 = ttl_vacc["vacc_2"]
    return abs(exp_vacc_1 - obs_vacc_1) + abs(exp_vacc_2 - obs_vacc_2)


# TODO This will need to be updated to do a stochastic simulation. The
# GillesPy2 package allows you to switch between stochastic and
# deterministic simulations of the same model:
# https://gillespy2.readthedocs.io/en/latest/tutorials/tut_toggle_switch/tut_toggle_switch.html
# This seems like it would be a sensible replacement for the current
# ODE based approach.
def sir_vacc(params: SIRParams, sir_0: SIRInitialConditions, ts) -> SIRSolution:
    y0 = [sir_0.s0_1, sir_0.s0_2, sir_0.i0_1, sir_0.i0_2, sir_0.r0_1, sir_0.r0_2]

    def deriv(y, t, params):
        s_1, s_2, i_1, i_2, r_1, r_2 = y
        inf_11 = params.beta_11 * i_1 * s_1 / pop_size_1
        inf_12 = params.beta_12 * i_1 * s_2 / pop_size_2
        inf_21 = params.beta_21 * i_2 * s_1 / pop_size_1
        inf_22 = params.beta_22 * i_2 * s_2 / pop_size_2
        ds_1 = -inf_11 - inf_21
        ds_2 = -inf_12 - inf_22
        di_1 = inf_11 + inf_21 - params.gamma * i_1
        di_2 = inf_12 + inf_22 - params.gamma * i_2
        dr_1 = params.gamma * i_1
        dr_2 = params.gamma * i_2
        return [ds_1, ds_2, di_1, di_2, dr_1, dr_2]

    result = scipy.integrate.odeint(
        deriv, y0, ts, args=(params,)
    )  # This needs to change probably.
    return SIRSolution(
        s1=result[:, 0],
        s2=result[:, 1],
        i1=result[:, 2],
        i2=result[:, 3],
        r1=result[:, 4],
        r2=result[:, 5],
        times=ts,
    )


# TODO This needs to be extended to include some sort of visualisation
# of the cost associated with infections and vaccnations.
def plot_SIRSolution(sir_sol: SIRSolution) -> None:
    plt.figure(figsize=(12, 8))
    plt.plot(sir_sol.times, sir_sol.s1, label="Susceptible 1")
    plt.plot(sir_sol.times, sir_sol.s2, label="Susceptible 2")
    plt.plot(sir_sol.times, sir_sol.i1, label="Infected 1")
    plt.plot(sir_sol.times, sir_sol.i2, label="Infected 2")
    plt.plot(sir_sol.times, sir_sol.r1, label="Recovered 1")
    plt.plot(sir_sol.times, sir_sol.r2, label="Recovered 2")
    plt.xlabel("Time (months)")
    plt.ylabel("Population")
    plt.legend()
    plt.title("Disease Spread and Vaccination Dynamics")
    plt.grid()
    plt.show()
    return None


# TODO Make a consistent naming between "objective" and "loss".
def objective_func_factory(
    params: SIRParams, pop_size_1: float, pop_size_2: float, a: float, b: float
) -> float:
    def objective(vacc_props: list) -> float:
        init_cond = initial_cond_from_vacc(
            vacc_props[0], vacc_props[1], pop_size_1, pop_size_2
        )
        sir_sol = sir_vacc(params, init_cond, ts)
        return (
            (1 - a - b) * loss_clinical_burden(sir_sol)
            + a * loss_equity_of_burden(sir_sol)
            + b * loss_equity_of_vaccination(sir_sol)
        )

    return objective


def optimal_initial_conditions(
    params: SIRParams, pop_size_1: float, pop_size_2: float, a: float, b: float
) -> SIRInitialConditions:
    objective = objective_func_factory(params, pop_size_1, pop_size_2, a, b)
    vacc_upper_bound_1 = 1 - (1 / pop_size_1)
    vacc_upper_bound_2 = 1 - (1 / pop_size_2)
    opt_result = scipy.optimize.minimize(
        objective,
        [0.5, 0.5],
        bounds=[(0, vacc_upper_bound_1), (0, vacc_upper_bound_2)],
        method="Nelder-Mead",
    )
    if opt_result.success:
        return {
            "opt_init_cond": initial_cond_from_vacc(
                opt_result.x[0], opt_result.x[1], pop_size_1, pop_size_2
            ),
            "obejctive_value": opt_result.fun,
        }
    else:
        raise ValueError("Optimization failed with message: " + opt_result.message)


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
optimal_init_cond = optimal_initial_conditions(params, pop_size_1, pop_size_2, 0, 0)
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
optimal_init_cond = optimal_initial_conditions(params, pop_size_1, pop_size_2, 0.5, 0.0)
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
optimal_init_cond = optimal_initial_conditions(params, pop_size_1, pop_size_2, 0.0, 0.5)
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
optimal_init_cond = optimal_initial_conditions(params, pop_size_1, pop_size_2, 0.2, 0.2)
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
optimal_init_cond = optimal_initial_conditions(params, pop_size_1, pop_size_2, 0.3, 0.3)
result_optimal = sir_vacc(params, optimal_init_cond["opt_init_cond"], ts)
print("Total Infections: " + str(result_optimal.total_infections()))
print("Total Vaccinated: " + str(result_optimal.total_vaccinated()))
print("Loss -- clinical burden: " + str(loss_clinical_burden(result_optimal)))
print("Loss -- equity of burden: " + str(loss_equity_of_burden(result_optimal)))
print(
    "Loss -- equity of vaccination: " + str(loss_equity_of_vaccination(result_optimal))
)
