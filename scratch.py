
import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
from dataclasses import dataclass


@dataclass
class SIRParams:
    beta_11: float
    beta_12: float
    beta_21: float
    beta_22: float
    gamma: float

@dataclass
class SIRInitialConditions:
    s0_1: float
    s0_2: float
    i0_1: float
    i0_2: float
    r0_1: float
    r0_2: float

@dataclass
class SIRSolution:
    s1: float
    s2: float
    i1: float
    i2: float
    r1: float
    r2: float
    times: float



def total_infections(sir_sol: SIRSolution) -> dict:
    ttl_vacc = total_vaccinated(sir_sol)
    return {"inf_in_1": sir_sol.i1[-1] + sir_sol.r1[-1] - ttl_vacc["vacc_1"],
            "inf_in_2": sir_sol.i2[-1] + sir_sol.r2[-1] - ttl_vacc["vacc_2"]}


def total_vaccinated(sir_sol: SIRSolution) -> dict:
    return {"vacc_1": sir_sol.r1[0],
            "vacc_2": sir_sol.r2[0]}


def initial_cond_from_vacc(vacc_prop_1: float,
                           vacc_prop_2: float,
                           pop_size_1: float,
                           pop_size_2: float) -> SIRInitialConditions:
    return SIRInitialConditions(
        pop_size_1 * (1 - vacc_prop_1) - 1,
        pop_size_2 * (1 - vacc_prop_2) - 1,
        1,
        1,
        pop_size_1 * vacc_prop_1,
        pop_size_2 * vacc_prop_2)


def clinical_burden(sir_sol: SIRSolution) -> float:
    ttl_infs = total_infections(sir_sol)
    return ttl_infs["inf_in_1"] + ttl_infs["inf_in_2"]


def vaccinated_propotions(sir_sol: SIRSolution) -> dict:
    ttl_vacc = total_vaccinated(sir_sol)
    ttl_pop = {"pop_1": sir_sol.s1[0] + sir_sol.i1[0] + sir_sol.r1[0],
               "pop_2": sir_sol.s2[0] + sir_sol.i2[0] + sir_sol.r2[0]}
    return {"vacc_prop_1": ttl_vacc["vacc_1"] / ttl_pop["pop_1"],
            "vacc_prop_2": ttl_vacc["vacc_2"] / ttl_pop["pop_2"]}




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

    result = scipy.integrate.odeint(deriv, y0, ts, args=(params,))
    return SIRSolution(s1=result[:, 0], s2=result[:, 1], i1=result[:, 2], i2=result[:, 3], r1=result[:, 4], r2=result[:, 5], times=ts)


def plot_SIRSolution(sir_sol: SIRSolution) -> None:
    plt.figure(figsize=(12, 8))
    plt.plot(sir_sol.times, sir_sol.s1, label='Susceptible 1')
    plt.plot(sir_sol.times, sir_sol.s2, label='Susceptible 2')
    plt.plot(sir_sol.times, sir_sol.i1, label='Infected 1')
    plt.plot(sir_sol.times, sir_sol.i2, label='Infected 2')
    plt.plot(sir_sol.times, sir_sol.r1, label='Recovered 1')
    plt.plot(sir_sol.times, sir_sol.r2, label='Recovered 2')
    plt.xlabel('Time (months)')
    plt.ylabel('Population')
    plt.legend()
    plt.title('Disease Spread and Vaccination Dynamics')
    plt.grid()
    plt.show()
    return None



def objective_func_factory(params: SIRParams,
                           pop_size_1: float,
                           pop_size_2: float) -> float:
    def objective(vacc_props: list) -> float:
        init_cond = initial_cond_from_vacc(vacc_props[0], vacc_props[1], pop_size_1, pop_size_2)
        sir_sol = sir_vacc(params, init_cond, ts)
        return clinical_burden(sir_sol)
    return objective



def optimal_initial_conditions(params: SIRParams,
                               pop_size_1: float,
                               pop_size_2: float) -> SIRInitialConditions:
     objective = objective_func_factory(params, pop_size_1, pop_size_2)
     vacc_upper_bound_1 = 1 - (1 / pop_size_1)
     vacc_upper_bound_2 = 1 - (1 / pop_size_2)
     optimal_vacc_props = scipy.optimize.minimize(
          objective,
          [0.91, 0.91],
          bounds=[(0, vacc_upper_bound_1), (0, vacc_upper_bound_2)])
     return initial_cond_from_vacc(optimal_vacc_props.x[0], optimal_vacc_props.x[1], pop_size_1, pop_size_2)


# ================================

params = SIRParams(0.40, 0.35, 0.35, 0.30, 0.2)
vacc_prop_1 = 0.1
vacc_prop_2 = 0.1
pop_size_1 = 500
pop_size_2 = 500
init_cond = initial_cond_from_vacc(vacc_prop_1, vacc_prop_2, pop_size_1, pop_size_2)
ts = np.arange(0, final_time, time_step)
result = sir_vacc(params, init_cond, ts)
print("======================================================================")
print("Demo Results:")
print("Total Infections: " + str(total_infections(result)))
print("Total Vaccinated: " + str(total_vaccinated(result)))
print("Vaccination Proportions: " + str(vaccinated_propotions(result)))
print("Clinical Burden: " + str(clinical_burden(result)))
print("======================================================================")
print("Results with optimal vaccination proportions for ethics: (1,0,0):")
optimal_init_cond = optimal_initial_conditions(params, pop_size_1, pop_size_2)
result_optimal = sir_vacc(params, optimal_init_cond, ts)
print("Total Infections: " + str(total_infections(result_optimal)))
print("Total Vaccinated: " + str(total_vaccinated(result_optimal)))
print("Vaccination Proportions: " + str(vaccinated_propotions(result_optimal)))
print("Clinical Burden: " + str(clinical_burden(result_optimal)))
