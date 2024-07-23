import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
from dataclasses import dataclass
import numpy as np
import pandas as pd
import gillespy2

@dataclass
class SIRParams:
    beta_11: float
    beta_12: float
    beta_21: float
    beta_22: float
    gamma: float
    

@dataclass
class BurdenParams:
    perc_hosp_inf: float
    days_hosp_inf_1: float
    days_hosp_inf_2: float
    perc_hosp_vacc_1: float
    perc_hosp_vacc_2: float
    days_hosp_vacc_1: float
    days_hosp_vacc_2: float
    
@dataclass
class OptParams:
    model_type: str
    no_runs: int 
    initial_vacc_1: float
    initial_vacc_2: float
    stat_type: str
    opt_method: str
    
    


@dataclass
class SIRInitialConditions:
    s0_1: int
    s0_2: int
    i0_1: int
    i0_2: int
    r0_1: int
    r0_2: int


# TODO This needs to include a slot for the cost of an infection.
# TODO This needs to include a slot for the cost of vaccination.
@dataclass
class SIRSolution:
    s1: [float]
    s2: [float]
    i1: [float]
    i2: [float]
    r1: [float]
    r2: [float]
    times: [float]

    def total_vaccinated(self) -> dict:
        return {
            "vacc_1": self.r1[0],
            "vacc_2": self.r2[0],
            "total": self.r1[0] + self.r2[0],
        }

    def total_infections(self) -> dict:
        return {
            "inf_in_1": self.i1[-1] + self.r1[-1] - self.r1[0],
            "inf_in_2": self.i2[-1] + self.r2[-1] - self.r2[0],
            "total": self.i1[-1]
            + self.r1[-1]
            - self.r1[0]
            + self.i2[-1]
            + self.r2[-1]
            - self.r2[0],
        }

    def total_population(self) -> dict:
        return {
            "pop_1": self.s1[0] + self.i1[0] + self.r1[0],
            "pop_2": self.s2[0] + self.i2[0] + self.r2[0],
            "total": self.s1[0]
            + self.i1[0]
            + self.r1[0]
            + self.s2[0]
            + self.i2[0]
            + self.r2[0],
        }


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
def loss_clinical_burden(sir_sols: [SIRSolution],
                         disease_burden_params:BurdenParams) -> float:
    
    loss_clinical_burdens = []
    for sir_sol in sir_sols:
        ttl_infs = sir_sol.total_infections()
        ttl_vacc = sir_sol.total_vaccinated()
        net_days_hosp_for_inf = \
            disease_burden_params.perc_hosp_inf * \
                (disease_burden_params.days_hosp_inf_1 * ttl_infs["inf_in_1"]+
                 disease_burden_params.days_hosp_inf_2 * ttl_infs["inf_in_2"])
                
        net_days_hosp_for_vacc = \
            (disease_burden_params.days_hosp_vacc_1 * ttl_vacc["vacc_1"] * \
             disease_burden_params.perc_hosp_vacc_1) + \
             (disease_burden_params.days_hosp_vacc_2 * ttl_vacc["vacc_2"] * \
              disease_burden_params.perc_hosp_vacc_2)
        
        loss_clinical_burden = net_days_hosp_for_inf + net_days_hosp_for_vacc
        loss_clinical_burdens.append(loss_clinical_burden)
    return np.array(loss_clinical_burdens)
    # return ttl_infs["total"] + 0.5 * ttl_vacc["total"]


# TODO This needs to be updated use the stochastic burden of the cost
# of an infection. Recall the $ C_{I}^{i} $ is random per infection
# (zero-inflated), although this should be calculated in the stochastic
# simulation.
def loss_equity_of_burden(sir_sols: [SIRSolution],
                          disease_burden_params:BurdenParams) -> float:
    # ttl_infs = sir_sol.total_infections()
    # ttl_pop = sir_sol.total_population()
    # exp_burden_1 = ttl_infs["total"] * (ttl_pop["pop_1"] / ttl_pop["total"])
    # exp_burden_2 = ttl_infs["total"] * (ttl_pop["pop_2"] / ttl_pop["total"])
    # obs_burden_1 = ttl_infs["inf_in_1"]
    # obs_burden_2 = ttl_infs["inf_in_2"]
    loss_equity_of_burdens = []
    for sir_sol in sir_sols:
        ttl_infs = sir_sol.total_infections()
        ttl_pop = sir_sol.total_population()
        obs_burden_1 = disease_burden_params.perc_hosp_inf * \
                       disease_burden_params.days_hosp_inf_1 * \
                       ttl_infs["inf_in_1"]
        obs_burden_2 = disease_burden_params.perc_hosp_inf * \
                       disease_burden_params.days_hosp_inf_2 * \
                       ttl_infs["inf_in_2"]
        total_inf_burden = obs_burden_1 + obs_burden_2
        exp_burden_1 = total_inf_burden * (ttl_pop["pop_1"] / ttl_pop["total"])
        exp_burden_2 = total_inf_burden * (ttl_pop["pop_2"] / ttl_pop["total"])
        
        loss_equity_of_burden = abs(exp_burden_1 - obs_burden_1) + abs(exp_burden_2 - obs_burden_2)
        loss_equity_of_burdens.append(loss_equity_of_burden)
    return np.array(loss_equity_of_burdens)


# TODO This needs to make use of the cost of vaccination (which is
# stored in the SIRSolution object). Recall that $ C_{V}^{i} $ is
# random per vaccination cost (zero-inflated) and should be calculated
# in the stochastic simulation.
def loss_equity_of_vaccination(sir_sols: [SIRSolution],
                               disease_burden_params:BurdenParams) -> float:
    # ttl_vacc = sir_sol.total_vaccinated()
    # ttl_pop = sir_sol.total_population()
    # exp_vacc_1 = ttl_vacc["total"] * (ttl_pop["pop_1"] / ttl_pop["total"])
    # exp_vacc_2 = ttl_vacc["total"] * (ttl_pop["pop_2"] / ttl_pop["total"])
    # obs_vacc_1 = ttl_vacc["vacc_1"]
    # obs_vacc_2 = ttl_vacc["vacc_2"]
    loss_equity_of_vaccinations = []
    for sir_sol in sir_sols:
        ttl_vacc = sir_sol.total_vaccinated()
        ttl_pop = sir_sol.total_population()
        obs_vacc_1 = disease_burden_params.perc_hosp_vacc_1 * \
                     disease_burden_params.days_hosp_vacc_1 * \
                     ttl_vacc["vacc_1"] 
        obs_vacc_2 = disease_burden_params.perc_hosp_vacc_2 * \
                     disease_burden_params.days_hosp_vacc_2 * \
                     ttl_vacc["vacc_2"] 
                     
        #obs_vacc_1 = 0.002 * 6.0 * 2 * ttl_vacc["vacc_1"]
        #obs_vacc_2 = 0.002 * 6.0 * ttl_vacc["vacc_2"]
        total_vacc_burden = obs_vacc_1 + obs_vacc_2
        exp_vacc_1 = total_vacc_burden * (ttl_pop["pop_1"] / ttl_pop["total"])
        exp_vacc_2 = total_vacc_burden * (ttl_pop["pop_2"] / ttl_pop["total"])
        
        loss_equity_of_vaccination = abs(exp_vacc_1 - obs_vacc_1) + abs(exp_vacc_2 - obs_vacc_2)
        loss_equity_of_vaccinations.append(loss_equity_of_vaccination)
    return np.array(loss_equity_of_vaccinations)


# TODO This will need to be updated to do a stochastic simulation. The
# GillesPy2 package allows you to switch between stochastic and
# deterministic simulations of the same model:
# https://gillespy2.readthedocs.io/en/latest/tutorials/tut_toggle_switch/tut_toggle_switch.html
# This seems like it would be a sensible replacement for the current
# ODE based approach.
def sir_vacc(params: SIRParams,
             sir_0: SIRInitialConditions, ts) -> [SIRSolution]:
    y0 = [sir_0.s0_1, sir_0.s0_2, sir_0.i0_1, sir_0.i0_2, sir_0.r0_1, sir_0.r0_2]

    pop_size_1 = sir_0.s0_1 + sir_0.i0_1 + sir_0.r0_1
    pop_size_2 = sir_0.s0_2 + sir_0.i0_2 + sir_0.r0_2

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
    return [SIRSolution(
        s1=result[:, 0],
        s2=result[:, 1],
        i1=result[:, 2],
        i2=result[:, 3],
        r1=result[:, 4],
        r2=result[:, 5],
        times=ts,
    )]

def sir_vacc_SSA(params: SIRParams, sir_0: SIRInitialConditions,
                 opt_params: OptParams, ts) -> [SIRSolution]:
    y0 = [sir_0.s0_1, sir_0.s0_2, sir_0.i0_1, sir_0.i0_2, sir_0.r0_1, sir_0.r0_2]

    pop_size_1 = sir_0.s0_1 + sir_0.i0_1 + sir_0.r0_1
    pop_size_2 = sir_0.s0_2 + sir_0.i0_2 + sir_0.r0_2
    
    class SIR2(gillespy2.Model): 
        def __init__(self, y0, params, pop_size_1, pop_size_2,t):
            # First call the gillespy2.Model initializer. 
            gillespy2.Model.__init__(self, name='SIR2')
            # Define compartments.
            s_1 = gillespy2.Species(name='s_1', initial_value=y0[0], mode='continuous')
            s_2 = gillespy2.Species(name='s_2', initial_value=y0[1], mode='continuous')
            i_1 = gillespy2.Species(name='i_1', initial_value=y0[2], mode='continuous')
            i_2 = gillespy2.Species(name='i_2', initial_value=y0[3], mode='continuous')
            r_1 = gillespy2.Species(name='r_1', initial_value=y0[4], mode='continuous')
            r_2 = gillespy2.Species(name='r_2', initial_value=y0[5], mode='continuous')
            
            self.add_species([s_1, s_2, i_1, i_2, r_1, r_2])
            
            #Define parameters
            beta_11 = gillespy2.Parameter(name='beta_11', 
                                         expression=params.beta_11/pop_size_1)
            beta_12 = gillespy2.Parameter(name='beta_12', 
                                          expression=params.beta_12/pop_size_2)
            beta_21 = gillespy2.Parameter(name='beta_21', 
                                          expression=params.beta_21/pop_size_1)
            beta_22 = gillespy2.Parameter(name='beta_22', 
                                          expression=params.beta_22/pop_size_2)
            gamma = gillespy2.Parameter(name='gamma', 
                                          expression=params.gamma)
            
            self.add_parameter([beta_11, beta_12, beta_21, beta_22, gamma])
            
            # Define derivatives.
            inf_11 = gillespy2.Reaction(name="inf_11", 
                                        rate=beta_11, 
                                 reactants={s_1:1,i_1:1}, products={i_1:2}) 
            inf_12 = gillespy2.Reaction(name="inf_12", 
                                        rate=beta_12, 
                                 reactants={s_2:1,i_1:1}, products={i_1:1,i_2:1}) 
            inf_21 = gillespy2.Reaction(name="inf_21", 
                                        rate=beta_21, 
                                 reactants={s_1:1,i_2:1}, products={i_1:1,i_2:1}) 
            inf_22 = gillespy2.Reaction(name="inf_22", 
                                        rate=beta_22, 
                                 reactants={s_2:1,i_2:1}, products={i_2:2}) 
            rec_1 = gillespy2.Reaction(name="rec_1", 
                                        rate=gamma, 
                                 reactants={i_1:1}, products={r_1:1}) 
            rec_2 = gillespy2.Reaction(name="rec_2", 
                                        rate=gamma, 
                                 reactants={i_2:1}, products={r_2:1}) 
            self.add_reaction([inf_11,inf_12,inf_21,inf_22,rec_1, rec_2])
            self.timespan(t)
    
    #gillespy plots
    model = SIR2(y0, params, pop_size_1, pop_size_2, ts)
    #results = model.run(number_of_trajectories=1,
    #                    algorithm="ODE")
    
    results = model.run(number_of_trajectories=opt_params.no_runs,
                        algorithm= opt_params.model_type)#"Tau-Hybrid")
    
    
    
    sir_sols = [] 
    
    for result in results:
        sir_sols.append(SIRSolution(
            s1=result["s_1"],
            s2=result["s_2"],
            i1=result["i_1"],
            i2=result["i_2"],
            r1=result["r_1"],
            r2=result["r_2"],
            times=ts,))
    
    return sir_sols





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
        params: SIRParams, disease_burden_params:BurdenParams, 
        opt_params: OptParams,
        ts, 
        pop_size_1: float, pop_size_2: float, a: float, b: float
) -> float:
    def objective(vacc_props: list) -> float:
        init_cond = initial_cond_from_vacc(
            vacc_props[0], vacc_props[1], pop_size_1, pop_size_2
        )
        if opt_params.model_type == "ODE":
            sir_sol = sir_vacc(params, init_cond, ts)
            return (
                (1 - a - b) * loss_clinical_burden(sir_sol, disease_burden_params)[0]
                + a * loss_equity_of_burden(sir_sol, disease_burden_params)[0]
                + b * loss_equity_of_vaccination(sir_sol, disease_burden_params)[0]
            )
        elif opt_params.model_type == "SSA":
            sir_sols = sir_vacc_SSA(params, init_cond, opt_params, ts)
            
            loss_clinical_burden1 = loss_clinical_burden(sir_sols, disease_burden_params)
            loss_equity_of_burden1 = loss_equity_of_burden(sir_sols, disease_burden_params)
            loss_equity_of_vaccination1 = loss_equity_of_vaccination(sir_sols, disease_burden_params)
            
            objectives = (
                            (1 - a - b) * loss_clinical_burden1
                            + a * loss_equity_of_burden1
                            + b * loss_equity_of_vaccination1
                        )
            if opt_params.stat_type == "mean":
                cur_objective =  np.mean(objectives)
            elif opt_params.stat_type == "median":
                cur_objective =  np.median(objectives)
                
            if not isinstance(cur_objective, float):
                print("here")         
            return cur_objective
        
    return objective

# TODO: SSA model optimization doesnt work with method="Nelder-Mead"
# The optimization for SSA model needs to be fixed
def optimal_initial_conditions(
        params: SIRParams,
        disease_burden_params:BurdenParams,
        opt_params: OptParams,
        ts, pop_size_1: float, 
        pop_size_2: float, a: float, b: float
) -> SIRInitialConditions:
    objective = objective_func_factory(params, disease_burden_params,
                                       opt_params,
                                       ts, pop_size_1, pop_size_2, a, b)
    vacc_upper_bound_1 = 1 - (1 / pop_size_1)
    vacc_upper_bound_2 = 1 - (1 / pop_size_2)
    
    
    opt_result = scipy.optimize.minimize(
        objective,
        [opt_params.initial_vacc_1, opt_params.initial_vacc_2],
        bounds=[(0, vacc_upper_bound_1), (0, vacc_upper_bound_2)],
        method=opt_params.opt_method,
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
