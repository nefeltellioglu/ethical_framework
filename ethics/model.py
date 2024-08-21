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
    # Infection burden parameters
    perc_hosp_inf: float
    days_hosp_inf_1: float
    days_hosp_inf_2: float
    # Protection from disease among (unprotected) vaccinated.
    # TODO Something about this seems a bit fishy to me, maybe we need
    # a clearer variable name. Is it a proportion???
    vacc_protection_from_disease_1: float
    vacc_protection_from_disease_2: float
    # Vaccination burden parameter
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


@dataclass
class SIROutcome:
    inf_1_no_vac: int
    inf_1_vu: int
    inf_1_vp: int
    total_vac_1: int
    inf_2_no_vac: int
    inf_2_vu: int
    inf_2_vp: int
    total_vac_2: int


@dataclass
class SIRInitialCondition:
    s0_1: int
    s0_2: int
    i0_1: int
    i0_2: int
    r0_1: int
    r0_2: int
    s0_1_vp: int
    s0_1_vu: int
    s0_2_vp: int
    s0_2_vu: int
    i0_1_vu: int
    i0_2_vu: int
    r0_1_vu: int
    r0_2_vu: int

    def pop_size(self, pop: int) -> int:
        """
        Returns the total population size for a given population.
        """
        match pop:
            case 1:
                return self.s0_1 + self.i0_1 + self.r0_1 + self.s0_1_vp + self.s0_1_vu + self.i0_1_vu + self.r0_1_vu
            case 2:
                return self.s0_2 + self.i0_2 + self.r0_2 + self.s0_2_vp + self.s0_2_vu + self.i0_2_vu + self.r0_2_vu
            case _:
                raise ValueError(f"Invalid population: {pop}.")

    @staticmethod
    def integer_initial_conditions(vacc_prop_1: float,
                                   vacc_prop_2: float,
                                   pop_size_1: int,
                                   pop_size_2: int,
                                   vacc_protection_1: float,
                                   vacc_protection_2: float):
        """
        Returns initial conditions for the SIR model with integer
        values given a population size and vaccination proportion
        assuming that initially there is a single infected individual
        in each population.

        :param vacc_prop_1: Proportion of population 1 vaccinated.
        :param vacc_prop_2: Proportion of population 2 vaccinated.
        :param pop_size_1: Population size of population 1.
        :param pop_size_2: Population size of population 2.
        :param vacc_protection_1: Proportion of population 1 vaccinated that are protected.
        :param vacc_protection_2: Proportion of population 2 vaccinated that are protected.
        :return: Initial conditions for the SIR model.

        Note: "Protected" means that these individuals cannot get
        infected. This does not yet account for the possibility of an
        imperfect vaccination.

        Note: The integer values are important so that this can be
        used as the initial condition for a CTMC model.
        """
        num_vac_1 = int(pop_size_1 * vacc_prop_1)
        num_vac_1_protected = int(num_vac_1 * vacc_protection_1)
        num_vac_2 = int(pop_size_2 * vacc_prop_2)
        num_vac_2_protected = int(num_vac_2 * vacc_protection_2)
        return SIRInitialCondition(
            pop_size_1 - num_vac_1 - 1,
            pop_size_2 - num_vac_2 - 1,
            1,
            1,
            0,
            0,
            num_vac_1_protected,
            num_vac_1 - num_vac_1_protected,
            num_vac_2_protected,
            num_vac_2 - num_vac_2_protected,
            0,
            0,
            0,
            0
        )


# TODO There is no dynamic vaccination. Double check if that's needed.
# TODO Check if more stats are necessary to be collected from fancy model
# TODO There should be a method to construct `SIRSolution' objects
#      from an initial condition and a parameters object.
@dataclass
class SIRSolution:
    s1: [float]
    s2: [float]
    i1: [float]
    i2: [float]
    r1: [float]
    r2: [float]
    s1_vp: [float]
    s1_vu: [float]
    s2_vp: [float]
    s2_vu: [float]
    i1_vu: [float]
    i2_vu: [float]
    r1_vu: [float]
    r2_vu: [float]
    
    times: [float]

    def total_vaccinated(self) -> dict:
        return {
            "vacc_1": self.s1_vp[0] + self.s1_vu[0],
            "vacc_2": self.s2_vp[0] + self.s2_vu[0],
            "vacc_1_vu": self.s1_vu[0],
            "vacc_2_vu": self.s2_vu[0],
            "total": self.s1_vp[0] + self.s1_vu[0] + \
                     self.s2_vp[0] + self.s2_vu[0],
        }

    def total_infections(self) -> dict:
        
        
                    
        inf_in_1_novacc = self.i1[-1] + self.r1[-1] - self.r1[0] 
        inf_in_2_novacc = self.i2[-1] + self.r2[-1] - self.r2[0]
        
        inf_in_1_vu = self.i1_vu[-1] + self.r1_vu[-1] - self.r1_vu[0]
        inf_in_2_vu = self.i2_vu[-1] + self.r2_vu[-1] - self.r2_vu[0]
        
        inf_in_1 = inf_in_1_novacc + inf_in_1_vu
        
        inf_in_2 = inf_in_2_novacc + inf_in_2_vu
        return {
            "inf_in_1_novacc": inf_in_1_novacc,
            "inf_in_2_novacc": inf_in_2_novacc,
            "inf_in_1_vu": inf_in_1_vu,
            "inf_in_2_vu": inf_in_2_vu,
            "inf_in_1": inf_in_1,
            "inf_in_2": inf_in_2,
            "total": inf_in_1 + inf_in_2,
        }

    def total_population(self) -> dict:
        
        pop_1 = self.s1[0] + self.i1[0] + self.r1[0] + self.s1_vp[0] + \
                self.s1_vu[0] + self.i1_vu[0] + self.r1_vu[0] 
                
        pop_2 = self.s2[0] + self.i2[0] + self.r2[0] + self.s2_vp[0] + \
                self.s2_vu[0] + self.i2_vu[0] + self.r2_vu[0] 
        return {
            "pop_1": pop_1,
            "pop_2": pop_2,
            "total": pop_1 + pop_2,
        }






# TODO This needs to include the cost of vaccination in the same units
# as clinical burden.
# TODO: Nefel introduced eff against severe outcome -> Can someone double check?
def loss_clinical_burden(sir_sols: [SIRSolution],
                         disease_burden_params:BurdenParams) -> float:
    
    loss_clinical_burdens = []
    for sir_sol in sir_sols:
        ttl_infs = sir_sol.total_infections()
        ttl_vacc = sir_sol.total_vaccinated()
        net_days_hosp_for_inf_no_vacc = \
            disease_burden_params.perc_hosp_inf * \
                (disease_burden_params.days_hosp_inf_1 * ttl_infs["inf_in_1_novacc"]+
                 disease_burden_params.days_hosp_inf_2 * ttl_infs["inf_in_2_novacc"])
        
        net_days_hosp_for_inf_vacc = \
            disease_burden_params.perc_hosp_inf * \
                (disease_burden_params.days_hosp_inf_1 * (1 - disease_burden_params.vacc_protection_from_disease_1) * ttl_infs["inf_in_1_vu"] +
                 disease_burden_params.days_hosp_inf_2 * (1 - disease_burden_params.vacc_protection_from_disease_2) * ttl_infs["inf_in_2_vu"])

        net_days_hosp_for_adverse = \
            (disease_burden_params.days_hosp_vacc_1 * ttl_vacc["vacc_1"] * \
             disease_burden_params.perc_hosp_vacc_1) + \
             (disease_burden_params.days_hosp_vacc_2 * ttl_vacc["vacc_2"] * \
              disease_burden_params.perc_hosp_vacc_2)
        
        loss_clinical_burden = net_days_hosp_for_inf_no_vacc + net_days_hosp_for_inf_vacc + net_days_hosp_for_adverse
        loss_clinical_burdens.append(loss_clinical_burden)
    return loss_clinical_burdens
    # return ttl_infs["total"] + 0.5 * ttl_vacc["total"]


# TODO This needs to be updated use the stochastic burden of the cost
# of an infection. Recall the $ C_{I}^{i} $ is random per infection
# (zero-inflated), although this should be calculated in the stochastic
# simulation.
# TODO: Nefel introduced eff against severe outcome -> Can someone double check?
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
        obs_burden_1_novacc = disease_burden_params.perc_hosp_inf * \
                       disease_burden_params.days_hosp_inf_1 * \
                       ttl_infs["inf_in_1_novacc"]
        obs_burden_1_vu = disease_burden_params.perc_hosp_inf * \
                       disease_burden_params.days_hosp_inf_1 * \
                       (1 - disease_burden_params.vacc_protection_from_disease_1 )* \
                       ttl_infs["inf_in_1_vu"]
        obs_burden_1 = obs_burden_1_novacc + obs_burden_1_vu
        
        obs_burden_2_novacc = disease_burden_params.perc_hosp_inf * \
                       disease_burden_params.days_hosp_inf_2 * \
                       ttl_infs["inf_in_2_novacc"]
        obs_burden_2_vu = disease_burden_params.perc_hosp_inf * \
                       disease_burden_params.days_hosp_inf_2 * \
                       (1 - disease_burden_params.vacc_protection_from_disease_2) * \
                       ttl_infs["inf_in_2_vu"]
        obs_burden_2 = obs_burden_2_novacc + obs_burden_2_vu               
                      
        total_inf_burden = obs_burden_1 + obs_burden_2
        exp_burden_1 = total_inf_burden * (ttl_pop["pop_1"] / ttl_pop["total"])
        exp_burden_2 = total_inf_burden * (ttl_pop["pop_2"] / ttl_pop["total"])
        
        loss_equity_of_burden = abs(exp_burden_1 - obs_burden_1) + abs(exp_burden_2 - obs_burden_2)
        loss_equity_of_burdens.append(loss_equity_of_burden)
    return loss_equity_of_burdens


# TODO This needs to make use of the cost of vaccination (which is
# stored in the SIRSolution object). Recall that $ C_{V}^{i} $ is
# random per vaccination cost (zero-inflated) and should be calculated
# in the stochastic simulation.
# TODO: Alex, Cam - No change here when I introduced eff against severe outcome -> Is this correct?
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
    return loss_equity_of_vaccinations


def sir_vacc(params: SIRParams,
             sir_0: SIRInitialCondition, ts) -> [SIRSolution]:
    y0 = [sir_0.s0_1, sir_0.s0_2, sir_0.i0_1, sir_0.i0_2, sir_0.r0_1, sir_0.r0_2,
          sir_0.s0_1_vp, sir_0.s0_1_vu, sir_0.s0_2_vp, sir_0.s0_2_vu, 
          sir_0.i0_1_vu, sir_0.i0_2_vu, sir_0.r0_1_vu, sir_0.r0_2_vu,]


    pop_size_1 = sir_0.s0_1 + sir_0.i0_1 + sir_0.r0_1 + \
                 sir_0.s0_1_vp + sir_0.s0_1_vu + sir_0.i0_1_vu + sir_0.r0_1_vu
    pop_size_2 = sir_0.s0_2 + sir_0.i0_2 + sir_0.r0_2 + \
                 sir_0.s0_2_vp + sir_0.s0_2_vu + sir_0.i0_2_vu + sir_0.r0_2_vu
    
    def deriv(y, t, params):
        s_1, s_2, i_1, i_2, r_1, r_2, s_1_vp, s_1_vu, s_2_vp, s_2_vu, i_1_vu,  i_2_vu, r_1_vu, r_2_vu = y
        
        #unvaccinated  I infect unvaccinated S
        inf_11 = params.beta_11 * i_1 * s_1 / pop_size_1
        inf_12 = params.beta_12 * i_1 * s_2 / pop_size_2
        inf_21 = params.beta_21 * i_2 * s_1 / pop_size_1
        inf_22 = params.beta_22 * i_2 * s_2 / pop_size_2
        
        #vaccinated unprotected I (Ivu) infect unvaccinated S (S)
        inf_1vu1 = params.beta_11 * i_1_vu * s_1 / pop_size_1
        inf_1vu2 = params.beta_12 * i_1_vu * s_2 / pop_size_2
        inf_2vu1 = params.beta_21 * i_2_vu * s_1 / pop_size_1
        inf_2vu2 = params.beta_22 * i_2_vu * s_2 / pop_size_2
        
        #unvaccinated I infect vaccinated unprotected S (Svu)
        inf_11vu = params.beta_11 * i_1 * s_1_vu / pop_size_1
        inf_12vu = params.beta_12 * i_1 * s_2_vu / pop_size_2
        inf_21vu = params.beta_21 * i_2 * s_1_vu / pop_size_1
        inf_22vu = params.beta_22 * i_2 * s_2_vu / pop_size_2
        
        #vaccinated unprotected I (Ivu) infect vaccinated unprotected S (Svu)
        inf_1vu1vu = params.beta_11 * i_1_vu * s_1_vu / pop_size_1
        inf_1vu2vu = params.beta_12 * i_1_vu * s_2_vu / pop_size_2
        inf_2vu1vu = params.beta_21 * i_2_vu * s_1_vu / pop_size_1
        inf_2vu2vu = params.beta_22 * i_2_vu * s_2_vu / pop_size_2
        
        
        ds_1 = -inf_11 - inf_21 - inf_1vu1 - inf_2vu1
        ds_2 = -inf_12 - inf_22 - inf_1vu2 - inf_2vu2
        di_1 = inf_11 + inf_21 + inf_1vu1 + inf_2vu1 - params.gamma * i_1
        di_2 = inf_12 + inf_22 + inf_1vu2 + inf_2vu2 - params.gamma * i_2
        dr_1 = params.gamma * i_1
        dr_2 = params.gamma * i_2
        ds_1_vp = 0
        ds_1_vu = -inf_11vu - inf_21vu - inf_1vu1vu - inf_2vu1vu 
        ds_2_vp = 0
        ds_2_vu = -inf_12vu - inf_22vu - inf_1vu2vu - inf_2vu2vu
        di_1_vu = inf_11vu + inf_21vu + inf_1vu1vu + inf_2vu1vu - params.gamma * i_1_vu  
        di_2_vu = inf_12vu + inf_22vu + inf_1vu2vu + inf_2vu2vu - params.gamma * i_2_vu  
        dr_1_vu = params.gamma * i_1_vu  
        dr_2_vu = params.gamma * i_2_vu  
        
        return [ds_1, ds_2, di_1, di_2, dr_1, dr_2, 
                ds_1_vp ,ds_1_vu, ds_2_vp,  ds_2_vu, di_1_vu, di_2_vu,
                dr_1_vu, dr_2_vu]

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
        
        s1_vp=result[:, 6],
        s1_vu=result[:, 7],
        s2_vp=result[:, 8],
        s2_vu=result[:, 9],
        i1_vu=result[:, 10],
        i2_vu=result[:, 11],
        r1_vu=result[:, 12],
        r2_vu=result[:, 13],
        
        times=ts,
    )]



# TODO: introduce fancy model
def sir_vacc_SSA(params: SIRParams, sir_0: SIRInitialCondition,
                 opt_params: OptParams, ts) -> [SIRSolution]:
    y0 = [sir_0.s0_1, sir_0.s0_2, sir_0.i0_1, sir_0.i0_2, sir_0.r0_1, sir_0.r0_2,
          sir_0.s0_1_vp, sir_0.s0_1_vu, sir_0.s0_2_vp, sir_0.s0_2_vu, 
          sir_0.i0_1_vu, sir_0.i0_2_vu, sir_0.r0_1_vu, sir_0.r0_2_vu,]


    pop_size_1 = sir_0.s0_1 + sir_0.i0_1 + sir_0.r0_1 + \
                 sir_0.s0_1_vp + sir_0.s0_1_vu + sir_0.i0_1_vu + sir_0.r0_1_vu
    pop_size_2 = sir_0.s0_2 + sir_0.i0_2 + sir_0.r0_2 + \
                 sir_0.s0_2_vp + sir_0.s0_2_vu + sir_0.i0_2_vu + sir_0.r0_2_vu
    
    class SIR2(gillespy2.Model): 
        def __init__(self, y0, params, pop_size_1, pop_size_2,t):
            # First call the gillespy2.Model initializer. 
            gillespy2.Model.__init__(self, name='SIR2')
            # Define compartments.
            s_1 = gillespy2.Species(name='s_1', initial_value=y0[0])
            s_2 = gillespy2.Species(name='s_2', initial_value=y0[1])
            i_1 = gillespy2.Species(name='i_1', initial_value=y0[2])
            i_2 = gillespy2.Species(name='i_2', initial_value=y0[3])
            r_1 = gillespy2.Species(name='r_1', initial_value=y0[4])
            r_2 = gillespy2.Species(name='r_2', initial_value=y0[5])
            s_1_vp = gillespy2.Species(name='s_1_vp', initial_value=y0[6])
            s_1_vu = gillespy2.Species(name='s_1_vu', initial_value=y0[7])
            s_2_vp = gillespy2.Species(name='s_2_vp', initial_value=y0[8])
            s_2_vu = gillespy2.Species(name='s_2_vu', initial_value=y0[9])
            i_1_vu = gillespy2.Species(name='i_1_vu', initial_value=y0[10])
            i_2_vu = gillespy2.Species(name='i_2_vu', initial_value=y0[11])
            r_1_vu = gillespy2.Species(name='r_1_vu', initial_value=y0[12])
            r_2_vu = gillespy2.Species(name='r_2_vu', initial_value=y0[13])
            
            
            
            
            self.add_species([s_1, s_2, i_1, i_2, r_1, r_2, 
                              s_1_vp, s_1_vu, s_2_vp, s_2_vu, 
                              i_1_vu,  i_2_vu, r_1_vu, r_2_vu])
            
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
            
            #unvaccinated  I infect unvaccinated S
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
            
            
            #vaccinated unprotected I (Ivu) infect unvaccinated S (S)
            inf_1vu1 = gillespy2.Reaction(name="inf_1vu1", 
                                        rate=beta_11, 
                                 reactants={s_1:1,i_1_vu:1}, products={i_1:1, i_1_vu:1}) 
            inf_1vu2 = gillespy2.Reaction(name="inf_1vu2", 
                                        rate=beta_12, 
                                 reactants={s_2:1,i_1_vu:1}, products={i_2:1, i_1_vu:1}) 
            inf_2vu1 = gillespy2.Reaction(name="inf_2vu1", 
                                        rate=beta_21, 
                                 reactants={s_1:1,i_2_vu:1}, products={i_1:1 ,i_2_vu:1}) 
            inf_2vu2 = gillespy2.Reaction(name="inf_2vu2", 
                                        rate=beta_22, 
                                 reactants={s_2:1,i_2_vu:1}, products={i_2:1,i_2_vu:1}) 
            
            #unvaccinated I infect vaccinated unprotected S (Svu)
            inf_11vu = gillespy2.Reaction(name="inf_11vu", 
                                        rate=beta_11, 
                                 reactants={s_1_vu:1,i_1:1}, products={i_1_vu:1,i_1:1}) 
            inf_12vu = gillespy2.Reaction(name="inf_12vu", 
                                        rate=beta_12, 
                                 reactants={s_2_vu:1,i_1:1}, products={i_2_vu:1,i_1:1}) 
            inf_21vu = gillespy2.Reaction(name="inf_21vu", 
                                        rate=beta_21, 
                                 reactants={s_1_vu:1,i_2:1}, products={i_1_vu:1,i_2:1}) 
            inf_22vu = gillespy2.Reaction(name="inf_22vu", 
                                        rate=beta_22, 
                                 reactants={s_2_vu:1,i_2:1}, products={i_2_vu:1,i_2:1}) 
        
            
            #vaccinated unprotected I (Ivu) infect vaccinated unprotected S (Svu)
            inf_1vu1vu = gillespy2.Reaction(name="inf_1vu1vu", 
                                        rate=beta_11, 
                                 reactants={s_1_vu:1,i_1_vu:1}, products={i_1_vu:2}) 
            inf_1vu2vu = gillespy2.Reaction(name="inf_1vu2vu", 
                                        rate=beta_12, 
                                 reactants={s_2_vu:1,i_1_vu:1}, products={i_2_vu:1,i_1_vu:1}) 
            inf_2vu1vu = gillespy2.Reaction(name="inf_2vu1vu", 
                                        rate=beta_21, 
                                 reactants={s_1_vu:1,i_2_vu:1}, products={i_1_vu:1,i_2_vu:1}) 
            inf_2vu2vu = gillespy2.Reaction(name="inf_2vu2vu", 
                                        rate=beta_22, 
                                 reactants={s_2_vu:1,i_2_vu:1}, products={i_2_vu:2}) 
            
            #recovery from vu
            rec_1_vu = gillespy2.Reaction(name="rec_1_vu", 
                                        rate=gamma, 
                                 reactants={i_1_vu:1}, products={r_1_vu:1}) 
            rec_2_vu = gillespy2.Reaction(name="rec_2_vu", 
                                        rate=gamma, 
                                 reactants={i_2_vu:1}, products={r_2_vu:1}) 
            
            
            self.add_reaction([inf_11,inf_12,inf_21,inf_22,rec_1, rec_2,
                              inf_1vu1, inf_1vu2, inf_2vu1,inf_2vu2,
                              inf_11vu,inf_12vu, inf_21vu, inf_22vu,
                              inf_1vu1vu,inf_1vu2vu, inf_2vu1vu, inf_2vu2vu,
                              rec_1_vu, rec_2_vu
                               ])
            self.timespan(t)
    
    model = SIR2(y0, params, pop_size_1, pop_size_2, ts)
    
    results = model.run(number_of_trajectories=opt_params.no_runs,
                        algorithm= opt_params.model_type)
    
    #sanity check with ODE model results: DONE
    #results = model.run(number_of_trajectories=opt_params.no_runs,
    #                    algorithm= "ODE")#"Tau-Hybrid")
    
    
    
    sir_sols = [] 
    
    for result in results:
        sir_sols.append(SIRSolution(
            s1=result["s_1"],
            s2=result["s_2"],
            i1=result["i_1"],
            i2=result["i_2"],
            r1=result["r_1"],
            r2=result["r_2"],
            
            
            s1_vp=result["s_1_vp"],
            s1_vu=result["s_1_vu"],
            s2_vp=result["s_2_vp"],
            s2_vu=result["s_2_vu"],
            i1_vu=result["i_1_vu"],
            i2_vu=result["i_2_vu"],
            r1_vu=result["r_1_vu"],
            r2_vu=result["r_2_vu"],
            
            times=ts,))
    
    return sir_sols





# TODO This needs to be extended to include some sort of visualisation
# of the cost associated with infections and vaccnations.
# TODO: introduce fancy model
def plot_SIRSolution(sir_sol: SIRSolution) -> None:
    plt.figure(figsize=(12, 8))
    plt.plot(sir_sol.times, sir_sol.s1, label="Susceptible 1")
    plt.plot(sir_sol.times, sir_sol.s2, label="Susceptible 2")
    plt.plot(sir_sol.times, sir_sol.i1, label="Infected 1")
    plt.plot(sir_sol.times, sir_sol.i2, label="Infected 2")
    plt.plot(sir_sol.times, sir_sol.r1, label="Recovered 1")
    plt.plot(sir_sol.times, sir_sol.r2, label="Recovered 2")
    plt.plot(sir_sol.times, sir_sol.s1_vp, label="Susceptible VP 1")
    plt.plot(sir_sol.times, sir_sol.s1_vu, label="Susceptible VU 1")
    plt.plot(sir_sol.times, sir_sol.s2_vp, label="Susceptible VP 2")
    plt.plot(sir_sol.times, sir_sol.s2_vu, label="Susceptible VU 2")
    plt.plot(sir_sol.times, sir_sol.i1_vu, label="Infected VU 1")
    plt.plot(sir_sol.times, sir_sol.i2_vu, label="Infected VU 2")
    plt.plot(sir_sol.times, sir_sol.r1_vu, label="Recovered VU 1")
    plt.plot(sir_sol.times, sir_sol.r2_vu, label="Recovered VU 2")
    
    
    plt.xlabel("Time (days)")
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
        pop_size_1: float, pop_size_2: float,
        vacc_protection_1: float, vacc_protection_2: float,
        a: float, b: float
) -> float:
    def objective(vacc_props: list) -> float:
        init_cond = SIRInitialCondition.integer_initial_conditions(
            vacc_props[0], vacc_props[1], pop_size_1, pop_size_2,
            vacc_protection_1, vacc_protection_2
        )

        if opt_params.model_type == "ODE":
            sir_sol = sir_vacc(params, init_cond, ts)
            #plot_SIRSolution(sir_sol[0])
            
            return (
                (1 - a - b) * loss_clinical_burden(sir_sol, disease_burden_params)[0]
                + a * loss_equity_of_burden(sir_sol, disease_burden_params)[0]
                + b * loss_equity_of_vaccination(sir_sol, disease_burden_params)[0]
            )
        elif opt_params.model_type in ["SSA", "Tau-Hybrid"]:
            sir_sols = sir_vacc_SSA(params, init_cond, opt_params, ts)
            #plot_SIRSolution(sir_sols[0])
            loss_clinical_burden1 = loss_clinical_burden(sir_sols, disease_burden_params)
            loss_equity_of_burden1 = loss_equity_of_burden(sir_sols, disease_burden_params)
            loss_equity_of_vaccination1 = loss_equity_of_vaccination(sir_sols, disease_burden_params)
            
            objectives = (
                            (1 - a - b) * loss_clinical_burden1
                            + a * loss_equity_of_burden1
                            + b * loss_equity_of_vaccination1
                        )
            if opt_params.stat_type == "mean":
                cur_objective =  np.nanmean(objectives)
            elif opt_params.stat_type == "median":
                cur_objective =  np.nanmedian(objectives)
                
            return cur_objective

    return objective

#TODO: introduce grid search for stochastic runs
def optimal_initial_conditions(
        params: SIRParams,
        disease_burden_params:BurdenParams,
        opt_params: OptParams,
        ts,
        pop_size_1: float, pop_size_2: float,
        vacc_protection_1: float, vacc_protection_2: float,
        a: float, b: float
) -> SIRInitialCondition:
    objective = objective_func_factory(params, disease_burden_params,
                                       opt_params,
                                       ts, pop_size_1, pop_size_2,
                                       vacc_protection_1, vacc_protection_2,
                                       a, b)
    vacc_upper_bound_1 = 1 - (1 / pop_size_1)
    vacc_upper_bound_2 = 1 - (1 / pop_size_2)
    opt_result = scipy.optimize.minimize(
        objective,
        [opt_params.initial_vacc_1, opt_params.initial_vacc_2],
        bounds=[(0, vacc_upper_bound_1), (0, vacc_upper_bound_2)],
        method="Nelder-Mead",
    )
    if opt_result.success:
        return {
            "opt_init_cond": SIRInitialCondition.integer_initial_conditions(
                opt_result.x[0], opt_result.x[1], pop_size_1, pop_size_2,
                vacc_protection_1, vacc_protection_2
            ),
            "obejctive_value": opt_result.fun,
        }
    else:
        raise ValueError("Optimization failed with message: " + opt_result.message)


def loss(
    outcome: SIROutcome,
    ic: SIRInitialCondition,
    burden_params: BurdenParams,
    a: float,
    b: float,
) -> float:
    """
    Burden is the sum of infection burden and vaccination burden.

    Args:
    outcome: `SIROutcome`
    ic: `SIRInitialCondition`
    burden_params: `BurdenParams`
    a: The a parameter for the loss function.
    b: The b parameter for the loss function.
    """
    pop_size_1, pop_size_2 = ic.pop_size(1), ic.pop_size(2)
    total_pop = float(pop_size_1 + pop_size_2)

    # import pdb; pdb.set_trace()
    obs_vb_1 = (
        outcome.total_vac_1
        * burden_params.perc_hosp_vacc_1
        * burden_params.days_hosp_vacc_1
    )
    obs_vb_2 = (
        outcome.total_vac_2
        * burden_params.perc_hosp_vacc_2
        * burden_params.days_hosp_vacc_2
    )
    obs_ib_1_no_vac = (
        outcome.inf_1_no_vac
        * burden_params.perc_hosp_inf
        * burden_params.days_hosp_inf_1
    )
    obs_ib_2_no_vac = (
        outcome.inf_2_no_vac
        * burden_params.perc_hosp_inf
        * burden_params.days_hosp_inf_2
    )
    obs_ib_1_vu = (
        outcome.inf_1_vu
        * burden_params.perc_hosp_inf
        * burden_params.days_hosp_inf_1
        * (1 - burden_params.vacc_protection_from_disease_1)
    )
    obs_ib_2_vu = (
        outcome.inf_2_vu
        * burden_params.perc_hosp_inf
        * burden_params.days_hosp_inf_2
        * (1 - burden_params.vacc_protection_from_disease_2)
    )
    obs_cb_1 = obs_ib_1_no_vac + obs_ib_1_vu + obs_vb_1
    obs_cb_2 = obs_ib_2_no_vac + obs_ib_2_vu + obs_vb_2

    exp_cb_1 = (obs_cb_1 + obs_cb_2) * (pop_size_1 / total_pop)
    exp_cb_2 = (obs_cb_1 + obs_cb_2) * (pop_size_2 / total_pop)
    exp_vb_1 = (obs_vb_1 + obs_vb_2) * (pop_size_1 / total_pop)
    exp_vb_2 = (obs_vb_1 + obs_vb_2) * (pop_size_2 / total_pop)

    loss_total_clinical_burden = obs_cb_1 + obs_cb_2 + obs_vb_1 + obs_vb_2
    loss_equity_of_clinical_burden = abs(exp_cb_1 - obs_cb_1) + abs(exp_cb_2 - obs_cb_2)
    loss_equity_of_vaccination_burden = abs(exp_vb_1 - obs_vb_1) + abs(
        exp_vb_2 - obs_vb_2
    )

    return (
        (1 - a - b) * loss_total_clinical_burden
        + a * loss_equity_of_clinical_burden
        + b * loss_equity_of_vaccination_burden
    )
