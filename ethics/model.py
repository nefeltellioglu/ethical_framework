import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
from dataclasses import dataclass
import numpy as np
import pandas as pd
import typing as tp


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
    prop_hosp_inf_1: float
    prop_hosp_inf_2: float
    days_hosp_inf_1: float
    days_hosp_inf_2: float
    # Protection from disease among (unprotected) vaccinated.
    # TODO Something about this seems a bit fishy to me, maybe we need
    # a clearer variable name. Is it a proportion???
    vacc_protection_from_disease_1: float
    vacc_protection_from_disease_2: float
    # Vaccination burden parameter
    prop_hosp_vacc_1: float
    prop_hosp_vacc_2: float
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
        if pop == 1:
            return (
                self.s0_1
                + self.i0_1
                + self.r0_1
                + self.s0_1_vp
                + self.s0_1_vu
                + self.i0_1_vu
                + self.r0_1_vu
            )
        elif pop == 2:
            return (
                self.s0_2
                + self.i0_2
                + self.r0_2
                + self.s0_2_vp
                + self.s0_2_vu
                + self.i0_2_vu
                + self.r0_2_vu
            )
        else:
            raise ValueError(f"Invalid population: {pop}.")

    def total_number_vaccinated(self) -> int:
        """
        Returns the total number of vaccinated individuals in both populations.
        """
        num_vac_1, num_vac_2 = self.number_vaccinated_by_group()
        return num_vac_1 + num_vac_2

    def number_vaccinated_by_group(self) -> tp.Tuple[int, int]:
        """
        Returns the number of vaccinated individuals in each population.
        """
        return (self.s0_1_vp + self.s0_1_vu +
                self.i0_1_vu + self.r0_1_vu,
                self.s0_2_vp + self.s0_2_vu +
                self.i0_2_vu + self.r0_2_vu)

    @staticmethod
    def integer_initial_conditions(
        vacc_prop_1: float,
        vacc_prop_2: float,
        pop_size_1: int,
        pop_size_2: int,
        vacc_protection_1: float,
        vacc_protection_2: float,
        initial_inf_rec: dict,
    ):
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
            s0_1=pop_size_1
            - num_vac_1
            - initial_inf_rec["i0_1"]
            - initial_inf_rec["r0_1"],
            s0_2=pop_size_2
            - num_vac_2
            - initial_inf_rec["i0_2"]
            - initial_inf_rec["r0_2"],
            i0_1=initial_inf_rec["i0_1"],
            i0_2=initial_inf_rec["i0_2"],
            r0_1=initial_inf_rec["r0_1"],
            r0_2=initial_inf_rec["r0_2"],
            s0_1_vp=num_vac_1_protected,
            s0_1_vu=num_vac_1
            - num_vac_1_protected
            - initial_inf_rec["i0_1_vu"]
            - initial_inf_rec["r0_1_vu"],
            s0_2_vp=num_vac_2_protected,
            s0_2_vu=num_vac_2
            - num_vac_2_protected
            - initial_inf_rec["i0_2_vu"]
            - initial_inf_rec["r0_2_vu"],
            i0_1_vu=initial_inf_rec["i0_1_vu"],
            i0_2_vu=initial_inf_rec["i0_2_vu"],
            r0_1_vu=initial_inf_rec["r0_1_vu"],
            r0_2_vu=initial_inf_rec["r0_2_vu"],
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
            "total": self.s1_vp[0] + self.s1_vu[0] + self.s2_vp[0] + self.s2_vu[0],
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

        pop_1 = (
            self.s1[0]
            + self.i1[0]
            + self.r1[0]
            + self.s1_vp[0]
            + self.s1_vu[0]
            + self.i1_vu[0]
            + self.r1_vu[0]
        )

        pop_2 = (
            self.s2[0]
            + self.i2[0]
            + self.r2[0]
            + self.s2_vp[0]
            + self.s2_vu[0]
            + self.i2_vu[0]
            + self.r2_vu[0]
        )
        return {
            "pop_1": pop_1,
            "pop_2": pop_2,
            "total": pop_1 + pop_2,
        }


def sir_vacc(params: SIRParams, sir_0: SIRInitialCondition, ts) -> [SIRSolution]:
    y0 = [
        sir_0.s0_1,
        sir_0.s0_2,
        sir_0.i0_1,
        sir_0.i0_2,
        sir_0.r0_1,
        sir_0.r0_2,
        sir_0.s0_1_vp,
        sir_0.s0_1_vu,
        sir_0.s0_2_vp,
        sir_0.s0_2_vu,
        sir_0.i0_1_vu,
        sir_0.i0_2_vu,
        sir_0.r0_1_vu,
        sir_0.r0_2_vu,
    ]

    pop_size_1 = (
        sir_0.s0_1
        + sir_0.i0_1
        + sir_0.r0_1
        + sir_0.s0_1_vp
        + sir_0.s0_1_vu
        + sir_0.i0_1_vu
        + sir_0.r0_1_vu
    )
    pop_size_2 = (
        sir_0.s0_2
        + sir_0.i0_2
        + sir_0.r0_2
        + sir_0.s0_2_vp
        + sir_0.s0_2_vu
        + sir_0.i0_2_vu
        + sir_0.r0_2_vu
    )

    def deriv(y, t, params):
        (
            s_1,
            s_2,
            i_1,
            i_2,
            r_1,
            r_2,
            s_1_vp,
            s_1_vu,
            s_2_vp,
            s_2_vu,
            i_1_vu,
            i_2_vu,
            r_1_vu,
            r_2_vu,
        ) = y

        # # unvaccinated  I infect unvaccinated S
        inf_11 = params.beta_11 * i_1 * s_1 / pop_size_1
        inf_12 = params.beta_12 * i_1 * s_2 / pop_size_2
        inf_21 = params.beta_21 * i_2 * s_1 / pop_size_1
        inf_22 = params.beta_22 * i_2 * s_2 / pop_size_2

        # vaccinated unprotected I (Ivu) infect unvaccinated S (S)
        inf_1vu1 = params.beta_11 * i_1_vu * s_1 / pop_size_1
        inf_1vu2 = params.beta_12 * i_1_vu * s_2 / pop_size_2
        inf_2vu1 = params.beta_21 * i_2_vu * s_1 / pop_size_1
        inf_2vu2 = params.beta_22 * i_2_vu * s_2 / pop_size_2

        # unvaccinated I infect vaccinated unprotected S (Svu)
        inf_11vu = params.beta_11 * i_1 * s_1_vu / pop_size_1
        inf_12vu = params.beta_12 * i_1 * s_2_vu / pop_size_2
        inf_21vu = params.beta_21 * i_2 * s_1_vu / pop_size_1
        inf_22vu = params.beta_22 * i_2 * s_2_vu / pop_size_2

        # vaccinated unprotected I (Ivu) infect vaccinated unprotected S (Svu)
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

        return [
            ds_1,
            ds_2,
            di_1,
            di_2,
            dr_1,
            dr_2,
            ds_1_vp,
            ds_1_vu,
            ds_2_vp,
            ds_2_vu,
            di_1_vu,
            di_2_vu,
            dr_1_vu,
            dr_2_vu,
        ]

    result = scipy.integrate.odeint(
        deriv, y0, ts, args=(params,)
    )  # This needs to change probably.
    return [
        SIRSolution(
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
        )
    ]


def loss_terms(
    outcome: SIROutcome,
    ic: SIRInitialCondition,
    burden_params: BurdenParams,
) -> (float, float, float):
    """
    The three loss terms: total clinical burden, loss of equity in
    clinical burden, and loss of equity in vaccination burden.
    Clinical burden includes both the infection burden and vaccination
    burden.

    Args:
    outcome: `SIROutcome`
    ic: `SIRInitialCondition`
    burden_params: `BurdenParams`

    Details:
    This function provides the terms needed to compute the following
    loss function:

    (
         (1 - a - b) * loss_total_clinical_burden
         + a * loss_equity_of_clinical_burden
         + b * loss_equity_of_vaccination_burden
    )
    """
    pop_size_1, pop_size_2 = ic.pop_size(1), ic.pop_size(2)
    total_pop = float(pop_size_1 + pop_size_2)

    obs_vb_1 = (
        outcome.total_vac_1
        * burden_params.prop_hosp_vacc_1
        * burden_params.days_hosp_vacc_1
    )
    obs_vb_2 = (
        outcome.total_vac_2
        * burden_params.prop_hosp_vacc_2
        * burden_params.days_hosp_vacc_2
    )
    obs_ib_1_no_vac = (
        outcome.inf_1_no_vac
        * burden_params.prop_hosp_inf_1
        * burden_params.days_hosp_inf_1
    )
    obs_ib_2_no_vac = (
        outcome.inf_2_no_vac
        * burden_params.prop_hosp_inf_2
        * burden_params.days_hosp_inf_2
    )
    obs_ib_1_vu = (
        outcome.inf_1_vu
        * burden_params.prop_hosp_inf_1
        * burden_params.days_hosp_inf_1
        * (1 - burden_params.vacc_protection_from_disease_1)
    )
    obs_ib_2_vu = (
        outcome.inf_2_vu
        * burden_params.prop_hosp_inf_2
        * burden_params.days_hosp_inf_2
        * (1 - burden_params.vacc_protection_from_disease_2)
    )

    obs_vb_tot = obs_vb_1 + obs_vb_2

    obs_cbi_1 = obs_ib_1_no_vac + obs_ib_1_vu
    obs_cbi_2 = obs_ib_2_no_vac + obs_ib_2_vu
    obs_cbi_tot = obs_cbi_1 + obs_cbi_2

    exp_cbi_1 = obs_cbi_tot * (pop_size_1 / total_pop)
    exp_cbi_2 = obs_cbi_tot * (pop_size_2 / total_pop)
    exp_vb_1 = obs_vb_tot * (pop_size_1 / total_pop)
    exp_vb_2 = obs_vb_tot * (pop_size_2 / total_pop)

    loss_total_clinical_burden = obs_cbi_1 + obs_cbi_2 + obs_vb_1 + obs_vb_2
    loss_equity_of_infection_burden = abs(exp_cbi_1 - obs_cbi_1) + abs(exp_cbi_2 - obs_cbi_2)
    loss_equity_of_vaccination_burden = abs(exp_vb_1 - obs_vb_1) + abs(exp_vb_2 - obs_vb_2)

    return (
        loss_total_clinical_burden,
        loss_equity_of_infection_burden,
        loss_equity_of_vaccination_burden,
    )

def global_loss(
    loss_total_clinical_burden: float,
    loss_equity_of_infection_burden: float,
    loss_equity_of_vaccination_burden: float,
    weight_EI:float,
    weight_EV:float
) -> (float):
    loss = (1 - weight_EI - weight_EV) * loss_total_clinical_burden + \
           weight_EI * loss_equity_of_infection_burden + \
           weight_EV * loss_equity_of_vaccination_burden
    return loss





### functions for calculating burdens from outcomes: 


#### ====================================================================
#define function to calculate total burden from SIROutcome object
# burden from adverse vaccination reactions (group 1)

# count of infections
def count_infections_group_1(sir: SIROutcome) -> float:
    return(sir.inf_1_no_vac + sir.inf_1_vu)

def count_infections_group_2(sir: SIROutcome) -> float:
    return(sir.inf_2_no_vac + sir.inf_2_vu)

def count_vaccinations_group_1(sir: SIROutcome) -> float:
    return sir.total_vac_1

def count_vaccinations_group_2(sir: SIROutcome) -> float:
    return sir.total_vac_2



def burden_adverse_group_1(
    sir: SIROutcome, dbp: BurdenParams
    ) -> float:
    return (dbp.days_hosp_vacc_1
            * sir.total_vac_1
            * dbp.prop_hosp_vacc_1)

# burden from adverse vaccination reactions (group 2)
def burden_adverse_group_2(
    sir: SIROutcome, dbp: BurdenParams
    ) -> float:
    return (dbp.days_hosp_vacc_2
            * sir.total_vac_2
            * dbp.prop_hosp_vacc_2)

#burden from infections in unvaccinated people (group 1)
def burden_infections_group_1_noVacc(
    sir: SIROutcome, dbp: BurdenParams
    ) -> float:
    return (dbp.prop_hosp_inf_1 *
            dbp.days_hosp_inf_1 * sir.inf_1_no_vac)

# burden from infections in unvaccinated people (group 2)
def burden_infections_group_2_noVacc(
    sir: SIROutcome, dbp: BurdenParams
    ) -> float:
    return (dbp.prop_hosp_inf_2 *
            dbp.days_hosp_inf_2 * sir.inf_2_no_vac)

#burden from infections in vaccinated people (group 1)
def burden_infections_group_1_Vacc(
    sir: SIROutcome, dbp: BurdenParams
    ) -> float:
    return (dbp.prop_hosp_inf_1 *
            (1 - dbp.vacc_protection_from_disease_1) *
            dbp.days_hosp_inf_1 *
            sir.inf_1_vu )

#burden from infections in vaccinated people (group 2)
def burden_infections_group_2_Vacc(
    sir: SIROutcome, dbp: BurdenParams
    ) -> float:
    return (dbp.prop_hosp_inf_2 *
            (1 - dbp.vacc_protection_from_disease_2) *
            dbp.days_hosp_inf_2 *
            sir.inf_2_vu )

# total infection burden group 1
def total_burden_infections_group_1(
    sir: SIROutcome, dbp: BurdenParams
    ) -> float:
    tot= (burden_infections_group_1_noVacc(sir, dbp) +
            burden_infections_group_1_Vacc(sir, dbp))
    return (tot)

def total_burden_infections_group_2(
    sir: SIROutcome, dbp: BurdenParams
    ) -> float:
    tot= (burden_infections_group_2_noVacc(sir, dbp) +
            burden_infections_group_2_Vacc(sir, dbp))
    return (tot)


def total_vaccinations(sir: SIROutcome) -> float:
    return (sir.total_vac_1 + sir.total_vac_2)

# aggregate burden components
def total_burden_infections(
    sir: SIROutcome, dbp: BurdenParams
    ) -> float:
    tot_1 = total_burden_infections_group_1(sir, dbp)
    tot_2 = total_burden_infections_group_2(sir, dbp)
    return (tot_1 + tot_2)

def total_burden_adverse(sir: SIROutcome, dbp: BurdenParams) -> float:
    return(burden_adverse_group_1(sir, dbp) +
           burden_adverse_group_2(sir, dbp))

def total_clinical_burden(sir:SIROutcome, dbp: BurdenParams) -> float:
    return (total_burden_infections(sir, dbp) +
            total_burden_adverse(sir, dbp))

# per-capita burdens: 
def pop_1(ic: SIRInitialCondition) -> int:
    return ic.pop_size(1)

def pop_2(ic: SIRInitialCondition) -> int:
    return ic.pop_size(2)

def total_pop(ic: SIRInitialCondition) -> int:
    return (pop_1(ic) + pop_2(ic))

def adverse_per_capita_1(sir: SIROutcome, 
                         dbp: BurdenParams, 
                         ic: SIRInitialCondition) -> float:
    return(burden_adverse_group_1(sir, dbp) / pop_1(ic))

def adverse_per_capita_2(sir: SIROutcome, 
                         dbp: BurdenParams, 
                         ic: SIRInitialCondition) -> float:
    return(burden_adverse_group_2(sir, dbp) / pop_2(ic)) 

def infection_burden_per_capita_1(sir: SIROutcome, 
                                  dbp: BurdenParams, 
                                  ic: SIRInitialCondition) -> float:
    return(total_burden_infections_group_1(sir, dbp) / pop_1(ic))

def infection_burden_per_capita_2(sir: SIROutcome, 
                                  dbp: BurdenParams, 
                                  ic: SIRInitialCondition) -> float:
    return(total_burden_infections_group_2(sir, dbp) / pop_2(ic)) 