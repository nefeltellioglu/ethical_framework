
import ethics.model as em

#### ====================================================================
#define function to calculate total burden from SIROutcome object
# burden from adverse vaccination reactions (group 1)
def burden_adverse_group_1(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    return (dbp.days_hosp_vacc_1
            * sir.total_vac_1
            * dbp.prop_hosp_vacc_1)

# burden from adverse vaccination reactions (group 2)
def burden_adverse_group_2(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    return (dbp.days_hosp_vacc_2
            * sir.total_vac_2
            * dbp.prop_hosp_vacc_2)

#burden from infections in unvaccinated people (group 1)
def burden_infections_group_1_noVacc(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    return (dbp.prop_hosp_inf_1 *
            dbp.days_hosp_inf_1 * sir.inf_1_no_vac)

# burden from infections in unvaccinated people (group 2)
def burden_infections_group_2_noVacc(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    return (dbp.prop_hosp_inf_2 *
            dbp.days_hosp_inf_2 * sir.inf_2_no_vac)

#burden from infections in vaccinated people (group 1)
def burden_infections_group_1_Vacc(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    return (dbp.prop_hosp_inf_1 *
            (1 - dbp.vacc_protection_from_disease_1) *
            dbp.days_hosp_inf_1 *
            sir.inf_1_vu )

#burden from infections in vaccinated people (group 2)
def burden_infections_group_2_Vacc(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    return (dbp.prop_hosp_inf_2 *
            (1 - dbp.vacc_protection_from_disease_2) *
            dbp.days_hosp_inf_2 *
            sir.inf_2_vu )

# total infection burden group 1
def total_burden_infections_group_1(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    tot= (burden_infections_group_1_noVacc(sir, dbp) +
            burden_infections_group_1_Vacc(sir, dbp))
    return (tot)

def total_burden_infections_group_2(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    tot= (burden_infections_group_2_noVacc(sir, dbp) +
            burden_infections_group_2_Vacc(sir, dbp))
    return (tot)


def total_vaccinations(sir: em.SIROutcome) -> float:
    return (sir.total_vac_1 + sir.total_vac_2)

# aggregate burden components
def total_burden_infections(
    sir: em.SIROutcome, dbp: em.BurdenParams
    ) -> float:
    tot_1 = total_burden_infections_group_1(sir, dbp)
    tot_2 = total_burden_infections_group_2(sir, dbp)
    return (tot_1 + tot_2)

def total_burden_adverse(sir: em.SIROutcome, dbp: em.BurdenParams) -> float:
    return(burden_adverse_group_1(sir, dbp) +
           burden_adverse_group_2(sir, dbp))

def total_clinical_burden(sir:em.SIROutcome, dbp: em.BurdenParams) -> float:
    return (total_burden_infections(sir, dbp) +
            total_burden_adverse(sir, dbp))

