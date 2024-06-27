# Using Conmat to generate some reasonable contact matrices

install.packages("conmat", repos = "https://idem-lab.r-universe.dev")

library(socialmixr)
library(conmat)
polymod_contact_data <- get_polymod_contact_data(setting = "all")
polymod_survey_data <- get_polymod_population()

set.seed(2024 - 09 - 06)
contact_model <- fit_single_contact_model(
  contact_data = polymod_contact_data,
  population = polymod_survey_data
)

nsw <- abs_age_state("NSW")

set.seed(2024 - 09 - 06)
synthetic_contact_nsw <- predict_contacts(
  model = contact_model,
  population = nsw,
  age_breaks = c(0, 70, Inf)
)

synthetic_contact_nsw

# Make proportions
synthetic_contact_nsw$contacts <- synthetic_contact_nsw$contacts/
  sum(synthetic_contact_nsw$contacts)
synthetic_contact_nsw

