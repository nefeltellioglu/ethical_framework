# Using Conmat to generate some reasonable contact matrices

#install.packages("conmat", repos = "https://idem-lab.r-universe.dev")

library(socialmixr)
library(conmat)
# NUMBER of contacts, and number of people in FROM group
polymod_contact_data <- get_polymod_contact_data(setting = "all") 
# population total in age group
polymod_survey_data <- get_polymod_population()

# Model of the contact RATE
set.seed(2024 - 09 - 06)
contact_model <- fit_single_contact_model(
  contact_data = polymod_contact_data,
  population = polymod_survey_data
)

# Source population to predict for (population by age)
nsw <- abs_age_state("NSW")

set.seed(2024 - 09 - 06)
synthetic_contact_nsw <- predict_contacts(
  model = contact_model,
  population = nsw,
  age_breaks = c(0, 70, Inf)
)

# Synthetic contact matrix - 
# "The contact column is the predicted number of contacts from the specified age group to the other one."
synthetic_contact_nsw

# Convert to our pop props (lazy and not quite right)

sevpl<- sum(nsw$population[nsw$lower.age.limit>=70])
sevmi <- sum(nsw$population[nsw$lower.age.limit<70])

synthetic_contact_nsw[1,3] <- synthetic_contact_nsw[1,3]/(sevmi/(sevmi+sevpl))
synthetic_contact_nsw[2,3] <- synthetic_contact_nsw[2,3]/(sevpl/(sevmi+sevpl))
synthetic_contact_nsw[3,3] <- synthetic_contact_nsw[3,3]/(sevmi/(sevmi+sevpl))
synthetic_contact_nsw[4,3] <- synthetic_contact_nsw[4,3]/(sevpl/(sevmi+sevpl))
synthetic_contact_nsw

synthetic_contact_nsw[1,3] <- synthetic_contact_nsw[1,3]*0.9
synthetic_contact_nsw[2,3] <- synthetic_contact_nsw[2,3]*0.1
synthetic_contact_nsw[3,3] <- synthetic_contact_nsw[3,3]*0.9
synthetic_contact_nsw[4,3] <- synthetic_contact_nsw[4,3]*0.1
synthetic_contact_nsw

