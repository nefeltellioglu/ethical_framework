# try out implementation using only socialmixr: 



library(socialmixr)

polymod_contact_matrix = contact_matrix(polymod, 
                                        countries = "United Kingdom", 
                                        age.limits = c(0, 70), 
                                        symmetric = TRUE,
                                        per.capita = TRUE)




per_capita_matrix = polymod_contact_matrix[2]

matrix_sum = sum(per_capita_matrix[[1]])

per_capita_matrix[[1]] = per_capita_matrix[[1]] * (1/matrix_sum)



