% final size at critical threshold


N1 = 4395000;
N2 = 605000;

R0 = 1 ;

gamma = 0.096;

c = [0.38, 0.14; 0.14, 0.34];

beta = R0 * 2 * gamma / (c(1, 1) + c(2, 2) + ...
                (c(1, 1)^2 - 2*c(2, 2) * c(1, 1) + ...
                c(2, 2)^2 + 4 * c(1, 2)*c(2, 2)) ^0.5);

b = beta * c;

rec_0_1 = 0;
rec_0_2 = 0;
I_0_1 = 1;
I_0_2 = 1;
sus_0_1 = N1 - I_0_1 - rec_0_1;
sus_0_2 = N2 - I_0_2 - rec_0_2 ;

% arbitrary guesses 
R_final_1 = 10;
R_final_2 = 10;

n_iter = 100 ;

R_final_1_vec = NaN(1, n_iter);
R_final_2_vec = NaN(1, n_iter);

R_final_1_vec(1) = R_final_1;
R_final_2_vec(1) = R_final_2;

for i = 2:n_iter 

    R_final_1_vec(i) = N1 - sus_0_1 * exp(-1 * 1/(gamma * N1) *...
                       (b(1, 1) * (R_final_1_vec(i-1) - rec_0_1) + ...
                        b(2, 1) * (R_final_2_vec(i-1) - rec_0_2)));

    R_final_2_vec(i) = N2 - sus_0_2 * exp(-1 * 1/(gamma * N2) *...
                       (b(1, 2) * (R_final_1_vec(i-1) - rec_0_1) + ...
                        b(2, 2) * (R_final_2_vec(i-1) - rec_0_2)));

end 

figure(1)
plot(R_final_1_vec)

figure(2)
plot(R_final_2_vec)

FSz = R_final_2_vec(n_iter) + R_final_1_vec(n_iter)