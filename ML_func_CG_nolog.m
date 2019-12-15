function mlfunc=ML_func_CG_nolog(params, B, Weights, bincenter, gamma_nums, Nw, t_strt, t_end, Dobcktrck)

pause_num = length(gamma_nums);

k_array = params(1:pause_num+1);

q_array = params(pause_num+2:2*pause_num+1);

sigma_ln = params(2*pause_num+2);

t_bar = params(2*pause_num+3);

Model_at_bincenters = Calculate_Model_v2(bincenter, q_array, k_array, gamma_nums, sigma_ln, t_bar, Nw, B, t_strt, t_end, Dobcktrck);

Model_at_bincenters(Model_at_bincenters == 0 | isnan(Model_at_bincenters) | Model_at_bincenters == Inf) = 10^(-300);

mlfunc = - dot(Weights, log(Model_at_bincenters));