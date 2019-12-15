function mlfunc=ML_func_CG_v2(params, Weights, bincenter, gamma_nums, Nw, t_strt, t_end, Dobcktrck)

pause_num = length(gamma_nums);

[q_array, k_array, sigma_ln, t_bar, B] = MLparams2modelparams_v2(params, pause_num);

Model_at_bincenters = Calculate_Model_v2(bincenter, q_array, k_array, gamma_nums, sigma_ln, t_bar, Nw, B, t_strt, t_end, Dobcktrck);

Model_at_bincenters(Model_at_bincenters == 0 | isnan(Model_at_bincenters) | Model_at_bincenters == Inf) = 10^(-300);

mlfunc = - dot(Weights, log(Model_at_bincenters));