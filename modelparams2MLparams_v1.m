function params = modelparams2MLparams(q_array, k_array, t_bar, sigma_ln)

k_param = -log10(q_array);

q_param = -log10(k_array);

params = [k_param, q_param, t_bar, sigma_ln];