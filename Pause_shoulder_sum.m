function y = Pause_shoulder_sum(t, k_array, q_array, p_array, tau_array, Nw, gamma_nums, Dobcktrck, B)

terms = 1:length(gamma_nums);

y = sum(arrayfun(@(i) Pause_shoulder(t, k_array(i), q_array(i), p_array(i), tau_array(i), Nw, gamma_nums(i), Dobcktrck==i, B), terms));