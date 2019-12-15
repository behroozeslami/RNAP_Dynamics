function y = Integrate_Pause_shoulder_sum(t1, t2, k_array, q_array, p_array, tau_array, Nw, gamma_nums, Dobcktrck, B)

terms = 1:length(gamma_nums);

y = sum(arrayfun(@(i) Integrate_Pause_shoulder(t1, t2, k_array(i), q_array(i), p_array(i), tau_array(i), Nw, gamma_nums(i), Dobcktrck==i, B), terms));