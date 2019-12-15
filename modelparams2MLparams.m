function params = modelparams2MLparams(q_array, k_array, t_bar, sigma_ln)

p = 1-sum(q_array);

q_ratio = q_array/p;

pq = [0,-log10(q_ratio)];

pk = [0,-log10(k_array)];

k_param = pk(2:length(pk)) - pk(1:length(pk)-1);

q_param =  pq(2:length(pq)) - pq(1:length(pq)-1);

params = [k_param, q_param, t_bar, sigma_ln];