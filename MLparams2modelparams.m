function [q_array, k_array, sigma_ln, t_bar] = MLparams2modelparams(params, pause_num)

k_params = params(1:(pause_num+1));

q_ratio_params = params((pause_num+2):(2*pause_num +1));

t_bar = params(2*pause_num+2);

sigma_ln = params(2*pause_num+3);

pk = cumsum(k_params);

pq_ratio = cumsum(q_ratio_params);

k_array = 10.^(-pk);

q_ratio = 10.^(-pq_ratio);

sum_q_ratio = sum([1, q_ratio]);

q_array = q_ratio/sum_q_ratio;