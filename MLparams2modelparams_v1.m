function [q_array, k_array, sigma_ln, t_bar] = MLparams2modelparams(params, pause_num)

pk = params(1:(pause_num+1));

pq = params((pause_num+2):(2*pause_num +1));

t_bar = params(2*pause_num+2);

sigma_ln = params(2*pause_num+3);

k_array = 10.^(-pk);

q_array = 10.^(-pq);
