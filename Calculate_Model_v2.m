 function y = Calculate_Model_v2(t_array, q_array, k_array, gamma_nums, sigma_ln, t_bar, Nw, B, t_strt, t_end, Dobcktrck)

pause_k_array = k_array(2:(length(q_array)+1));

q_array_tot = [1-sum(q_array),q_array];

p_array = cumsum([1-sum(q_array), q_array]);

p_array = p_array(1:length(q_array));

k_array_for_tau = k_array;

if Dobcktrck > 0
    
    k = pause_k_array(Dobcktrck);

    k_f = 2*k/(1+B);

    k_b = B*k_f;
    
    k_array_for_tau(Dobcktrck+1) = abs(k_f-k_b);
    
end

qoverk = q_array_tot ./k_array_for_tau;

n=1:length(q_array);

tau_array = arrayfun(@(i) sum(qoverk(1:i))/p_array(i), n);

mu = log(t_bar) + sigma_ln^2;

tail_func = @(t) (t >= t_bar)*Pause_shoulder_sum(t, pause_k_array, q_array, p_array, tau_array, Nw, gamma_nums, Dobcktrck, B) + ...
            ((t_bar <= t) & (t < (Nw-1)*tau_array(1)))*nchoosek(Nw,1)*q_array(1)*p_array(1)^(Nw-1)*exp(pause_k_array(1)*(Nw-1)*tau_array(1))*Gamma(t,pause_k_array(1),1);

Int_tail = Integrate_Pause_shoulder_sum(t_bar, t_end, pause_k_array, q_array, p_array, tau_array, Nw, gamma_nums, Dobcktrck, B) + ...
                (t_bar < (Nw-1)*tau_array(1))*nchoosek(Nw,1)*q_array(1)*p_array(1)^(Nw-1)*exp(pause_k_array(1)*(Nw-1)*tau_array(1))*IntegrateGamma(t_bar,(Nw-1)*tau_array(1),pause_k_array(1),1);
            
Qc = arrayfun(tail_func, t_bar)/lognpdf(t_bar,mu,sigma_ln);

cutoff_func = @(t) (t < t_bar)*Qc*lognpdf(t,mu,sigma_ln);

Int_cutoff = Qc*0.5*(erf((log(t_bar)-mu)/(sqrt(2)*sigma_ln)) - erf((log(t_strt)-mu)/(sqrt(2)*sigma_ln)));

Q = 2*(1 - Int_tail - Int_cutoff)/(erf((log(t_end)-mu)/(sqrt(2)*sigma_ln)) - erf((log(t_strt)-mu)/(sqrt(2)*sigma_ln)));

peak_func = @(t) Q*lognpdf(t,mu,sigma_ln);

Model = @(t) (tail_func(t) + cutoff_func(t) + peak_func(t));

if Q < 0
    y = 10^(-300)*ones([1, length(t_array)]);
else
    y = arrayfun(Model, t_array);
    
    %y(isnan(y)) = 10^(-300);
end




    

