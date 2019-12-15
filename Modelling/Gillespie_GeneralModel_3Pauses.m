function [t_array, x0_array, s_array] = Gillespie_GeneralModel_3Pauses(Params, t0, Trace_time)

Dobcktrck = 0;

B = 1;

k_array = Params(1:4);

q_array = Params(5:7);

[t_array, x0_array] = Multistate_Gillespie(k_array,q_array,Dobcktrck,t0,B,Trace_time);

s_array = zeros([1, length(t_array)]);