function [M,L,U,STD] = Conf_int_4_combined(x,alpha)

M = mean(x);

L = quantile(x,alpha);

U =  quantile(x,1-alpha);

STD = std(x);