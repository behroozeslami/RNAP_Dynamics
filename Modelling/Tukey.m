function [Low_Bond, High_Bond] = Tukey(x, k)

    Q1 = quantile(x, 0.25);
    
    Q3 = quantile(x, 0.75);
    
    delta_Q = Q3-Q1;
    
    Low_Bond = Q1-k*delta_Q;
    
    High_Bond = Q3+k*delta_Q;