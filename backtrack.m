function y = backtrack(t, k, B)

k_f = 2*k/(1+B);

k_b = B*k_f;

k_bt = sqrt(k_f*k_b);

if t < 0 
    y = 0;
else
    y = k_f*exp(-(k_f + k_b)*t)* besseli(1,2*k_bt*t)/(k_bt*t);
end

if isnan(y) || y == Inf || y == -Inf
    y = 0;
end