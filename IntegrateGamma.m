function I = IntegrateGamma(t1,t2,k,Nw)

if t1 < 0
    t1 = 0;
end

n_array = 1:Nw;

I = sum(arrayfun(@(n) (Gamma(t1,k,n)-Gamma(t2,k,n))/k, n_array));

if t2 < 0
    I = 0;
end