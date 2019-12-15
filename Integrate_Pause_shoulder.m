function y = Integrate_Pause_shoulder(t1, t2, k, q, p, tau, Nw, gamma_num, bcktrck_flag, B)

if bcktrck_flag
    
    norm_bt = (B <= 1) + (B > 1)*(1/B);
    
    n_array = 1;
    
    y = sum(arrayfun(@(n) nchoosek(Nw,n)*q^n*p^(Nw-n)*norm_bt, n_array));
    
else
    
    n_array = 1:gamma_num;

    y = sum(arrayfun(@(n) nchoosek(Nw,n)*q^n*p^(Nw-n)*IntegrateGamma(t1-(Nw-n)*tau,t2-(Nw-n)*tau,k,n), n_array));
    
end

