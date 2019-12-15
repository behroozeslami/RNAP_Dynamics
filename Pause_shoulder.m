function y = Pause_shoulder(t, k, q, p, tau, Nw, gamma_num, bcktrck_flag, B)

if bcktrck_flag
    
    n_array = 1;
    
    y = sum(arrayfun(@(n) nchoosek(Nw,n)*q^n*p^(Nw-n)*(t>=(Nw-n)*tau)*backtrack(t-(Nw-n)*tau,k,B), n_array));
    
else
    
    n_array = 1:gamma_num;

    y = sum(arrayfun(@(n) nchoosek(Nw,n)*q^n*p^(Nw-n)*(t>=(Nw-n)*tau)*Gamma(t-(Nw-n)*tau,k,n), n_array));
    
end
