function P = P_F_formula(x,delta,Pmin_ref, K_ref,c)

f0 = 12;

f = x;

f_ref = 7.5;

df =(f-f_ref)/f0;

Pmin = (1+(Pmin_ref.^(-1)-1).*exp(df*delta)).^-1;

K = (Pmin./Pmin_ref).*K_ref.*exp(-df*(1-delta));

P = 1-(1-Pmin)./(1+K/c);
