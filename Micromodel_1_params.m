
function Params = Micromodel_1_params(k_array, q_array)

q_array(2) = q_array(2)*exp(-k_array(3)/k_array(1));

q_array(3) = q_array(3)*exp(-k_array(4)/k_array(1));

kel = k_array(1);

P = sum(q_array);

P2 = (q_array(2)+q_array(3))*(1-P)/P;

k1_0 = k_array(2)/(1-P);

kp2 = P2*k1_0;

ke1 = (1-P2)*k1_0;

P3 = q_array(3)*(1-P)/(P*P2);

kp3 = P3*k_array(3);

ke2 = (1-P3)*k_array(3);

ke3 = k_array(4);

Params = [kel, P, ke1, ke2, ke3 kp2 kp3];