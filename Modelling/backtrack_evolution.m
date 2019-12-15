[Param, txt] = xlsread('Parameters/Micro_Model_1_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx');

ke2 = Param(4);

ke3 = Param(5);

kp2 = Param(6);

kp3 = Param(7);

M = [-kp2,0,0;kp2,-(kp3+ke2),0;0,kp3,-ke3];

Prob0 = [1,0,0]';

t_array = 0.1:0.1:600;

Prob1 = zeros([1,length(t_array)]);

Prob2 = zeros([1,length(t_array)]);

Prob3 = zeros([1,length(t_array)]);

for i=1:length(t_array)
    
    t = t_array(i);
    
    Prob = expm(M*t)*Prob0;
    
    Prob1(i) = Prob(1);
    
    Prob2(i) = Prob(2);
    
    Prob3(i) = Prob(3);
end

BT = Prob2+Prob3;

fraction = Prob3./BT;

figure()

semilogx(t_array/60,fraction,'r', 'LineWidth',1.5)

xlabel('time (min)')

ylabel('P3/(P2+P3)')
