
function [t_array, x0_array, s_array] = Gillespie_MicroModel_2(Params, t0, Trace_time)

f0 = 12;

kel = Params(1);

P = Params(2);

ke1 = Params(3);

kp20 = Params(4);

kf0 = Params(5);

kb0 = Params(6);

delta = Params(7);

f = Params(8);

Wall = Params(9);

kp2 = kp20*exp((7.5-f)*(1-delta)/f0);

kb = kb0*exp((7.5-f)*(1-delta)/f0);
    
kf = kf0*exp((-7.5+f)*delta/f0);

max_index = floor(Trace_time/t0)+1;

T_tot = (max_index-1)*t0;

t_array = 0:t0:T_tot;

x0_array=zeros([1,max_index]);

s_array=zeros([1,max_index]);

s = 0;

n = 0;

T = 0;

BT_length = 0;

index = 1;

while T < T_tot
    
        if s == 0
            
           if rand < (1-P)
         
                new_n = n + 1;
            
                new_s = s;
                
                T = T - log(rand)/kel;
            
            else
            
                new_n = n;
            
                new_s = s - 1;
            
           end

        end
            
        if s == -1
             
            k = ke1 + kp2;
        
            T = T - log(rand)/k;
        
            p = ke1/k;
        
            if rand < p
         
                new_n = n;
            
                new_s = s + 1;
            
            else
            
                new_n = n-1;
            
                new_s = s - 1;
                
                BT_length = BT_length - 1;
            
            end
                
        end
        
        if s <-1
             
            k = kf + kb;
            
            if BT_length == -Wall
                
                k = kf;
                
            end
        
            T = T - log(rand)/k;
        
            p = kf/k;
        
            if rand < p
         
                new_n = n+1;
            
                new_s = s+1;
                
                BT_length = BT_length + 1;
            
            else
            
                new_n = n-1;
            
                new_s = s - 1;
                
                BT_length = BT_length - 1;
            
            end
                
        end
        
        new_index = floor(T/t0) + 1;
     
        if new_index > max_index
        
            new_index = max_index;
        
        end
    
        x0_array(index:new_index) = n;
     
        s_array(index:new_index) = s;
     
        index = new_index + 1;
     
        n = new_n;
     
        s = new_s;
    
end