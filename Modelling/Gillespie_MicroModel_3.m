
function [t_array, x0_array, s_array] = Gillespie_MicroModel_3(Params, t0, Trace_time)

kel = Params(1);

P = Params(2);

ke1 = Params(3);

ke2 = Params(4);

ke3 = Params(5);

q = Params(6);

kp3 = Params(7);

max_index = floor(Trace_time/t0)+1;

T_tot = (max_index-1)*t0;

t_array = 0:t0:T_tot;

x0_array=zeros([1,max_index]);

s_array=zeros([1,max_index]);

s = 0;

n = 0;

T = 0;

index = 1;

while T < T_tot
    
        if s == 0
            
           if rand < (1-P)
         
                new_n = n + 1;
            
                new_s = s;
                
                T = T - log(rand)/kel;
            
           else
                
            if rand < (1-q)
            
                    new_n = n;
            
                    new_s = -1;
            else
                
                new_n = n;
                
                new_s = -2;
                
            end
            
           end

        end
            
        if s == -1
             
            k = ke1;
        
            T = T - log(rand)/k;
         
            new_n = n;
            
            new_s = s + 1;
                
        end
        
        if s == -2
             
            k = ke2 + kp3;
        
            T = T - log(rand)/k;
        
            p = ke2/k;
        
            if rand < p
         
                new_n = n;
            
                new_s = 0;
            
            else
            
                new_n = n;
            
                new_s = s - 1;
            
            end
                
        end
        
        if s == -3
             
            k = ke3;
        
            T = T - log(rand)/k;
                
            new_n = n;
            
            new_s = 0;
         
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