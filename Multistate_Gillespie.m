function [t_array, x0_array]=Multistate_Gillespie(k_array,q_array,Dobcktrck,t0,B,Trace_time)

max_index = floor(Trace_time/t0) + 1;

x0_array=zeros([1,max_index]);

n = 1;

index = 1;

t = 0;

q_array_tot = [(1-sum(q_array)), q_array];

P_array = zeros([1, length(q_array_tot)+1]);

for i=1:length(q_array_tot)
    
    P_array(i+1)=P_array(i)+q_array_tot(i);
end

t_array = t0 * (0:max_index-1);

while 1
    
    random = rand;
    
    k = k_array(1);
    
    k_index = 1;
    
    for i=2:length(q_array_tot)
        
        if (P_array(i) < random && random <= P_array(i+1))
            
            k = k_array(i);
            
            k_index = i;
        end
    end
    
    if (Dobcktrck > 0 && k_index == Dobcktrck+1)
        
        k_f = 2*k/(1+B);

        k_b = B*k_f;
        
        k_bt_step = k_f + k_b;
        
        p_b = k_f/k_bt_step;
        
        next_n = n + 1;
        
        while n < next_n
            
            t = t -log(rand)/k_bt_step;
            
            new_index = floor(t/t0) + 1;
     
            if new_index > max_index
                
                x0_array(index:max_index) = n;
                
                break
                
            end
            
            x0_array(index:new_index) = n;
     
            index = new_index + 1;
            
            if rand < p_b
                n = n + 1;
            else
                n = n - 1;
            end
            
        end
        
        if new_index > max_index
            
            break
            
        end
        
    else
        
     t = t -log(rand)/k;
     
     new_index = floor(t/t0) + 1;
     
     if new_index > max_index
         
         x0_array(index:max_index) = n;
         
         break
         
     end
     
     x0_array(index:new_index) = n;
     
     index = new_index + 1;
     
     n = n + 1;
        
    end
    
end
