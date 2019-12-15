function [DWT_array]=find_Dwell_time(Tc_array)

DWT_array=zeros([1,length(Tc_array)-1]);

for i=1:length(DWT_array)
    
    DWT_array(i)=Tc_array(i+1)-Tc_array(i);
    
end

