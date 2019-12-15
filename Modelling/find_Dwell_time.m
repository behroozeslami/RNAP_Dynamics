function [DWT_array]=find_Dwell_time(Tc_array)

len = length(Tc_array);

DWT_array = Tc_array(2:len) - Tc_array(1:len-1);


