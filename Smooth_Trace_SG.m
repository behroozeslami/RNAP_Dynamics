function [xs_array]=Smooth_Trace_SG(x_array,W, k)

xs_array = sgolayfilt(x_array,k,2*W+1);
