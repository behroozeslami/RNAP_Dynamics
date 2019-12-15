function [new_bins, new_bar_pos] = remove_empty_bins(bins, DWT_hist)

new_bins = bins(1:length(bins)-1);

new_bins = new_bins(DWT_hist > 0);

Nbins = length(new_bins);

new_bar_pos=zeros([1,Nbins-1]);

for i=1:Nbins-1
    
    new_bar_pos(i)=sqrt(new_bins(i)*new_bins(i+1));
    
end

