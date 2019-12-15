function index_2_keep = Exclude_Blacklist(selected_trace_number, Blacklist)

index_2_keep = true([1, length(selected_trace_number)]);

for j=1:length(Blacklist)
    
    bool_vec = selected_trace_number == Blacklist(j);
    
    if sum(bool_vec)
        
        index_2_keep(bool_vec) = false;
        
    end
    
end
