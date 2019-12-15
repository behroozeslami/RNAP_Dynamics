function sum_LL_grid = Calculate_LogLikelihood_Combine_from_SingleTrace(grid,selected_trace_indices,trace_indices,param_num,nboot)

[L1,L2]=size(grid);

combined_grid = zeros([L1,length(selected_trace_indices)]);

sub_grid = grid(:,param_num+1:L2);

sum_LL_grid = zeros([L1,nboot]);

for n=1:nboot

    for i=1:length(selected_trace_indices)

        trace_index = selected_trace_indices(i);

        temp_grid = sub_grid(:,(trace_indices==trace_index));

        combined_grid(:,i) = temp_grid(:,n);
    end
    
    sum_LL_grid(:,n) = sum(combined_grid,2); 
end

