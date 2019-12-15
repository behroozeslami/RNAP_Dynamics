function [new_col] = remove_outliers_and_average(col, param_num, grid)

new_col = zeros([length(col), 1]);

param_vals = grid(:,1:param_num);

mean_best_vals = zeros([1, param_num]);



for j=1:param_num
    
    param_val = param_vals(:,j);
    
    best_vals = cell2mat(arrayfun(@(I) param_val(I,1)*ones([1,col(I,1)]), 1:length(col), 'UniformOutput', false));
    
    [Low_Bond, High_Bond] = Tukey(best_vals, 1.5);
    
    noOut_best_vals = best_vals(best_vals >= Low_Bond & best_vals <= High_Bond);
    
    mean_best_vals(j) = mean(noOut_best_vals);
    
end

d = arrayfun(@(row) dot((param_vals(row,:)- mean_best_vals), (param_vals(row,:)- mean_best_vals)), 1:length(param_vals));

[M, I] = min(d);

new_col(I,1) = 1; 


