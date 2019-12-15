function MLgrid = make_2Dgrid(x,y,grid_name)

MLgrid =zeros([length(x)*length(y), 2]);

n = 0;

for i=1:length(x)
    
    for j=1:length(y)
        
        n = n+1;
        
        MLgrid(n,1) = x(i);
        
        MLgrid(n,2) = y(j);
        
    end
    
end

save([grid_name '.txt'], 'MLgrid', '-ascii') 