function MLgrid = make_3Dgrid(x,y,z,grid_name)

MLgrid =zeros([length(x)*length(y)*length(z), 3]);

n = 0;

for i=1:length(x)
    
    for j=1:length(y)
        
        for k=1:length(z)
            
            n = n+1;
        
            MLgrid(n,1) = x(i);
        
            MLgrid(n,2) = y(j);
            
            MLgrid(n,3) = z(k);
            
        end
        
    end
    
end

save([grid_name '.txt'], 'MLgrid', '-ascii') 