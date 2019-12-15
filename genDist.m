function x = genDist(P, bins)

CDF = [0; cumsum(P)];

r = rand;

for j=1:length(P)
    
    if (r > CDF(j)) && (r < CDF(j+1))
        
        x = bins(j);
        
    end
    
end
