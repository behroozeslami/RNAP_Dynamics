function MLgrid = make_1Dgrid(x,grid_name)

MLgrid = transpose(x);

save([grid_name '.txt'], 'MLgrid', '-ascii') 