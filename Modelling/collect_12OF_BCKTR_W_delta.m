grid_name = '12OF_BT_wall_delta_grid.txt';

New_grid2D_name = '12OF_BT_wall_delta_grid_2.txt';

data_name = '12OF_BT';

save_grid_name = 'ML_grid_12OF_BT';

simul_name = 'Simul_DWTpdf_12OF_BT';

proc_name = 'Processivity_12OF_BT';

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Nw = 4;

t_strt = 0.05;

t_end = 10000;

output_foldername = 'Select_OF_12.5pN_1mM_NTP';

str_DWT_bins = [output_foldername '/' 'DWT_pdf_bincenters' '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];

grid = load(grid_name, '-ascii');

ndat = 1;

B = load(str_DWT_bins, '-mat');

bins = B.centers;

nbins = length(bins);

L = length(grid);

Simul_grid = [grid, zeros([L,nbins])];

Proc_grid = [grid, zeros([L,2])];

grid = [grid, zeros([L,ndat])];

param_num = 2;

Param_names = {'W', 'delta'};

Fit_folder = [output_foldername, '/Fit_Results_'];

for i=1:param_num
    
    Param_name = Param_names(i);
    
    Fit_folder = [Fit_folder, Param_name{1}];
    
    if i < param_num
        
        Fit_folder = [Fit_folder, '-'];
        
    end
    
end

Fit_folder = [Fit_folder,  '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw)];

j = 0;

g = [];

for n=1:L
    
    str_id = [];
    
    for i=1:param_num
    
        Param_name = Param_names(i);
        
        str_id = [str_id, '_', Param_name{1}, '=', num2str(grid(n,i))];
        
    end
    
    str_load = [Fit_folder '/Fit_Results' str_id '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.txt'];
    
    str_load_2 =[Fit_folder '/Simul_DWTpdf' str_id '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.txt'];
    
    str_load_3 =[Fit_folder '/Processivity' str_id '.txt'];
    
    if ~(exist(str_load, 'file') == 2)
        
        j = j+1;
        
        g(j,:) = grid(n,1:param_num);
        
    else
        
        grid(n,:) = load(str_load,'-ascii');
        
        Simul_grid(n,:) = load(str_load_2,'-ascii');
        
        Proc_grid(n,:) = load(str_load_3,'-ascii');
        
    end
    
end

%grid(:,param_num+1:ndat+param_num) = cell2mat(arrayfun(@(col) grid(:,col)-min(grid(:,col)), param_num+1:ndat+param_num, 'UniformOutput', 0));

str_id2 = [];

for i=1:param_num
    
    Param_name = Param_names(i);
    
    str_id2 = [str_id2, Param_name{1}];
    
    if i < param_num
        
        str_id2 = [str_id2, '-'];
        
    end
    
end

folder_out = [output_foldername '/Analysis_' data_name '_' str_id2];

if ~(exist(folder_out,'dir')==7)
    
    mkdir(folder_out)
    
end

str_out = [folder_out '/' save_grid_name '_' str_id2 '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.mat'];

save(str_out, 'grid', '-mat');

str_out_2 = [folder_out '/' simul_name '_' str_id2 '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '.mat'];

save(str_out_2, 'Simul_grid', '-mat');

str_out_3 = [folder_out '/' proc_name '_' str_id2 '.mat'];

save(str_out_3, 'Proc_grid', '-mat');

g = unique(g,'rows');

save(New_grid2D_name,'g','-ascii')

disp(g)

exit;

