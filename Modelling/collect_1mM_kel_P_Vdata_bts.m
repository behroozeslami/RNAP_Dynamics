grid_name = '1mM_kel_P_grid.txt';

New_grid2D_name = '1mM_kel_P_grid_2.txt';

data_name = '1mM';

save_grid_name = 'V_ML_grid_1mM';

simul_name = 'Simul_Vpdf_1mM';

fit_indices = [1,2];

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Ts = (2*W+1);

T = Ts;

Nw = 10;

Vmin_fit = -20;

Vmax_fit = 30;

output_foldername =  'Select_AF_7.5pN_1mM_NTP';

grid = load(grid_name, '-ascii');

ndat = 100;

nbins = 99;

L = length(grid);

Simul_grid = [grid, zeros([L,nbins])];

grid = [grid, zeros([L,ndat])];

param_num = length(fit_indices);

[Param_dat, txt] = xlsread('Parameters/Micro_Model_1_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx');

Param_names = {};

for i=1:param_num
    
    Param_name = txt(1,fit_indices(i)+1);
    
    Param_names(i) = Param_name;
    
end

Fit_folder = [output_foldername, '/V_Fit_Results_'];

for i=1:param_num
    
    Param_name = Param_names(i);
    
    Fit_folder = [Fit_folder, Param_name{1}];
    
    if i < param_num
        
        Fit_folder = [Fit_folder, '-'];
        
    end
    
end

Fit_folder = [Fit_folder,  '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'T=' num2str(T*t0) 's'];

j = 0;

g = [];

for n=1:L
    
    str_id = [];
    
    for i=1:param_num
    
        Param_name = Param_names(i);
        
        str_id = [str_id, '_', Param_name{1}, '=', num2str(grid(n,i))];
        
    end
    
    str_load = [Fit_folder '/V_Fit_Results' str_id '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'T=' num2str(T*t0) 's' '_' 'limits=' num2str(Vmin_fit) '-' num2str(Vmax_fit) '.txt'];
    
    str_load_2 =[Fit_folder '/Simul_Vpdf' str_id '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'T=' num2str(T*t0) 's' '.txt'];
    
    if ~(exist(str_load, 'file') == 2)
        
        j = j+1;
        
        g(j,:) = grid(n,1:param_num);
        
    else
        
        grid(n,:) = load(str_load,'-ascii');
        
        Simul_grid(n,:) = load(str_load_2,'-ascii');
        
    end
    
end

grid(:,param_num+1:ndat+param_num) = cell2mat(arrayfun(@(col) grid(:,col)-min(grid(:,col)), param_num+1:ndat+param_num, 'UniformOutput', 0));

str_id2 = [];

for i=1:param_num
    
    Param_name = Param_names(i);
    
    str_id2 = [str_id2, Param_name{1}];
    
    if i < param_num
        
        str_id2 = [str_id2, '-'];
        
    end
    
end

folder_out = [output_foldername '/Analysis_Velocity_' data_name '_bts_' str_id2];

if ~(exist(folder_out,'dir')==7)
    
    mkdir(folder_out)
    
end

str_out = [folder_out '/' save_grid_name '_' str_id2 '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'T=' num2str(T*t0) 's' '_' 'limits=' num2str(Vmin_fit) '-' num2str(Vmax_fit) '.mat'];

save(str_out, 'grid', '-mat');

str_out_2 = [folder_out '/' simul_name '_' str_id2 '_' 'Ts='  num2str((2*W+1)*t0) 's' '_' 'T=' num2str(T*t0) '.mat'];

save(str_out_2, 'Simul_grid', '-mat');

g = unique(g,'rows');

save(New_grid2D_name,'g','-ascii')

disp(g)

exit;

