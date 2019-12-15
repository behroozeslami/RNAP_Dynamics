[Param_dat, txt] = xlsread('Parameters/Micro_Model_1_Parameters_Select_AF_7.5pN_1mM_NTP.xlsx');

Params = Param_dat(1,:);

errors = Param_dat(2,:);

force_indices = 6;

kb = Params(force_indices(1));

err = errors(force_indices(1));

delta = 0.7;

f0 = 12;

F_array = [-12.5, -10, -9, -7.5, -5, 5, 7.5, 10 12.5];

kp2_array = kb*exp((7.5-F_array)*(1-delta)/f0);

err_array = err*exp((7.5-F_array)*(1-delta)/f0);

kp2_dat = [F_array', kp2_array', err_array']';

folder_name = 'Force_dependency';

if ~(exist(folder_name,'dir')==7)
    
    mkdir(folder_name)
    
end

str_save = [folder_name '/' 'kp2-F'];

row_header(1:3,1) = {'Force (pN)', 'kp2', 'error'};

xlswrite([str_save '.xlsx'],kp2_dat,'Sheet1','B1');

xlswrite([str_save '.xlsx'],row_header,'Sheet1','A1'); 
