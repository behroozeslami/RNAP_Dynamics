P_array = zeros([1,5]);

C_array = [1, 10, 100, 500, 1000];

LE_array = zeros([1, 5]);

UE_array = zeros([1, 5]);

EM_array = zeros([1, 5]);

cam_freq = 25;

t0 = 1/cam_freq;

W = 12;

Nw = 4;

slash = '/';

% 1uM

t_strt = 0.05;

t_end = 10000;

output_foldername = 'NEW_AF_7.5pN_1uM_NTP';

save_folder = ['Analysis_1uM_SingleTraces_kel-P' slash '1uM_LV_1'];

data_name = '1uM_LV_1';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(1) = Fit_dat(1,2);

LE_array(1) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(1) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(1) = Fit_dat(4,2);

% 10uM

t_strt = 0.05;

t_end = 10000;

output_foldername = 'NEW_AF_7.5pN_10uM_NTP';

save_folder = ['Analysis_10uM_SingleTraces_kel-P' slash '10uM_LV_1'];

data_name = '10uM_LV_1';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(2) = Fit_dat(1,2);

LE_array(2) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(2) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(2) = Fit_dat(4,2);

% 100uM

t_strt = 0.05;

t_end = 10000;

output_foldername = 'NEW_AF_7.5pN_100uM_NTP';

save_folder = ['Analysis_100uM_SingleTraces_kel-P' slash '100uM_MV_2'];

data_name = '100uM_MV_2';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(3) = Fit_dat(1,2);

LE_array(3) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(3) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(3) = Fit_dat(4,2);

% 500uM

t_strt = 0.05;

t_end = 10000;

output_foldername = 'EXPORT_AF_7.5pN_0.5mM_NTP';

save_folder = ['Analysis_500uM_SingleTraces_kel-P' slash '0.5mM_MV_2'];

data_name = '0.5mM_MV_2';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(4) = Fit_dat(1,2);

LE_array(4) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(4) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(4) = Fit_dat(4,2);

% 1mM

t_strt = 0.05;

t_end = 10000;

output_foldername = 'Export_AF_7.5pN_1mM_NTP';

save_folder = ['Analysis_1mM_SingleTraces_kel-P' slash '1mM_MV_2'];

data_name = '1mM_MV_2';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

P_array(5) = Fit_dat(1,2);

LE_array(5) = Fit_dat(1,2)-Fit_dat(2,2);

UE_array(5) = Fit_dat(3,2)-Fit_dat(1,2);

EM_array(5) = Fit_dat(4,2);

x = C_array(1:5)';

y = P_array(1:5)';

MM = fittype('1-pf_min./(1+K./x)');

[MM_fit,gof,fitinfo] = fit(x,y,MM,'Startpoint',[10, 1]);

PM = MM_fit.pf_min;

KM = MM_fit.K;

P_predict = arrayfun(@(C) 1-PM/(1+KM/C), 0.5:0.01:1000);

PL = figure();

plot(0.5:0.01:1000, P_predict, 'k', 'LineWidth', 2)

hold on 

errorbar(C_array, P_array, LE_array, UE_array, 'rO','MarkerSize',6, 'MarkerFaceColor', 'r')

hold off

xlabel('[NTP] (uM)')

ylabel('P') 

xlim([-5,1005])

folder_name = 'NTP_dependency';

if ~(exist(folder_name,'dir')==7)
    
    mkdir(folder_name)
    
end

str_save = [folder_name '/' 'P1-NTP'];

saveas(PL,[str_save,'.fig'], 'fig')

saveas(PL,[str_save,'.jpg'], 'jpg')






