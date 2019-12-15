kel_array = zeros([1,5]);

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

save_folder = ['Analysis_1uM_SingleTraces_kel-P' slash '1uM_SM'];

data_name = '1uM_SM';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

kel_array(1) = Fit_dat(1,1);

LE_array(1) = Fit_dat(1,1)-Fit_dat(2,1);

UE_array(1) = Fit_dat(3,1)-Fit_dat(1,1);

EM_array(1) = Fit_dat(4,1);

% 10uM

t_strt = 0.05;

t_end = 10000;

output_foldername = 'NEW_AF_7.5pN_10uM_NTP';

save_folder = ['Analysis_10uM_SingleTraces_kel-P' slash '10uM_SM'];

data_name = '10uM_SM';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

kel_array(2) = Fit_dat(1,1);

LE_array(2) = Fit_dat(1,1)-Fit_dat(2,1);

UE_array(2) = Fit_dat(3,1)-Fit_dat(1,1);

EM_array(2) = Fit_dat(4,1);

% 100uM

t_strt = 0.05;

t_end = 10000;

output_foldername = 'NEW_AF_7.5pN_100uM_NTP';

save_folder = ['Analysis_100uM_SingleTraces_kel-P' slash '100uM_SM'];

data_name = '100uM_SM';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

kel_array(3) = Fit_dat(1,1);

LE_array(3) = Fit_dat(1,1)-Fit_dat(2,1);

UE_array(3) = Fit_dat(3,1)-Fit_dat(1,1);

EM_array(3) = Fit_dat(4,1);

% 500uM

t_strt = 0.05;

t_end = 10000;

output_foldername = 'EXPORT_AF_7.5pN_0.5mM_NTP';

save_folder = ['Analysis_500uM_SingleTraces_kel-P' slash '500uM_SM'];

data_name = '500uM_SM';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

kel_array(4) = Fit_dat(1,1);

LE_array(4) = Fit_dat(1,1)-Fit_dat(2,1);

UE_array(4) = Fit_dat(3,1)-Fit_dat(1,1);

EM_array(4) = Fit_dat(4,1);

% 1mM

t_strt = 0.05;

t_end = 10000;

output_foldername = 'Export_AF_7.5pN_1mM_NTP';

save_folder = ['Analysis_1mM_SingleTraces_kel-P' slash '1mM_SM'];

data_name = '1mM_SM';

folder_grid = [output_foldername slash save_folder];

str_FitResult = [folder_grid slash 'Fit_outcome' '_' data_name '_' 'kel' '-' 'P' '_' 'Ts=' num2str((2*W+1)*t0) 's' '_' 'Nw=' num2str(Nw) '_' 'limits=' num2str(t_strt) '-' num2str(t_end) '.xlsx'];

Fit_dat = xlsread(str_FitResult);

kel_array(5) = Fit_dat(1,1);

LE_array(5) = Fit_dat(1,1)-Fit_dat(2,1);

UE_array(5) = Fit_dat(3,1)-Fit_dat(1,1);

EM_array(5) = Fit_dat(4,1);

x = C_array(1:5)';

y = kel_array(1:5)';

MM = fittype('Vmax./(1+K./x)');

[MM_fit,gof,fitinfo] = fit(x,y,MM,'Startpoint',[50, 25]);

VM = MM_fit.Vmax;

KM = MM_fit.K;

kel_predict = arrayfun(@(C) VM/(1+KM/C), 0.5:0.01:1000);

PL = figure();

plot(0.5:0.01:1000, kel_predict, 'k', 'LineWidth', 2)

hold on 

errorbar(C_array, kel_array, LE_array, UE_array, 'rO','MarkerSize',6, 'MarkerFaceColor', 'r')

errorbar(C_array, kel_array, EM_array, 'k.')

hold off

xlabel('[NTP] (uM)')

ylabel('k_{el} (bp/s)') 

xlim([-5,1005])

folder_name = 'NTP_dependency_NewSel';

if ~(exist(folder_name,'dir')==7)
    
    mkdir(folder_name)
    
end

str_save = [folder_name '/' 'kel-NTP'];

saveas(PL,[str_save,'.fig'], 'fig')

saveas(PL,[str_save,'.jpg'], 'jpg')






