delta = 1;

Wall = 5;

f = -12.5;

kel = 18.5;

P = 0.29;

Trace_time = 2100;

trace_num = 5;

cam_freq = 25;

t0 = 1/cam_freq;

mainfoldername = 'Select_OF_12.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

if ~(exist(output_foldername,'dir')==7)
    
    mkdir(output_foldername)
    
end

[Param_dat, txt] = xlsread('Parameters/Backtrack_Model_Parameters_raw_Select_AF_7.5pN_1mM_NTP.xlsx');

Params = Param_dat(1,:);

Param_name_1 = {'W'};

Param_name_2 = {'delta'};

trace_folder = [output_foldername slash 'Simulated_traces_12OF'];

if ~(exist(trace_folder,'dir')==7)
    
    mkdir(trace_folder)
    
end

str_trace = [trace_folder slash 'Simulated_traces_Diff_BT'  '_' Param_name_1{1} '=' num2str(Wall) '_' Param_name_2{1} '=' num2str(delta)];

Params(1) = kel;

Params(2) = P;

Params = [Params, delta, f, Wall];

Gillespie_Func = @(Trace_time) Gillespie_MicroModel_2(Params, t0, Trace_time);

%cmap = colormap(jet(trace_num));

cmap={'r', 'b', 'g', 'k', 'm'};

max_index = floor(Trace_time/t0)+1;

T_tot = (max_index-1)*t0;

t_array = 0:t0:T_tot;

trace_dat = zeros([length(t_array), 6]);

trace_dat(:,1) = t_array';

PL = figure();

for j=1:trace_num

[t_array, x_array, s_array] = feval(Gillespie_Func,Trace_time);

trace_dat(:,j+1) = x_array';

plot(t_array, x_array, '-', 'Color', cmap{j})

hold on

end

hold off

xlabel('Time (s)')

ylabel('Position (bp)')

saveas(PL,[str_trace,'.jpg'], 'jpg')

saveas(PL,[str_trace,'.fig'], 'fig')

col_name = {'Time (s)', 'Trace_1 (bp)', 'Trace_2 (bp)', 'Trace_3 (bp)', 'Trace_4 (bp)', 'Trace_5 (bp)'};

xlswrite([str_trace,'.xlsx'],trace_dat,'Sheet1','A2');

xlswrite([str_trace,'.xlsx'],col_name,'Sheet1','A1');


