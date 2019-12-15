mainfoldername = 'NewSel_AF_7.5pN_1mM_NTP';

traces = load([mainfoldername,'/' 'Trace_indices.mat'], '-mat');

traces = traces.selected_trace_number;

str_blacklist = [mainfoldername '/' 'blacklist.txt'];

blacklist = load(str_blacklist, '-ascii');

traces = traces(~ismember(traces, blacklist));

subfolder_name = '1mM_indv';

L = length(traces);

kel = zeros([1,L]);

P = zeros([1,L]);

for n=1:L
    
    disp(n)

    selected_trace_indices = traces(n);

    data_name = ['1mM_' num2str(selected_trace_indices)];

    Fit_dat = Analyze_1mM_kel_P_NewSel_SingleTraces_bts(data_name, subfolder_name, selected_trace_indices);

    kel(n) = Fit_dat(1,1);

    P(n) = Fit_dat(1,2);
    
end

corr = corrcoef([kel',P']);

c = corr(1,2);

save([mainfoldername,'/' 'Ptotal.mat'], 'P', '-mat')

save([mainfoldername,'/' 'kel.mat'], 'kel', '-mat')

save([mainfoldername,'/' 'selected_traces.mat'], 'traces', '-mat')

PL1 = figure();

plot(kel,P,'r.','MarkerSize',20)

xlabel('k_{el}')

ylabel('P_{total}')

title(['correlation coef. = ', num2str(c)])

str_fig = [mainfoldername,'/' 'kel-P'];

saveas(PL1,[str_fig,'.fig'], 'fig')

saveas(PL1,[str_fig,'.jpg'], 'jpg')