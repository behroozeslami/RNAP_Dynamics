output_foldername = 'Export_AF_7.5pN_1mM_NTP_I';

str_trace_index = [output_foldername '/' 'Trace_indices' '.mat'];

selected_trace_number = load(str_trace_index, '-mat');

TraceIDs = selected_trace_number.selected_trace_number;

str_selected_trace_index = ['Modelling' '/' output_foldername '/' 'Traces_3pause_lowp1' '.txt'];

selected_TraceIDs = load(str_selected_trace_index, '-ascii');

rows=1:length(TraceIDs);

rows = rows(ismember(TraceIDs, selected_TraceIDs));

for j=1:length(rows)
    
    row = rows(j);
    
    [fit_dat_All] = Analyze_fit_pauses_Singletraces(row);
    
end