mainfoldername = 'Export_AF_7.5pN_1mM_NTP';

subfolder = '';

slash = '/';

output_foldername = mainfoldername;

if ~strcmp('',subfolder)
    
    output_foldername = [mainfoldername slash subfolder];
    
end

directory = ['..' slash 'Richard_Data' slash output_foldername];

Batch_name = 'LongPause';
%Batch_name = 'IntermediatePause';
%Batch_name = 'ShortPause';

New_output_foldername = [output_foldername '_' Batch_name];

New_output_path = ['..' slash 'Richard_Data' slash New_output_foldername];

if ~(exist(New_output_path,'dir')==7)
    
    mkdir(New_output_path)
    
end

selected_trace_number = load([output_foldername slash 'TraceID_' Batch_name '.txt'], '-ascii');

L = length(selected_trace_number);

for i=1:L
    
    filename = [directory slash 'RNAP_' num2str(selected_trace_number(i)) '.txt'];
    
    New_filename = [New_output_path slash 'RNAP_' num2str(selected_trace_number(i)) '.txt'];
    
    copyfile(filename, New_filename)
    
end
