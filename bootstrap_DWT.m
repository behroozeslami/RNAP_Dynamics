function DWT_BS = bootstrap_DWT(row, DWT_dat, dat_length, nboot)

Size = size(DWT_dat);

col_num = Size(2);

DWT_len = dat_length(row);

DWT_array = DWT_dat(row,1:DWT_len);

zero_mat = zeros([nboot,(col_num-DWT_len)]);

DWT_BS = [cell2mat(arrayfun(@(index) datasample(DWT_array',DWT_len), 1:nboot, 'UniformOutput', false))', zero_mat];