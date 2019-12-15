function x_noise = gen_noise(sigma, t0, Trace_time)

max_index = floor(Trace_time/t0) + 1;

x_noise = sigma*sqrt(-2*log(rand([1, max_index]))).*cos(2*pi*rand([1, max_index]));