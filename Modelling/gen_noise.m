function x_noise = gen_noise(sigma, points)

x_noise = sigma*sqrt(-2*log(rand([1, points]))).*cos(2*pi*rand([1, points]));
